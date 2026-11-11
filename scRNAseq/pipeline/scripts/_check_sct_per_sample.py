#!/usr/bin/env python
"""
Diagnostic (not a pipeline step): fit SCTransform separately on each of the
14 samples and compare the regularized per-gene model parameters (theta,
log_umi coefficient) across samples. If they're tightly clustered, pooling
for 5.3's SCT path is fine. If they diverge meaningfully, that's real
evidence for a per-sample-then-merge rewrite.

Mirrors 5.3_normalize_merged.py's run_sct_path() exactly (same
CreateSeuratObject/SCTransform call, same vars.to.regress, same cell-cycle
scoring), just looped per-sample instead of pooled, and extracting
SCTResults() feature.attrs instead of running RunPCA.
"""
import sys
import numpy as np
import pandas as pd
import scanpy as sc

import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

anndata2ri.activate()
anndata2ri.scipy2ri.activate()

merged_h5ad = sys.argv[1]

importr("Seurat")


def get_cc_genes():
    s_genes = list(ro.r("Seurat::cc.genes$s.genes"))
    g2m_genes = list(ro.r("Seurat::cc.genes$g2m.genes"))
    return s_genes, g2m_genes


adata = sc.read_h5ad(merged_h5ad)
print(f"loaded: {adata.shape}", flush=True)
print(f"obs columns: {list(adata.obs.columns)}", flush=True)

s_genes, g2m_genes = get_cc_genes()
s_genes = [g for g in s_genes if g in adata.var_names]
g2m_genes = [g for g in g2m_genes if g in adata.var_names]
ad_cc = adata.copy()
sc.pp.normalize_total(ad_cc)
sc.pp.log1p(ad_cc)
sc.tl.score_genes_cell_cycle(ad_cc, s_genes=s_genes, g2m_genes=g2m_genes)
adata.obs["S_score"] = ad_cc.obs["S_score"].to_numpy()
adata.obs["G2M_score"] = ad_cc.obs["G2M_score"].to_numpy()
del ad_cc  # full pooled-size copy, only needed transiently -- don't keep it resident for all 14 fits
print("cell-cycle scores computed", flush=True)

sample_ids = adata.obs["sample_id"].astype(str).unique().tolist()
print(f"{len(sample_ids)} samples: {sample_ids}", flush=True)

ro.r('library(Matrix)')
summary_rows = []
per_gene_tables = {}

for sid in sample_ids:
    sub = adata[adata.obs["sample_id"].astype(str) == sid].copy()
    n_cells = sub.n_obs
    counts = sub.layers.get("soupx_counts", sub.X)
    gene_names = pd.Index(sub.var_names).str.replace("_", "-")
    metadata = sub.obs[["pct_counts_mt", "pct_counts_ribo", "S_score", "G2M_score"]].copy()
    metadata.columns = ["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"]

    ro.globalenv["data_mat"] = counts.T.astype("float32")
    ro.globalenv["gene_names"] = gene_names
    ro.globalenv["cell_names"] = sub.obs_names
    ro.globalenv["metadata"] = metadata
    ro.globalenv["vars_to_regress"] = ro.StrVector(["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"])

    print(f"[{sid}] n_cells={n_cells} -- building + SCTransform...", flush=True)
    ro.r("""
        counts <- as(data_mat, "dgCMatrix")
        rownames(counts) <- gene_names
        colnames(counts) <- cell_names
        seu <- CreateSeuratObject(counts = counts, meta.data = metadata)
        seu <- SCTransform(seu, vst.flavor = "v2", assay = "RNA",
                            vars.to.regress = vars_to_regress, verbose = FALSE)
    """)
    fa = ro.r('SCTResults(seu, assay = "SCT", slot = "feature.attributes")')
    fa_df = ro.conversion.get_conversion().rpy2py(fa)
    if sid == sample_ids[0]:
        print(f"  feature.attrs columns: {list(fa_df.columns)}", flush=True)
    per_gene_tables[sid] = fa_df
    theta_col = [c for c in fa_df.columns if "theta" in c.lower()]
    umi_col = [c for c in fa_df.columns if "log_umi" in c.lower() or "log10_umi" in c.lower()]
    row = {"sample_id": sid, "n_cells": n_cells, "n_genes_fit": fa_df.shape[0]}
    if theta_col:
        row["median_theta"] = float(np.median(fa_df[theta_col[0]]))
    if umi_col:
        row["median_log_umi_coef"] = float(np.median(fa_df[umi_col[0]]))
    summary_rows.append(row)
    print(f"  {row}", flush=True)
    ro.r("rm(seu); gc()")

summary = pd.DataFrame(summary_rows)
print("\n=== per-sample summary ===", flush=True)
print(summary.to_string(index=False), flush=True)

# genes common to every sample's fit, for an apples-to-apples comparison
common_genes = set.intersection(*[set(t.index) for t in per_gene_tables.values()])
print(f"\n{len(common_genes)} genes fit in every sample", flush=True)
common_genes = sorted(common_genes)

theta_col_name = [c for c in next(iter(per_gene_tables.values())).columns if "theta" in c.lower()][0]
umi_col_name = [c for c in next(iter(per_gene_tables.values())).columns if "log_umi" in c.lower() or "log10_umi" in c.lower()][0]

theta_mat = pd.DataFrame({sid: t.loc[common_genes, theta_col_name] for sid, t in per_gene_tables.items()})
umi_mat = pd.DataFrame({sid: t.loc[common_genes, umi_col_name] for sid, t in per_gene_tables.items()})

theta_mat.to_csv("sct_per_sample_theta.csv")
umi_mat.to_csv("sct_per_sample_log_umi_coef.csv")

print("\n=== spread across samples, for genes fit everywhere (log1p(theta) since theta is highly right-skewed) ===", flush=True)
log_theta = np.log1p(theta_mat)
per_gene_range_logtheta = (log_theta.max(axis=1) - log_theta.min(axis=1))
per_gene_range_umi = (umi_mat.max(axis=1) - umi_mat.min(axis=1))
print("median per-gene range (log1p theta, across 14 samples):", per_gene_range_logtheta.median(), flush=True)
print("median per-gene range (log_umi coef, across 14 samples):", per_gene_range_umi.median(), flush=True)
print("median per-gene log_umi coef itself (typical magnitude, for scale):", umi_mat.median(axis=1).median(), flush=True)

print("\nWritten: sct_per_sample_theta.csv, sct_per_sample_log_umi_coef.csv")
