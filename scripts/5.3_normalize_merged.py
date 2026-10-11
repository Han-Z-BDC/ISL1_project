#!/usr/bin/env python
"""
Normalize the merged multi-sample RNA object, both paths -- matches 4.1's
single-sample design, generalized to the merged cross-sample object:
  - log path: scanpy-native normalize_total -> log1p -> HVG -> PCA. No
    scaling, no covariate regression (matches current scanpy convention,
    same reasoning as 4.1).
  - SCT path: Seurat::SCTransform + RunPCA via rpy2 -- R's involvement is
    DELIBERATELY minimized to just this normalization math. No
    IntegrateLayers, no Harmony, nothing else touches R -- the resulting
    PCA embedding is exported back to Python immediately, where
    5.4_harmony_integrate.py runs harmonypy (Python-native, not Seurat's
    HarmonyIntegration wrapper) on both this and the log path's PCA
    identically.

Produces ONE h5ad with obsm["X_pca_log"] and obsm["X_pca_sct"], ready for
5.4_harmony_integrate.py to run harmonypy on either/both.

Optional third path (--sct-persample): SCTransform fit SEPARATELY per
sample rather than once on the pooled object, then Seurat's own
SelectIntegrationFeatures()/PrepSCTIntegration() (not custom logic --
reusing Seurat's tested feature-reconciliation) to get a common gene set
across all 14 per-sample fits, then PCA on the combined residual matrix.
Added after measuring real per-sample divergence in SCTransform's
regularized theta (NB dispersion) on this data -- median per-sample theta
ranged 0.395-1.411 (3.6x) with a typical gene swinging ~2x across the 14
fits (see scripts/_check_sct_per_sample.py) -- confirming the concern
Seurat issue #6054 raises (SCTransform v2's own paper suggests varying
theta by sample) is real here, not hypothetical. NOTE: the OTHER half of
that same concern (a shared log_umi depth-scaling slope getting biased by
pooling) does NOT apply to vst.flavor="v2" specifically -- that slope is
FIXED by v2's own design (confirmed: log_umi coefficient was IDENTICAL,
literally ln(10), across all 14 per-sample fits in that same check), not
data-driven, so pooling can't bias something that isn't being estimated.

Usage:
  python 5.3_normalize_merged.py --merged-rna-h5ad merged_rna_raw.h5ad \
      --config configs/normalize_rna.yaml --outdir out/ [--sct-persample]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

anndata2ri.activate()
anndata2ri.scipy2ri.activate()


def get_cc_genes():
    importr("Seurat")
    s_genes = list(ro.r("Seurat::cc.genes$s.genes"))
    g2m_genes = list(ro.r("Seurat::cc.genes$g2m.genes"))
    return s_genes, g2m_genes


def run_log_path(adata, cfg: dict) -> np.ndarray:
    log_cfg = cfg["log"]
    ad_log = adata.copy()
    sc.pp.normalize_total(ad_log)
    sc.pp.log1p(ad_log)
    sc.pp.highly_variable_genes(ad_log, n_top_genes=log_cfg["n_top_genes"], flavor=log_cfg["hvg_flavor"])
    ad_log_hvg = ad_log[:, ad_log.var["highly_variable"]].copy()
    sc.tl.pca(ad_log_hvg, random_state=cfg["random_state"])
    return ad_log_hvg.obsm["X_pca"]


def run_sct_path(adata, cfg: dict) -> np.ndarray:
    """R's ONLY job here: SCTransform + RunPCA. No clustering, no
    integration, no UMAP -- just the normalization math that has no
    Python equivalent, then immediately back to Python."""
    sct_cfg = cfg["sct"]
    counts = adata.layers.get("soupx_counts", adata.X)
    gene_names = pd.Index(adata.var_names).str.replace("_", "-")

    metadata = adata.obs[["pct_counts_mt", "pct_counts_ribo", "S_score", "G2M_score"]].copy()
    metadata.columns = ["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"]

    ro.globalenv["data_mat"] = counts.T.astype("float32")
    ro.globalenv["gene_names"] = gene_names
    ro.globalenv["cell_names"] = adata.obs_names
    ro.globalenv["metadata"] = metadata
    ro.globalenv["vars_to_regress"] = ro.StrVector(["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"])
    ro.r(f"set.seed({cfg['random_state']})")

    importr("Seurat")
    print("  [run_sct_path] building Seurat object...", flush=True)
    ro.r("""
        library(Matrix)
        counts <- as(data_mat, "dgCMatrix")
        rownames(counts) <- gene_names
        colnames(counts) <- cell_names
        seu <- CreateSeuratObject(counts = counts, meta.data = metadata)
    """)
    print("  [run_sct_path] running SCTransform...", flush=True)
    ro.r(f"""
        seu <- SCTransform(seu, vst.flavor = "{sct_cfg['vst_flavor']}", assay = "RNA",
                            vars.to.regress = vars_to_regress)
    """)
    print("  [run_sct_path] RunPCA...", flush=True)
    ro.r(f"""
        DefaultAssay(seu) <- "SCT"
        seu <- RunPCA(seu, npcs = {sct_cfg['n_pcs']}, reduction.name = "sct.pca")
    """)

    r_cell_order = list(ro.r("Cells(seu)"))
    if r_cell_order != list(adata.obs_names):
        raise ValueError(
            f"Cell order changed inside R during SCTransform/RunPCA -- "
            f"positional transfer back to adata would misalign cells. "
            f"({len(r_cell_order)} R cells vs {adata.n_obs} in adata)"
        )
    pca_emb = np.asarray(ro.r('Embeddings(seu, reduction = "sct.pca")'))
    print("  [run_sct_path] done, R's job here is finished.", flush=True)
    return pca_emb


def run_sct_persample_path(adata, cfg: dict, n_integration_features: int = 3000) -> tuple:
    """SCTransform fit separately per sample_id, then Seurat's own
    SelectIntegrationFeatures()/PrepSCTIntegration() (not reinvented here)
    to get one common gene set with recomputed residuals across all 14
    fits, then PCA on the cbind'd residual matrix via irlba (same
    algorithm Seurat's own RunPCA uses under the hood).

    Deliberately does NOT use Seurat's merge()/IntegrateLayers on already-
    SCTransform'd objects -- scale.data (residuals) handling across merge()
    on SCT assays is not a first-class supported path (that's precisely
    why PrepSCTIntegration + anchor-based IntegrateData/IntegrateLayers
    exist). Extracting each object's own scale.data for the common
    anchor.features gene set and cbind'ing manually sidesteps that
    entirely -- explicit and verifiable rather than relying on merge()
    behavior this pipeline hasn't tested.

    Returns (pca_embedding, cell_names_in_pca_row_order) -- caller must
    align back to adata.obs_names explicitly, do not assume order.
    """
    sct_cfg = cfg["sct"]
    sample_ids = adata.obs["sample_id"].astype(str).unique().tolist()
    print(f"  [run_sct_persample_path] {len(sample_ids)} samples: {sample_ids}", flush=True)

    importr("Seurat")
    ro.r(f"set.seed({cfg['random_state']})")
    ro.r("library(Matrix)")
    ro.r("object_list <- list()")

    for sid in sample_ids:
        sub = adata[adata.obs["sample_id"].astype(str) == sid]
        counts = sub.layers.get("soupx_counts", sub.X)
        gene_names = pd.Index(sub.var_names).str.replace("_", "-")
        metadata = sub.obs[["pct_counts_mt", "pct_counts_ribo", "S_score", "G2M_score"]].copy()
        metadata.columns = ["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"]

        ro.globalenv["data_mat"] = counts.T.astype("float32")
        ro.globalenv["gene_names"] = gene_names
        ro.globalenv["cell_names"] = sub.obs_names
        ro.globalenv["metadata"] = metadata
        ro.globalenv["vars_to_regress"] = ro.StrVector(["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"])
        ro.globalenv["sid"] = sid

        print(f"  [run_sct_persample_path] [{sid}] n_cells={sub.n_obs} -- SCTransform...", flush=True)
        ro.r(f"""
            counts <- as(data_mat, "dgCMatrix")
            rownames(counts) <- gene_names
            colnames(counts) <- cell_names
            seu <- CreateSeuratObject(counts = counts, meta.data = metadata)
            seu <- SCTransform(seu, vst.flavor = "{sct_cfg['vst_flavor']}", assay = "RNA",
                                vars.to.regress = vars_to_regress, verbose = FALSE)
            DefaultAssay(seu) <- "SCT"
            object_list[[sid]] <- seu
            rm(seu); gc()
        """)

    print(f"  [run_sct_persample_path] SelectIntegrationFeatures (nfeatures={n_integration_features})...", flush=True)
    ro.r(f"""
        features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = {n_integration_features})
        object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = features)
    """)

    print("  [run_sct_persample_path] combining per-sample residuals for the common feature set...", flush=True)
    ro.r("""
        residual_mats <- lapply(object_list, function(o) GetAssayData(o, assay = "SCT", layer = "scale.data")[features, , drop = FALSE])
        combined_residuals <- do.call(cbind, residual_mats)
    """)
    combined_cell_names = list(ro.r("colnames(combined_residuals)"))
    n_features_used = int(ro.r("length(features)")[0])
    print(f"  [run_sct_persample_path] combined matrix: {n_features_used} common features x {len(combined_cell_names)} cells", flush=True)

    print(f"  [run_sct_persample_path] PCA (irlba, {sct_cfg['n_pcs']} PCs)...", flush=True)
    ro.r(f"""
        library(irlba)
        pca <- prcomp_irlba(t(combined_residuals), n = {sct_cfg['n_pcs']}, center = TRUE, scale. = FALSE)
    """)
    pca_emb = np.asarray(ro.r("pca$x"))
    print("  [run_sct_persample_path] done.", flush=True)
    return pca_emb, combined_cell_names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-rna-h5ad", required=True, help="5.2_merge_rna.py output: merged_rna_raw.h5ad")
    ap.add_argument("--config", required=True, help="Path to configs/normalize_rna.yaml (reused from step 4.1)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sct-persample", action="store_true",
                     help="Also fit SCTransform separately per sample (Seurat's own "
                          "SelectIntegrationFeatures/PrepSCTIntegration to reconcile features) "
                          "instead of only the pooled fit -- adds obsm['X_pca_sct_persample']. "
                          "Real per-sample divergence in SCTransform's theta was measured on "
                          "this data before adding this (see the module docstring).")
    ap.add_argument("--n-integration-features", type=int, default=3000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    adata = sc.read_h5ad(args.merged_rna_h5ad)
    print(f"loaded merged object: {adata.shape}", flush=True)

    s_genes, g2m_genes = get_cc_genes()
    s_genes = [g for g in s_genes if g in adata.var_names]
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]
    ad_for_cc = adata.copy()
    sc.pp.normalize_total(ad_for_cc)
    sc.pp.log1p(ad_for_cc)
    sc.tl.score_genes_cell_cycle(ad_for_cc, s_genes=s_genes, g2m_genes=g2m_genes)
    adata.obs["S_score"] = ad_for_cc.obs["S_score"].to_numpy()
    adata.obs["G2M_score"] = ad_for_cc.obs["G2M_score"].to_numpy()
    adata.obs["Phase"] = ad_for_cc.obs["phase"].to_numpy()
    print("cell-cycle scores computed", flush=True)

    print("starting log path (pure scanpy)...", flush=True)
    adata.obsm["X_pca_log"] = run_log_path(adata, cfg)
    print("log path done. starting SCT path (R, minimal)...", flush=True)
    adata.obsm["X_pca_sct"] = run_sct_path(adata, cfg)
    print("SCT path done.", flush=True)

    if args.sct_persample:
        print("starting SCT per-sample path (R, per-sample fits + Seurat's own feature reconciliation)...", flush=True)
        pca_emb, cell_order = run_sct_persample_path(adata, cfg, args.n_integration_features)
        if set(cell_order) != set(adata.obs_names):
            raise ValueError(
                f"run_sct_persample_path returned a different cell set than adata "
                f"({len(cell_order)} vs {adata.n_obs}) -- cannot align."
            )
        # combined_residuals' column order is sample-then-cell, not necessarily
        # adata.obs_names' order -- reindex explicitly rather than assume.
        pos = {c: i for i, c in enumerate(cell_order)}
        reorder = np.array([pos[c] for c in adata.obs_names])
        adata.obsm["X_pca_sct_persample"] = pca_emb[reorder]
        print("SCT per-sample path done.", flush=True)

    out_path = outdir / "merged_rna_normalized.h5ad"
    adata.write_h5ad(out_path)
    msg = f"Written: {out_path} -- obsm has X_pca_log {adata.obsm['X_pca_log'].shape}, X_pca_sct {adata.obsm['X_pca_sct'].shape}"
    if args.sct_persample:
        msg += f", X_pca_sct_persample {adata.obsm['X_pca_sct_persample'].shape}"
    print(msg)


if __name__ == "__main__":
    main()
