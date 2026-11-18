"""Step 7.1c (part 2): one day, log1p (7.1) vs three SCTransform v2 variants (7.1c part 1: SCT, SCT_A, SCT_B)
-- UMAPs side by side and metrics.

Same cells, same S - G2M regression, same 30 PCs / 15 neighbours / UMAP settings (seed 1234); only the
normalization differs. No clustering. Gene panels on both embeddings are coloured by the same unregressed log1p
values, so the panels compare the embeddings, not the expression scales.

Metrics (per method):
  pc_covariates.csv   Spearman rho of PC1-10 with log10 counts, log10 genes, %MT, %ribo, nuclear fraction,
                      S - G2M; and R^2 of each covariate (plus genotype 0/1) from a linear fit on the 30 PCs used
                      for the graph (how much of it the embedding carries)
  knn_mixing.csv      per sample: of each cell's neighbours with the same genotype, the share from the other
                      replicate, observed / expected (1 = replicates fully mixed); same-genotype neighbour share
                      observed / expected (> 1 = genotypes separate)
  method_overlap.csv  per SCT variant: variable genes shared with log1p (of 3,000), Jaccard; mean per-cell overlap
                      of its kNN sets with the log1p graph and with the base SCT graph
  hvg_only.csv        top 40 genes variable in base SCT only / log1p only
Writes ${OUT_DIR}/day_analysis/{day}/sct/ (next to part 1's files). Run: via 7.1c_day_sct.sb
"""
import importlib.util
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("day71", HERE / "7.1_day_merged.py")
d71 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d71)
cat_scatter, num_scatter = d71.cat_scatter, d71.num_scatter
SEED, N_PCS, N_NEIGHBORS = d71.SEED, d71.N_PCS, d71.N_NEIGHBORS
SAMPLE_PAL, GENO_PAL, PHASE_PAL = d71.SAMPLE_PAL, d71.GENO_PAL, d71.PHASE_PAL

DAY = sys.argv[1]
RES = Path(os.environ["OUT_DIR"]) / "day_analysis" / DAY
SCT = RES / "sct"
GENES = ["PDX1", "PTF1A", "NKX6-1", "SOX9", "ONECUT1", "ARX", "GCG", "PAX4", "PAX6", "INS", "SLC18A1", "FEV",
         "LMX1A", "SST", "GHRL", "NKX2-2", "NEUROG3", "CHGA", "ISL1", "AFP", "COL1A2", "TAGLN", "MKI67", "TOP2A"]
METHODS = ["log1p", "SCT", "SCT_A", "SCT_B"]
LABEL = {"log1p": "log1p", "SCT": "SCT (S-G2M)", "SCT_A": "SCT A (S-G2M + %MT)", "SCT_B": "SCT B (S-G2M, no MT- HVGs)"}


def neighbour_sets(conn_dist):
    d = conn_dist.tocsr()
    return [d.indices[d.indptr[i]:d.indptr[i + 1]] for i in range(d.shape[0])]


def main():
    a = sc.read_h5ad(RES / f"{DAY}_merged.h5ad")
    a.obsm["X_pca_log1p"] = a.obsm["X_pca"]; a.obsm["X_umap_log1p"] = a.obsm["X_umap"]
    nn = {"log1p": neighbour_sets(a.obsp["distances"])}
    for m in METHODS[1:]:
        pca = pd.read_csv(SCT / f"{m}_pca.csv.gz", index_col=0)
        assert (pca.index == a.obs_names).all()
        a.obsm[f"X_pca_{m}"] = pca.to_numpy()
        sc.pp.neighbors(a, use_rep=f"X_pca_{m}", n_pcs=N_PCS, n_neighbors=N_NEIGHBORS, random_state=SEED, key_added=m)
        sc.tl.umap(a, random_state=SEED, neighbors_key=m)
        a.obsm[f"X_umap_{m}"] = a.obsm["X_umap"].copy()
        nn[m] = neighbour_sets(a.obsp[f"{m}_distances"])

    o = a.obs
    cov = pd.DataFrame({"log10 counts": np.log10(o.total_counts), "log10 genes": np.log10(o.n_genes_by_counts),
                        "% MT": o.pct_counts_mt, "% ribo": o.pct_counts_ribo, "nuclear fraction": o.nuclear_fraction,
                        "S - G2M": o.S_minus_G2M}, index=o.index)
    geno01 = (o.genotype.astype(str) == "KO").astype(float).to_numpy()

    # ---- PC vs covariates
    rows = []
    for m in METHODS:
        X = a.obsm[f"X_pca_{m}"][:, :N_PCS]
        for c in cov:
            r = {"method": m, "covariate": c}
            for k in range(10):
                r[f"PC{k + 1}"] = round(spearmanr(X[:, k], cov[c]).statistic, 2)
            r["R2_30PCs"] = round(LinearRegression().fit(X, cov[c]).score(X, cov[c]), 3)
            rows.append(r)
        rows.append({"method": m, "covariate": "genotype (KO=1)",
                     "R2_30PCs": round(LinearRegression().fit(X, geno01).score(X, geno01), 3),
                     **{f"PC{k + 1}": round(spearmanr(X[:, k], geno01).statistic, 2) for k in range(10)}})
    pcc = pd.DataFrame(rows); pcc.to_csv(SCT / "pc_covariates.csv", index=False)
    print(pcc.to_string(index=False))

    # ---- kNN mixing
    samp = o["sample"].astype(str).to_numpy(); geno = o.genotype.astype(str).to_numpy()
    rep = o.replicate.astype(str).to_numpy()
    rows = []
    for m in METHODS:
        other_rep, same_geno = np.full(a.n_obs, np.nan), np.zeros(a.n_obs)
        for i, nb in enumerate(nn[m]):
            sg = nb[geno[nb] == geno[i]]
            same_geno[i] = len(sg) / len(nb)
            if len(sg):
                other_rep[i] = (rep[sg] != rep[i]).mean()
        for s in d71.SAMPLE_ORDER:
            idx = samp == s; gt = s.split("_")[0]
            n_gt = (geno == gt).sum()
            exp_other = ((geno == gt) & (samp != s)).sum() / (n_gt - 1)
            exp_same = (n_gt - 1) / (a.n_obs - 1)
            rows.append({"method": m, "sample": s, "cells": int(idx.sum()),
                         "other_rep_share": round(np.nanmean(other_rep[idx]), 3), "other_rep_expected": round(exp_other, 3),
                         "replicate_mixing_obs_over_exp": round(np.nanmean(other_rep[idx]) / exp_other, 2),
                         "same_genotype_share": round(same_geno[idx].mean(), 3), "same_genotype_expected": round(exp_same, 3),
                         "genotype_separation_obs_over_exp": round(same_geno[idx].mean() / exp_same, 2)})
    mix = pd.DataFrame(rows); mix.to_csv(SCT / "knn_mixing.csv", index=False)
    print(mix.to_string(index=False))

    # ---- overlap between methods
    def knn_j(m1, m2):
        return round(np.mean([len(np.intersect1d(x, y)) / len(np.union1d(x, y)) for x, y in zip(nn[m1], nn[m2])]), 3)
    h_log = set(a.var_names[a.var["highly_variable"]])
    hv = {m: set(pd.read_csv(SCT / f"{m}_variable_features.csv").gene) for m in METHODS[1:]}
    ov = pd.DataFrame([{"day": DAY, "cells": a.n_obs, "method": m, "hvg_shared_with_log1p": len(h_log & hv[m]),
                        "hvg_jaccard_log1p": round(len(h_log & hv[m]) / len(h_log | hv[m]), 3),
                        "hvg_shared_with_SCT": len(hv["SCT"] & hv[m]),
                        "knn_jaccard_vs_log1p": knn_j("log1p", m), "knn_jaccard_vs_SCT": knn_j("SCT", m)}
                       for m in METHODS[1:]])
    ov.to_csv(SCT / "method_overlap.csv", index=False); print(ov.to_string(index=False))
    vf = pd.read_csv(SCT / "SCT_variable_features.csv"); h_sct = hv["SCT"]
    rank_log = a.var.loc[list(h_log), "dispersions_norm"].sort_values(ascending=False, kind="stable")
    only_sct = vf[~vf.gene.isin(h_log)].sort_values("rank", kind="stable").gene.head(40).tolist()
    only_log = [g for g in rank_log.index if g not in h_sct][:40]
    pd.DataFrame({"SCT_only": pd.Series(only_sct), "log1p_only": pd.Series(only_log)}).to_csv(SCT / "hvg_only.csv", index=False)
    marker_hits = pd.DataFrame({"gene": GENES, "log1p_HVG": [g in h_log for g in GENES],
                                **{f"{m}_HVG": [g in hv[m] for g in GENES] for m in METHODS[1:]}})
    marker_hits.to_csv(SCT / "marker_genes_in_hvg.csv", index=False)

    # ---- figures: rows = methods
    U = {m: a.obsm[f"X_umap_{m}"] for m in METHODS}
    lab = LABEL
    R = len(METHODS)

    fig, axs = plt.subplots(R, 4, figsize=(26, 6.2 * R))
    for r, m in enumerate(METHODS):
        cat_scatter(axs[r, 0], U[m], samp, f"{lab[m]}: sample", SAMPLE_PAL)
        cat_scatter(axs[r, 1], U[m], geno, f"{lab[m]}: genotype", GENO_PAL)
        for c, (gt, gl) in enumerate((("ctr", "Control"), ("KO", "ISL1 KO")), start=2):
            cat_scatter(axs[r, c], U[m], samp, f"{lab[m]}: {gl} (other genotype grey)",
                        {k: v for k, v in SAMPLE_PAL.items() if k.startswith(gt)}, mask=geno == gt)
    fig.suptitle(f"{DAY}: samples and genotype ({a.n_obs:,} cells)"); fig.tight_layout()
    fig.savefig(SCT / "cmp_samples.png", dpi=95); plt.close(fig)

    fig, axs = plt.subplots(R, 5, figsize=(29, 5.5 * R))
    for r, m in enumerate(METHODS):
        for ax, (c, vmin, vmax) in zip(axs[r], [("log10 counts", None, None), ("log10 genes", None, None),
                                               ("% MT", 0, 10), ("% ribo", None, None), ("nuclear fraction", 0.1, 0.6)]):
            num_scatter(ax, U[m], cov[c], f"{lab[m]}: {c}", vmin=vmin, vmax=vmax)
    fig.suptitle(f"{DAY}: QC metrics"); fig.tight_layout(); fig.savefig(SCT / "cmp_qc.png", dpi=95); plt.close(fig)

    def gx(g):
        x = a[:, g].X
        return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.ravel(x)

    fig, axs = plt.subplots(R, 4, figsize=(24, 5.5 * R))
    for r, m in enumerate(METHODS):
        cat_scatter(axs[r, 0], U[m], o.phase, f"{lab[m]}: phase", PHASE_PAL)
        num_scatter(axs[r, 1], U[m], o.S_minus_G2M, f"{lab[m]}: S - G2M", cmap="coolwarm", vmin=-0.5, vmax=0.5)
        num_scatter(axs[r, 2], U[m], gx("MKI67"), f"{lab[m]}: MKI67", cmap="magma")
        num_scatter(axs[r, 3], U[m], gx("TOP2A"), f"{lab[m]}: TOP2A", cmap="magma")
    fig.suptitle(f"{DAY}: cell cycle (S - G2M regressed in both)"); fig.tight_layout()
    fig.savefig(SCT / "cmp_cellcycle.png", dpi=95); plt.close(fig)

    genes = [g for g in GENES if g in a.var_names]
    for m in METHODS:
        ncol = 6; nrow = int(np.ceil(len(genes) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.8 * nrow))
        for ax, g in zip(axs.flat, genes):
            num_scatter(ax, U[m], gx(g), g, cmap="magma")
        for ax in list(axs.flat)[len(genes):]:
            ax.axis("off")
        fig.suptitle(f"{DAY}, {lab[m]} UMAP: log1p expression"); fig.tight_layout()
        fig.savefig(SCT / f"genes_{m}.png", dpi=85); plt.close(fig)

    pv = {"PC": np.arange(1, 31), "log1p": a.uns["pca_variance_ratio"][:30] / a.uns["pca_variance_ratio"][:50].sum()}
    for m in METHODS[1:]:
        pv[m] = pd.read_csv(SCT / f"{m}_pca_variance.csv").variance_ratio[:30].to_numpy()
    pd.DataFrame(pv).round(4).to_csv(SCT / "pca_variance.csv", index=False)
    pd.DataFrame(np.hstack([U[m] for m in METHODS[1:]]), index=a.obs_names,
                 columns=[f"{m}_UMAP{k}" for m in METHODS[1:] for k in (1, 2)]).to_csv(SCT / "sct_umaps.csv.gz")
    print("wrote", SCT)


if __name__ == "__main__":
    main()
