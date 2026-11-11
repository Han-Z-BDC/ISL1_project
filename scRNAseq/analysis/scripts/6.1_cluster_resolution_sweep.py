"""Step 6.1: Leiden resolution sweep on the chosen embedding (log, unintegrated).

Embedding chosen by the user 2026-10-11 from 5.4b: log path, no integration (all 4 samples
share one OCM channel/library; replicates already mix). Kept for the 12-sample run
(2026-11-08: D15, D22, D40 x ctr/KO x rep1/rep2), where integrating on diff_batch would also
remove day differences (Diff37 = D15 + D22, Diff34 = D40). Neighbors/Leiden conventions match
the lab pipeline's RNA track (pipeline/configs/normalize_rna.yaml: 30 PCs, n_neighbors=15,
random_state=1234) and diff-multiome's _cluster_resolution_sweep.py (flavor="igraph",
n_iterations=2, directed=False, same resolution list). The UMAP is recomputed from this same
graph so clusters and UMAP agree.

Inputs:  ${OUT_DIR}/integration/merged_rna_normalized.h5ad (pipeline 5.3)
Outputs: ${OUT_DIR}/clustering/
  h5ad/clustered_log_unintegrated.h5ad   all cells; layers counts/soupx_counts/log1p_norm,
                                         obsm X_pca_log + X_umap, obs leiden_r* for every resolution
  cluster_sweep.csv                      barcode -> cluster at every resolution
  n_clusters.csv
  composition_r*.csv                     cluster x sample cell counts, % of each day-genotype's cells
  markers_r*.csv                         top 10 Wilcoxon markers per cluster (log1p_norm)
  lineage.csv                            each cluster's main parent at the next-lower resolution
  plots/umap_leiden_grid.png, umap_genotype.png (genotype / replicate / day), dotplot_r*.png,
        composition_r*.png

Run through Slurm: sbatch analysis/scripts/6.1_cluster_resolution_sweep.sb
"""
import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0]
N_PCS = 30
N_NEIGHBORS = 15
RANDOM_STATE = 1234

# grouped for the dot plots: label -> genes
MARKERS = {
    "beta": ["INS", "IAPP", "NKX6-1", "MAFA", "G6PC2"],
    "alpha": ["GCG", "ARX", "IRX2", "TTR"],
    "delta": ["SST", "HHEX"],
    "PP/eps": ["PPY", "GHRL"],
    "EC": ["TPH1", "FEV", "SLC18A1", "LMX1A"],
    "pan-endo": ["CHGA", "ISL1", "NEUROD1", "NKX2-2"],
    "endo prog": ["NEUROG3"],
    "duct/prog": ["SOX9", "KRT19", "CFTR", "PDX1"],
    "mesench": ["VIM", "COL1A1", "COL3A1"],
    "cycling": ["MKI67", "TOP2A"],
    "gut": ["CDX2"],
}


def key_for(res):
    return f"leiden_r{str(res).replace('.', 'p')}"


def label_umap(ax, u, cats, title):
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(cats.categories):
        m = cats == c
        ax.scatter(u[m, 0], u[m, 1], s=1.2, color=cmap(i % 20), linewidths=0)
    for c in cats.categories:
        m = cats == c
        x, y = np.median(u[m, 0]), np.median(u[m, 1])
        ax.text(x, y, c, fontsize=9, weight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    ax.set_title(title, fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--normalized-h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    for d in ("h5ad", "plots"):
        (out / d).mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(args.normalized_h5ad)
    print(f"loaded {a.shape}", flush=True)
    lay = "soupx_counts" if "soupx_counts" in a.layers else "counts"
    a.X = a.layers[lay].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.layers["log1p_norm"] = a.X.copy()
    a.obs["day"] = a.obs["stage"].astype(str).str.split("-").str[0]       # D15 / D22 / D40
    a.obs["genotype"] = a.obs["stage"].astype(str).str.split("-").str[1]  # ctr / KO

    a.obsm["X_pca_log30"] = a.obsm["X_pca_log"][:, :N_PCS]
    sc.pp.neighbors(a, use_rep="X_pca_log30", n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)
    sc.tl.umap(a, random_state=RANDOM_STATE)
    u = a.obsm["X_umap"]

    n_cl = {}
    for res in RESOLUTIONS:
        k = key_for(res)
        sc.tl.leiden(a, resolution=res, key_added=k, random_state=RANDOM_STATE,
                     flavor="igraph", n_iterations=2, directed=False)
        n_cl[res] = a.obs[k].nunique()
        print(f"  r={res}: {n_cl[res]} clusters", flush=True)
    keys = [key_for(r) for r in RESOLUTIONS]
    pd.DataFrame({"barcode": a.obs_names, **{k: a.obs[k].astype(str).to_numpy() for k in keys}}).to_csv(
        out / "cluster_sweep.csv", index=False)
    pd.DataFrame({"resolution": list(n_cl), "n_clusters": list(n_cl.values())}).to_csv(out / "n_clusters.csv", index=False)

    # UMAP grid
    ncols = 4
    fig, axes = plt.subplots(2, ncols, figsize=(6.2 * ncols, 11))
    for ax, res in zip(axes.ravel(), RESOLUTIONS):
        label_umap(ax, u, pd.Categorical(a.obs[key_for(res)]), f"resolution {res}: {n_cl[res]} clusters")
    fig.tight_layout()
    fig.savefig(out / "plots" / "umap_leiden_grid.png", dpi=85)
    plt.close(fig)

    # genotype / replicate reference on the same UMAP
    pal = {"ctr": "#3b6fb6", "KO": "#d0573f", "rep1": "#4c9a6a", "rep2": "#9b6bc0",
           "D15": "#e3a72f", "D22": "#2a9d8f", "D40": "#7a4fa3"}
    days = sorted(a.obs["day"].unique(), key=lambda d: int(d[1:]))
    order = np.random.default_rng(0).permutation(a.n_obs)
    fig, axes = plt.subplots(1, 5, figsize=(27.5, 5.2), sharex=True, sharey=True)
    g = a.obs["genotype"].to_numpy(); r = a.obs["tech_batch"].astype(str).to_numpy()
    axes[0].scatter(u[order, 0], u[order, 1], s=1.2, c=[pal[x] for x in g[order]], linewidths=0)
    axes[0].set_title("merged (ctr + KO)")
    for ax, gg in zip(axes[1:3], ("ctr", "KO")):
        m = g == gg
        ax.scatter(u[~m, 0], u[~m, 1], s=1.2, c="#d9d9d9", linewidths=0)
        ax.scatter(u[m, 0], u[m, 1], s=1.2, c=pal[gg], linewidths=0)
        ax.set_title(f"{gg} only (n={m.sum():,})")
    axes[3].scatter(u[order, 0], u[order, 1], s=1.2, c=[pal[x] for x in r[order]], linewidths=0)
    axes[3].set_title("replicate")
    d = a.obs["day"].to_numpy()
    axes[4].scatter(u[order, 0], u[order, 1], s=1.2, c=[pal[x] for x in d[order]], linewidths=0)
    axes[4].set_title("day")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
    for k, ax in [("ctr", axes[0]), ("KO", axes[0]), ("rep1", axes[3]), ("rep2", axes[3])] + [(x, axes[4]) for x in days]:
        ax.scatter([], [], c=pal[k], s=20, label=k)
    axes[0].legend(fontsize=9); axes[3].legend(fontsize=9); axes[4].legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "plots" / "umap_genotype.png", dpi=85)
    plt.close(fig)

    markers = {grp: [x for x in gs if x in a.var_names] for grp, gs in MARKERS.items()}
    markers = {grp: gs for grp, gs in markers.items() if gs}
    # day -> genotype -> replicate; groups = day-genotype (2 replicates each)
    groups = [f"{dd}-{gg}" for dd in days for gg in ("ctr", "KO")]
    sample_order = sorted(a.obs["sample_id"].astype(str).unique(),
                          key=lambda x: (groups.index(x.split("_")[1]), x.split("_")[2]))
    n_group = a.obs["stage"].astype(str).value_counts()
    lineage = []
    for i, res in enumerate(RESOLUTIONS):
        k = key_for(res)
        # dot plot of canonical markers
        dp = sc.pl.dotplot(a, markers, groupby=k, layer="log1p_norm", standard_scale="var",
                           return_fig=True, show=False)
        dp.savefig(out / "plots" / f"dotplot_{k}.png", dpi=85, bbox_inches="tight")
        plt.close("all")

        # composition: counts per sample, then each cluster's share of every day-genotype's cells
        # and, per day, the KO share normalised to genotype size
        comp = pd.crosstab(a.obs[k], a.obs["sample_id"]).reindex(columns=sample_order, fill_value=0)
        for grp in groups:
            cols = [c for c in sample_order if c.split("_")[1] == grp]
            comp[grp] = comp[cols].sum(1)
        comp["n"] = comp[groups].sum(1)
        for grp in groups:
            comp[f"pct_of_{grp}_cells"] = (100 * comp[grp] / n_group[grp]).round(1)
        for dd in days:
            pc, pk_ = comp[f"pct_of_{dd}-ctr_cells"], comp[f"pct_of_{dd}-KO_cells"]
            comp[f"{dd}_KO_share_norm"] = (pk_ / (pc + pk_)).round(2)
        comp.to_csv(out / f"composition_{k}.csv")

        fig, axs = plt.subplots(len(days), 1, figsize=(max(6, 0.45 * len(comp)), 2.6 * len(days)), sharex=True)
        x = np.arange(len(comp))
        for ax, dd in zip(np.atleast_1d(axs), days):
            ax.bar(x - 0.2, comp[f"pct_of_{dd}-ctr_cells"], 0.4, color=pal["ctr"], label=f"% of {dd} ctr cells")
            ax.bar(x + 0.2, comp[f"pct_of_{dd}-KO_cells"], 0.4, color=pal["KO"], label=f"% of {dd} KO cells")
            ax.set_ylabel(f"{dd}: % of cells"); ax.legend(fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels(comp.index, fontsize=8)
        ax.set_xlabel(f"cluster (resolution {res})")
        fig.tight_layout()
        fig.savefig(out / "plots" / f"composition_{k}.png", dpi=85)
        plt.close(fig)

        # top markers per cluster
        sc.tl.rank_genes_groups(a, groupby=k, method="wilcoxon", layer="log1p_norm", use_raw=False,
                                pts=True, key_added=f"rgg_{k}")
        df = sc.get.rank_genes_groups_df(a, group=None, key=f"rgg_{k}")
        df = df[df["logfoldchanges"] > 0].groupby("group", observed=True).head(10)
        df.to_csv(out / f"markers_{k}.csv", index=False)
        del a.uns[f"rgg_{k}"]

        if i > 0:
            pk = key_for(RESOLUTIONS[i - 1])
            ct = pd.crosstab(a.obs[k], a.obs[pk])
            for c, row in ct.iterrows():
                lineage.append({"resolution": res, "cluster": c, "n": int(row.sum()),
                                "parent_resolution": RESOLUTIONS[i - 1], "main_parent": row.idxmax(),
                                "frac_from_main_parent": round(row.max() / row.sum(), 2)})
        print(f"  r={res}: dotplot, composition, markers done", flush=True)
    pd.DataFrame(lineage).to_csv(out / "lineage.csv", index=False)

    a.X = a.layers["soupx_counts"] if "soupx_counts" in a.layers else a.layers["counts"]
    a.write_h5ad(out / "h5ad" / "clustered_log_unintegrated.h5ad")
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
