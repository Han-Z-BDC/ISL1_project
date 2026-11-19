"""Step 7.1g: split one 7.1b cluster on the day's SCT B graph (Leiden restricted to that cluster's cells).

Purpose (user, 2026-11-18): D22 r1.0 cluster 1 has two arms on the UMAP, an SST / HHEX delta-like arm and a
GHRL / ARX / NEUROG3 stem; split it without a gene threshold. The day's graph and UMAP are not refitted: Leiden
(igraph, 2 iterations, seed 1234, as in 7.1b) runs on the stored SCT_B graph restricted to the cluster's cells,
at SUB_RES resolutions. Gene values are the unregressed log1p data (X).

Writes ${OUT_DIR}/day_analysis/{day}/clustering/split_r{res}_c{cluster}/ (reruns overwrite):
  subclusters.csv         barcode -> subcluster at each sub-resolution (cells of the cluster only)
  gene_by_subcluster.csv  per sub-resolution x subcluster x gene: cells, ctr / KO cells, pct > 0, pct > 1, mean log1p
  markers_s{sr}.csv       top 30 Wilcoxon markers of each subcluster vs the rest of the cluster
  umap_split.png          zoom on the cluster, coloured by subcluster at each sub-resolution
  feature_genes.png       zoom on the cluster, log1p feature plots of the genes
Run: sbatch analysis/scripts/7.1g_split_cluster.sb D22 1.0 1 SST HHEX ...
"""
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

OUT_DIR = Path(os.environ["OUT_DIR"])
DAY, RES, CL, GENES = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
KEY = "leiden_r" + RES.replace(".", "p")
OUT = OUT_DIR / "day_analysis" / DAY / "clustering" / f"split_r{RES}_c{CL}"
SUB_RES = [0.1, 0.2, 0.3, 0.5]
SEED = 1234
GREY = "#e3e3e3"
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", ["#dcdcdc", "#fdbb84", "#e34a33", "#7f0000"])


def zoom(ax, u, m):
    lo, hi = u[m].min(axis=0), u[m].max(axis=0)
    pad = 0.15 * (hi - lo)
    ax.set_xlim(lo[0] - pad[0], hi[0] + pad[0]); ax.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
    ax.set_xticks([]); ax.set_yticks([])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in OUT.glob("*"):
        p.unlink()
    a = sc.read_h5ad(OUT_DIR / "day_analysis" / DAY / f"{DAY}_merged.h5ad")
    missing = [g for g in GENES if g not in a.var_names]
    if missing:
        sys.exit(f"genes not found: {missing}")
    a.obs[KEY] = a.obs[KEY].astype(str).astype("category")
    m = (a.obs[KEY] == CL).to_numpy()
    ko = a.obs["genotype"].astype(str).str.upper().str.startswith("KO").to_numpy()
    u = np.asarray(a.obsm["X_umap_SCT_B"])
    print(DAY, KEY, "cluster", CL, ":", m.sum(), "cells")

    subs = {}
    for sr in SUB_RES:
        k = f"sub_s{sr}"
        sc.tl.leiden(a, resolution=sr, key_added=k, restrict_to=(KEY, [CL]), neighbors_key="SCT_B",
                     random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
        subs[sr] = a.obs.loc[m, k].astype(str).str.split(",").str[1]
        print(f"s{sr}:", subs[sr].value_counts().sort_index().to_dict())
    pd.DataFrame({f"s{sr}": s for sr, s in subs.items()}).to_csv(OUT / "subclusters.csv")

    b = a[m].copy()
    x = b[:, GENES].X
    x = pd.DataFrame(np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x), index=b.obs_names, columns=GENES)
    kob = ko[m]
    rows = []
    for sr, s in subs.items():
        for c in sorted(s.unique(), key=int):
            mm = (s == c).to_numpy()
            for g in GENES:
                v = x.loc[mm, g]
                rows.append({"sub_res": sr, "subcluster": c, "gene": g, "cells": int(mm.sum()),
                             "ctr_cells": int((mm & ~kob).sum()), "KO_cells": int((mm & kob).sum()),
                             "pct_gt0": round(100 * (v > 0).mean(), 1), "pct_gt1": round(100 * (v > 1).mean(), 1),
                             "mean_log1p": round(v.mean(), 3)})
        b.obs["s"] = s.astype("category")
        if b.obs["s"].nunique() > 1:
            sc.tl.rank_genes_groups(b, "s", method="wilcoxon", use_raw=False)
            mk = sc.get.rank_genes_groups_df(b, None)
            mk = mk.groupby("group", observed=True, group_keys=False).head(30)
            mk.to_csv(OUT / f"markers_s{sr}.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "gene_by_subcluster.csv", index=False)

    cmap = plt.get_cmap("tab10")
    fig, axs = plt.subplots(1, len(SUB_RES), figsize=(5.5 * len(SUB_RES), 5.2), squeeze=False)
    for ax, (sr, s) in zip(axs.flat, subs.items()):
        ax.scatter(u[~m, 0], u[~m, 1], c=GREY, s=4, linewidths=0, rasterized=True)
        um = u[m]
        for i, c in enumerate(sorted(s.unique(), key=int)):
            mm = (s == c).to_numpy()
            ax.scatter(um[mm, 0], um[mm, 1], color=cmap(i % 10), s=10, linewidths=0, label=f"{c} ({mm.sum()})")
            xy = np.median(um[mm], axis=0)
            ax.text(xy[0], xy[1], c, fontsize=11, weight="bold", ha="center", va="center")
        ax.legend(fontsize=8, markerscale=1.5, loc="best")
        ax.set_title(f"{DAY} r{RES} cluster {CL}: sub-resolution {sr}"); zoom(ax, u, m)
    fig.tight_layout(); fig.savefig(OUT / "umap_split.png", dpi=100); plt.close(fig)

    ncol = 4
    nrow = int(np.ceil(len(GENES) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.7 * nrow), squeeze=False)
    xa = a[:, GENES].X
    xa = np.asarray(xa.todense()) if hasattr(xa, "todense") else np.asarray(xa)
    for j, (ax, g) in enumerate(zip(axs.flat, GENES)):
        v = xa[:, j]
        vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1
        o = np.argsort(v, kind="stable")
        s_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=8, cmap=FEAT_CMAP, vmin=0, vmax=vmax, linewidths=0, rasterized=True)
        plt.colorbar(s_, ax=ax, fraction=0.04, pad=0.01)
        ax.set_title(g); zoom(ax, u, m)
    for ax in list(axs.flat)[len(GENES):]:
        ax.axis("off")
    fig.suptitle(f"{DAY}: log1p, zoom on r{RES} cluster {CL}"); fig.tight_layout()
    fig.savefig(OUT / "feature_genes.png", dpi=90); plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
