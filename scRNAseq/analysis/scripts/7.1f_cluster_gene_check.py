"""Step 7.1f: where chosen genes sit across the clusters of one 7.1b resolution (after 7.1b).

Purpose (user, 2026-11-18): on D22, find the delta-like bridge (SST, HHEX) and the NEUROG3+ endocrine progenitors
at resolution 1.0 before deciding how to label them on top of the r0.2 labels.
Values are the unregressed log1p data (X), as in 7.1b. "pct > 0" = share of the cluster's cells with any
expression; "pct > 1" = share with log1p > 1 (clearly expressing).
Dominant hormone: for cells expressing any of INS, GCG, SST, GHRL, the hormone with the highest log1p value.

Writes ${OUT_DIR}/day_analysis/{day}/clustering/check_r{res}/ (reruns overwrite):
  gene_by_cluster.csv     per cluster x gene: cells, ctr / KO cells, pct > 0, pct > 1, mean log1p
  hormone_by_cluster.csv  per cluster: cells by dominant hormone (none = no hormone detected)
  umap_clusters.png       UMAP with cluster numbers (all cells, control, KO)
  feature_genes.png       UMAP feature plots of the genes, with cluster numbers drawn on
Run: sbatch analysis/scripts/7.1f_cluster_gene_check.sb D22 1.0 SST HHEX ...
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

OUT_DIR = Path(os.environ["OUT_DIR"])
DAY, RES, GENES = sys.argv[1], sys.argv[2], sys.argv[3:]
KEY = "leiden_r" + RES.replace(".", "p")
OUT = OUT_DIR / "day_analysis" / DAY / "clustering" / f"check_r{RES}"
HORMONES = ["INS", "GCG", "SST", "GHRL"]
SEED = 1234
GREY = "#d9d9d9"
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", ["#dcdcdc", "#fdbb84", "#e34a33", "#7f0000"])


def label_clusters(ax, u, cl):
    for c in sorted(cl.unique(), key=int):
        xy = np.median(u[(cl == c).to_numpy()], axis=0)
        ax.text(xy[0], xy[1], c, fontsize=9, weight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in OUT.glob("*"):
        p.unlink()
    a = ad.read_h5ad(OUT_DIR / "day_analysis" / DAY / f"{DAY}_merged.h5ad", backed="r")
    genes = list(dict.fromkeys(GENES + HORMONES))
    missing = [g for g in genes if g not in a.var_names]
    if missing:
        sys.exit(f"genes not found: {missing}")
    obs = a.obs[[KEY, "genotype"]].copy()
    u = np.asarray(a.obsm["X_umap_SCT_B"])
    x = a[:, genes].to_memory().X
    x = pd.DataFrame(np.asarray(x.todense()) if hasattr(x, "todense") else np.asarray(x),
                     index=obs.index, columns=genes)
    cl = obs[KEY].astype(str)
    cats = sorted(cl.unique(), key=int)
    ko = obs["genotype"].astype(str).str.upper().str.startswith("KO")

    rows = []
    for c in cats:
        m = (cl == c).to_numpy()
        for g in GENES:
            v = x.loc[m, g]
            rows.append({"cluster": c, "gene": g, "cells": int(m.sum()), "ctr_cells": int((m & ~ko).sum()),
                         "KO_cells": int((m & ko).sum()), "pct_gt0": round(100 * (v > 0).mean(), 1),
                         "pct_gt1": round(100 * (v > 1).mean(), 1), "mean_log1p": round(v.mean(), 3)})
    pd.DataFrame(rows).to_csv(OUT / "gene_by_cluster.csv", index=False)

    h = x[HORMONES]
    dom = h.idxmax(axis=1).where(h.max(axis=1) > 0, "none")
    hb = pd.crosstab(cl, dom).reindex(index=cats, columns=HORMONES + ["none"], fill_value=0)
    hb.insert(0, "cells", hb.sum(axis=1))
    hb.rename_axis("cluster").to_csv(OUT / "hormone_by_cluster.csv")

    pal = {c: matplotlib.colors.to_hex(plt.get_cmap("tab20")(i % 20)) for i, c in enumerate(cats)}
    fig, axs = plt.subplots(1, 3, figsize=(21, 7))
    for ax, (name, mask) in zip(axs, [("all cells", np.ones(len(cl), bool)), ("control", ~ko.to_numpy()),
                                       ("KO", ko.to_numpy())]):
        ax.scatter(u[~mask, 0], u[~mask, 1], c=GREY, s=2.5, linewidths=0, rasterized=True)
        idx = np.random.default_rng(SEED).permutation(np.where(mask)[0])
        ax.scatter(u[idx, 0], u[idx, 1], c=cl.map(pal).to_numpy()[idx], s=3.5, linewidths=0, rasterized=True)
        label_clusters(ax, u, cl)
        ax.set_title(f"{DAY} r{RES}: {name}"); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(OUT / "umap_clusters.png", dpi=100); plt.close(fig)

    ncol = 4
    nrow = int(np.ceil(len(GENES) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5.5 * ncol, 5.2 * nrow), squeeze=False)
    for ax, g in zip(axs.flat, GENES):
        v = x[g].to_numpy()
        vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1
        o = np.argsort(v, kind="stable")
        s = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=3, cmap=FEAT_CMAP, vmin=0, vmax=vmax, linewidths=0, rasterized=True)
        plt.colorbar(s, ax=ax, fraction=0.04, pad=0.01)
        label_clusters(ax, u, cl)
        ax.set_title(g); ax.set_xticks([]); ax.set_yticks([])
    for ax in list(axs.flat)[len(GENES):]:
        ax.axis("off")
    fig.suptitle(f"{DAY}: log1p, cluster numbers at r{RES}"); fig.tight_layout()
    fig.savefig(OUT / "feature_genes.png", dpi=90); plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
