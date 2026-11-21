"""Helper: genes' expression along every D15 Slingshot lineage (7.3a, PCs 1-30; 6 lineages), control vs KO, binned
as 7.3h (curve weight >= 0.5, 10 equal-width bins per lineage, curve-weighted mean log1p of the 7.1d X, 3-bin
smoothing, < 15 cells of a genotype = no line). Computed from the expression, so any gene works. Not a pipeline step.
All lineages of a gene share one y-axis (user, 2026-11-21: plot genes in all 6 D15 trajectories).
Usage (CELLRANK_ENV or TOOLS python): python plot_gene_lineages_D15.py GENE1,GENE2 ${OUT_DIR}/day_analysis/D15 <output folder>
Writes <output folder>/<gene>_D15_lineages.png
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

genes, DD, DEST = sys.argv[1].split(","), Path(sys.argv[2]), Path(sys.argv[3])
TR, SUB = DD / "trajectory", DD / "subcluster_endocrine"
W_MIN, NBIN, MIN_CELLS = 0.5, 10, 15
COL = {"pre_beta2": "#2166ac", "pre_alpha": "#b2182b", "pre_EC": "#1b7837", "ENP_PHOX2A+": "#e08214",
       "KO_only": "#762a83", "pre_delta": "#8c510a"}


def binned(x, pt, w, on, geno, g):
    edges = np.linspace(np.nanmin(pt[on]), np.nanmax(pt[on]), NBIN + 1)
    b = np.clip(np.digitize(pt, edges) - 1, 0, NBIN - 1)
    M, n = np.full(NBIN, np.nan), np.zeros(NBIN, int)
    for i in range(NBIN):
        m = on & (geno == g) & (b == i); n[i] = m.sum()
        if n[i]:
            M[i] = np.average(x[m], weights=w[m])
    S = pd.Series(np.where(n >= MIN_CELLS, M, np.nan)).rolling(3, center=True, min_periods=1).mean().to_numpy()
    S[n < MIN_CELLS] = np.nan
    return (edges[:-1] + edges[1:]) / 2, M, S, n


a = ad.read_h5ad(SUB / "D15_endocrine.h5ad", backed="r")
sl = pd.read_csv(TR / "slingshot_pc30_cells.csv", index_col=0)
idx = a.obs_names.get_indexer(sl.index)
assert (idx >= 0).all()
X = a[:, genes].X
X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
X = pd.DataFrame(X[idx], index=sl.index, columns=genes)
geno = sl.genotype.astype(str).to_numpy()
ln = pd.read_csv(TR / "slingshot_pc30_lineages.csv")
paths = {r.lineage: r.clusters.split(" -> ") for r in ln.itertuples()}
DEST.mkdir(parents=True, exist_ok=True)
for G in genes:
    x = X[G].to_numpy()
    fig, axs = plt.subplots(2, 3, figsize=(18, 10), sharey=True)
    for ax, (L, path) in zip(axs.ravel(), paths.items()):
        pt, w = sl[f"pt_{L}"].to_numpy(), sl[f"w_{L}"].fillna(0).to_numpy()
        on = (w >= W_MIN) & ~np.isnan(pt)
        c = COL.get(path[-1], "#444")
        for g, ls in [("ctr", "-"), ("KO", "--")]:
            mid, M, S, n = binned(x, pt, w, on, geno, g)
            ax.plot(mid, S, ls=ls, c=c, lw=2.4, label=f"{'control' if g == 'ctr' else 'KO'} ({int((on & (geno == g)).sum())})")
            ok = n > 0
            ax.scatter(mid[ok], M[ok], s=np.clip(n[ok], 1, 2000) / 6 + 8, c=c, alpha=0.9 if g == "ctr" else 0.35,
                       edgecolors="k" if g == "ctr" else "none", linewidths=0.4, zorder=3)
        ax.set_title(f"{L} -> {path[-1]}\n{' > '.join(path)}", fontsize=9)
        ax.set_xlabel("Slingshot pseudotime (D15, PCs 1-30)"); ax.legend(fontsize=8, frameon=False)
    for ax in axs[:, 0]:
        ax.set_ylabel(f"{G} mean log1p")
    fig.suptitle(f"{G} along the 6 D15 lineages, control (solid) vs KO (dashed); line = 3-bin smoothed, dot = bin mean "
                 f"(size = cells); bins with < {MIN_CELLS} cells of a genotype have no line", fontsize=11)
    fig.tight_layout(); fig.savefig(DEST / f"{G}_D15_lineages.png", dpi=90); plt.close(fig)
    print("wrote", DEST / f"{G}_D15_lineages.png")
