"""Helper: genes' expression along every lineage of a D15 + D22 Slingshot run (7.4a), control vs KO. Not a pipeline step.
Works for the full integration (8 lineages) and the alpha_beta_EC subset (4 lineages); any gene (user, 2026-11-21:
plot genes in the 8 D15 + D22 trajectories). Cells on a lineage = curve weight >= 0.5 and a label on the lineage's
path (as 7.4g-i); expression = the integration's counts (merged_counts.h5ad) normalised to 10,000 per cell (totals over
all genes), log1p; 10 equal-width pseudotime bins per lineage, curve-weighted mean, 3-bin smoothing, < 15 cells of a
genotype = no line; background shading = share of D22 cells in the bin. All lineages of a gene share one y-axis.
Usage (CELLRANK_ENV or TOOLS python):
  python plot_gene_lineages_D15_D22.py GENE1,GENE2 <integration folder> <output folder>
  integration folder = ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine[_alpha_beta_EC] (merged_counts.h5ad and
  trajectory/slingshot_scvi_*.csv)
Writes <output folder>/<gene>_D15_D22_lineages.png
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

genes, INT, DEST = sys.argv[1].split(","), Path(sys.argv[2]), Path(sys.argv[3])
TR = INT / "trajectory"
W_MIN, NBIN, MIN_CELLS = 0.5, 10, 15
COL = {"SC_beta": "#2166ac", "SC_alpha": "#b2182b", "SC_EC2": "#1b7837", "SC_EC1": "#5aae61", "SC_delta": "#8c510a",
       "D15 KO_only": "#762a83", "D15 ENP_PHOX2A+": "#e08214", "GHRL+": "#c51b7d"}

sl = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0)
a = ad.read_h5ad(INT / "merged_counts.h5ad")
a = a[sl.index]
tot = np.asarray(a.X.sum(1)).ravel()
X = a[:, genes].X
X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
X = pd.DataFrame(np.log1p(X / tot[:, None] * 1e4), index=sl.index, columns=genes)
geno = sl.genotype.astype(str).to_numpy()
d22 = (sl.stage.astype(str) == "D22").to_numpy()
ln = pd.read_csv(TR / "slingshot_scvi_lineages.csv")
paths = {r.lineage: r.clusters.split(" -> ") for r in ln.itertuples()}
nc = 4
nr = int(np.ceil(len(paths) / nc))
DEST.mkdir(parents=True, exist_ok=True)
for G in genes:
    x = X[G].to_numpy()
    fig, axs = plt.subplots(nr, nc, figsize=(22, 5.2 * nr), sharey=True, squeeze=False)
    for ax in axs.ravel()[len(paths):]:
        ax.axis("off")
    for ax, (L, path) in zip(axs.ravel(), paths.items()):
        pt, w = sl[f"pt_{L}"].to_numpy(), sl[f"w_{L}"].fillna(0).to_numpy()
        on = (w >= W_MIN) & sl.label.isin(path).to_numpy() & ~np.isnan(pt)
        edges = np.linspace(np.nanmin(pt[on]), np.nanmax(pt[on]), NBIN + 1)
        b = np.clip(np.digitize(pt, edges) - 1, 0, NBIN - 1)
        mid = (edges[:-1] + edges[1:]) / 2
        for i in range(NBIN):
            m = on & (b == i)
            if m.any():
                ax.axvspan(edges[i], edges[i + 1], color="#f2e6c9", alpha=float(d22[m].mean()), lw=0)
        c = COL.get(path[-1], "#444")
        for g, ls in [("ctr", "-"), ("KO", "--")]:
            M, n = np.full(NBIN, np.nan), np.zeros(NBIN, int)
            for i in range(NBIN):
                m = on & (geno == g) & (b == i); n[i] = m.sum()
                if n[i]:
                    M[i] = np.average(x[m], weights=w[m])
            S = pd.Series(np.where(n >= MIN_CELLS, M, np.nan)).rolling(3, center=True, min_periods=1).mean().to_numpy()
            S[n < MIN_CELLS] = np.nan
            ax.plot(mid, S, ls=ls, c=c, lw=2.4, label=f"{'control' if g == 'ctr' else 'KO'} ({int((on & (geno == g)).sum())})")
            ok = n > 0
            ax.scatter(mid[ok], M[ok], s=np.clip(n[ok], 1, 3000) / 12 + 8, c=c, alpha=0.9 if g == "ctr" else 0.35,
                       edgecolors="k" if g == "ctr" else "none", linewidths=0.4, zorder=3)
        short = [p.replace("D15 ", "") for p in path]
        ax.set_title(f"{L} -> {path[-1]}\n{' > '.join(short)}", fontsize=8.5)
        ax.set_xlabel("Slingshot pseudotime (shaded = share of D22 cells)"); ax.legend(fontsize=8, frameon=False)
    for ax in axs[:, 0]:
        ax.set_ylabel(f"{G} mean log1p (counts / 10,000)")
    fig.suptitle(f"{G} along the {len(paths)} D15 + D22 lineages ({INT.name}), control (solid) vs KO (dashed); line = "
                 f"3-bin smoothed, dot = bin mean (size = cells); bins with < {MIN_CELLS} cells of a genotype have no line",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(DEST / f"{G}_D15_D22_lineages.png", dpi=90); plt.close(fig)
    print("wrote", DEST / f"{G}_D15_D22_lineages.png")
