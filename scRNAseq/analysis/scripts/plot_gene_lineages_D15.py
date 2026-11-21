"""Helper: genes' expression along the D15 Slingshot lineages (beta, alpha, EC), control vs KO, binned as 7.3h
(curve weight >= 0.5, 10 equal-width bins, curve-weighted mean log1p of the 7.1d X, 3-bin smoothing, < 15 cells = no
line). Computed from the expression, so any gene works. Not a pipeline step.
Usage (TOOLS python, LD_LIBRARY_PATH=$(dirname $TOOLS)/lib): python plot_gene_lineages_D15.py GENE1,GENE2
  ${OUT_DIR}/day_analysis/D15 <output folder>
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
LIN = {"beta": "pre_beta2", "alpha": "pre_alpha", "EC": "pre_EC"}
COL = {"beta": "#2166ac", "alpha": "#b2182b", "EC": "#1b7837"}
a = ad.read_h5ad(SUB / "D15_endocrine.h5ad", backed="r")
lab = pd.read_csv(TR / "cellrank_velocity_fates.csv", index_col=0, keep_default_na=False)
idx = a.obs_names.get_indexer(lab.index)
X = a[:, genes].X
X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
X = pd.DataFrame(X, index=a.obs_names, columns=genes).iloc[idx]
geno = a.obs.genotype.astype(str).to_numpy()[idx]
sl = pd.read_csv(TR / "slingshot_pc30_cells.csv", index_col=0).reindex(lab.index)
ln = pd.read_csv(TR / "slingshot_pc30_lineages.csv")
paths = {r.lineage: r.clusters.split(" -> ") for r in ln.itertuples()}
for G in genes:
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.2), sharey=True)
    for ax, (lin, fate) in zip(axs, LIN.items()):
        L = [k for k, v in paths.items() if v[-1] == fate][0]
        pt, w = sl[f"pt_{L}"].to_numpy(), sl[f"w_{L}"].fillna(0).to_numpy()
        on = (w >= W_MIN) & ~np.isnan(pt)
        edges = np.linspace(np.nanmin(pt[on]), np.nanmax(pt[on]), NBIN + 1)
        b = np.clip(np.digitize(pt, edges) - 1, 0, NBIN - 1); mid = (edges[:-1] + edges[1:]) / 2
        for g, ls in [("ctr", "-"), ("KO", "--")]:
            M, n = np.full(NBIN, np.nan), np.zeros(NBIN, int)
            for i in range(NBIN):
                m = on & (geno == g) & (b == i); n[i] = m.sum()
                if n[i]:
                    M[i] = np.average(X[G].to_numpy()[m], weights=w[m])
            S = pd.Series(np.where(n >= MIN_CELLS, M, np.nan)).rolling(3, center=True, min_periods=1).mean().to_numpy()
            S[n < MIN_CELLS] = np.nan
            ax.plot(mid, S, ls=ls, c=COL[lin], lw=2.4, label="control" if g == "ctr" else "KO")
            ok = n > 0
            ax.scatter(mid[ok], M[ok], s=np.clip(n[ok], 1, 2000) / 6 + 8, c=COL[lin],
                       alpha=0.9 if g == "ctr" else 0.35, edgecolors="k" if g == "ctr" else "none", linewidths=0.4, zorder=3)
        ax.set_title(f"{lin}: {' > '.join(paths[L])}", fontsize=9)
        ax.set_xlabel("Slingshot pseudotime (D15, PCs 1-30)"); ax.legend(fontsize=9, frameon=False)
    axs[0].set_ylabel(f"{G} mean log1p (line = 3-bin smoothed; dot = bin mean, size = cells)")
    fig.suptitle(f"{G} along the D15 lineages, control (solid) vs KO (dashed); bins with < {MIN_CELLS} cells of a "
                 "genotype have no line", fontsize=11)
    fig.tight_layout(); DEST.mkdir(parents=True, exist_ok=True)
    fig.savefig(DEST / f"{G}_D15_lineages.png", dpi=95); plt.close(fig)
    print("wrote", DEST / f"{G}_D15_lineages.png")
