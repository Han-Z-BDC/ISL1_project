"""Helper: one gene's binned expression along the D15 + D22 lineages (7.4a, alpha_beta_EC subset), control vs KO,
from the 7.4i curves (gene_pattern_curves.csv, so only genes classified in 7.4i). Not a pipeline step.
Usage (TOOLS or CELLRANK_ENV python): python plot_gene_lineages_D15_D22.py <gene> <7.4 trajectory folder> <output folder>
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

G, TR, DEST = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
W_MIN, NBIN = 0.5, 10
LIN = {"beta": "SC_beta", "alpha": "SC_alpha", "EC2": "SC_EC2", "EC1": "SC_EC1"}
COL = {"beta": "#2166ac", "alpha": "#b2182b", "EC2": "#1b7837", "EC1": "#762a83"}
c = pd.read_csv(TR / "gene_pattern_curves.csv"); c = c[c.gene == G]
assert len(c), f"{G} not in gene_pattern_curves.csv"
sl = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0)
ln = pd.read_csv(TR / "slingshot_scvi_lineages.csv")
paths = {r.lineage: r.clusters.split(" -> ") for r in ln.itertuples()}
fig, axs = plt.subplots(1, 4, figsize=(22, 5.2), sharey=True)
for ax, (lin, fate) in zip(axs, LIN.items()):
    L = [k for k, v in paths.items() if v[-1] == fate][0]
    on = (sl[f"w_{L}"].fillna(0) >= W_MIN) & sl.label.isin(paths[L]) & sl[f"pt_{L}"].notna()
    pt = sl.loc[on, f"pt_{L}"]
    edges = np.linspace(pt.min(), pt.max(), NBIN + 1)
    b = np.clip(np.digitize(pt, edges) - 1, 0, NBIN - 1)
    d22 = pd.Series((sl.loc[on, "stage"] == "D22").to_numpy()).groupby(b).mean()
    for i, f in d22.items():
        ax.axvspan(edges[i], edges[i + 1], color="#f2e6c9", alpha=float(f), lw=0)
    for g, ls in [("ctr", "-"), ("KO", "--")]:
        d = c[(c.lineage == lin) & (c.genotype == g)].sort_values("bin")
        ax.plot(d.pseudotime, d.smoothed, ls=ls, c=COL[lin], lw=2.4, label=f"{'control' if g == 'ctr' else 'KO'}")
        ax.scatter(d.pseudotime, d["mean"], s=np.clip(d.cells, 1, 2000) / 12 + 8, c=COL[lin],
                   alpha=0.9 if g == "ctr" else 0.35, edgecolors="k" if g == "ctr" else "none", linewidths=0.4, zorder=3)
        for _, r in d[d.cells < 15].iterrows():
            ax.annotate(f"{r.cells}", (r.pseudotime, 0), fontsize=7, color="#888", ha="center")
    ax.set_title(f"{lin}: {' > '.join(p.replace('D15 ', '') for p in paths[L])}", fontsize=9)
    ax.set_xlabel("Slingshot pseudotime (shaded = share of D22 cells in the bin)")
    ax.legend(fontsize=9, frameon=False)
axs[0].set_ylabel(f"{G} mean log1p (line = 3-bin smoothed; dot = bin mean, size = cells)")
fig.suptitle(f"{G} along the D15 + D22 lineages (alpha_beta_EC subset), control (solid) vs KO (dashed); bins with "
             "< 15 cells of a genotype have no line", fontsize=11)
fig.tight_layout(); DEST.mkdir(parents=True, exist_ok=True)
fig.savefig(DEST / f"{G}_D15_D22_lineages.png", dpi=95)
print("wrote", DEST / f"{G}_D15_D22_lineages.png")
