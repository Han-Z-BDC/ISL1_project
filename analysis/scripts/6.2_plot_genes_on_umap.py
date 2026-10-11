"""Step 6.2: plot genes on the clustering UMAP, next to a chosen Leiden resolution.

Reads 6.1's clustered_log_unintegrated.h5ad (backed; only the needed genes are loaded) and
writes ${OUT_DIR}/clustering/plots/genes_<tag>.png: the clusters at --resolution, then one
panel per gene (log1p_norm, high-expressing cells drawn on top), and prints, for every gene,
how its expressing cells are distributed across the clusters.

Light (no fitting): fine on the login node.
  source study_config/config.sh
  $TOOLS/python analysis/scripts/6.2_plot_genes_on_umap.py --h5ad $OUT_DIR/clustering/h5ad/clustered_log_unintegrated.h5ad \
      --resolution 0.2 --genes SST HHEX GHRL PPY --tag delta --outdir $OUT_DIR/clustering
"""
import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--resolution", required=True)
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--highlight", nargs="*", default=[],
                    help="genes to also show as a highlight panel: cells with expression > 0 in colour, "
                         "drawn larger on a grey background (for rare, scattered populations)")
    args = ap.parse_args()
    key = f"leiden_r{args.resolution.replace('.', 'p')}"

    a = ad.read_h5ad(args.h5ad, backed="r")
    u = np.asarray(a.obsm["X_umap"])
    cl = pd.Categorical(a.obs[key].astype(str), categories=sorted(a.obs[key].astype(str).unique(), key=int))
    geno = a.obs["genotype"].astype(str).to_numpy()
    genes = [g for g in args.genes if g in a.var_names]
    hl = [g for g in args.highlight if g in genes]
    X = a[:, genes].to_memory().layers["log1p_norm"]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)

    n = len(genes) + len(hl) + 1
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 4.8), sharex=True, sharey=True)
    cmap = plt.get_cmap("tab20")
    ax = axes[0]
    for i, c in enumerate(cl.categories):
        m = cl == c
        ax.scatter(u[m, 0], u[m, 1], s=1.2, color=cmap(i % 20), linewidths=0)
        ax.text(np.median(u[m, 0]), np.median(u[m, 1]), c, fontsize=10, weight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    ax.set_title(f"Leiden resolution {args.resolution}")
    for ax, (j, g) in zip(axes[1:], enumerate(genes)):
        x = X[:, j]
        o = np.argsort(x)
        im = ax.scatter(u[o, 0], u[o, 1], c=x[o], s=1.5, cmap="viridis", linewidths=0)
        ax.set_title(g)
        plt.colorbar(im, ax=ax, fraction=0.04)
    for ax, g in zip(axes[1 + len(genes):], hl):
        x = X[:, genes.index(g)]
        m = x > 0
        ax.scatter(u[~m, 0], u[~m, 1], s=1.2, c="#dddddd", linewidths=0)
        ax.scatter(u[m, 0], u[m, 1], s=7, c="#c0392b", linewidths=0)
        ax.set_title(f"{g} > 0 highlighted (n={m.sum():,})")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    out = Path(args.outdir) / "plots" / f"genes_{args.tag}.png"
    fig.savefig(out, dpi=90)
    print("wrote", out)

    # where do the expressing cells sit?
    for j, g in enumerate(genes):
        pos = X[:, j] > 0
        hi = X[:, j] > np.quantile(X[pos, j], 0.5) if pos.any() else pos
        t = pd.DataFrame({"cluster": cl, "pos": pos, "hi": hi, "geno": geno})
        s = t.groupby("cluster", observed=False).agg(n=("pos", "size"), n_pos=("pos", "sum"), n_hi=("hi", "sum"))
        s["pct_cluster_pos"] = (100 * s["n_pos"] / s["n"]).round(1)
        s["share_of_hi_cells"] = (100 * s["n_hi"] / max(hi.sum(), 1)).round(1)
        print(f"\n{g}: {pos.sum()} cells >0, {hi.sum()} above median of expressing cells; "
              f"hi cells ctr={int((hi & (geno == 'ctr')).sum())} KO={int((hi & (geno == 'KO')).sum())}")
        print(s.to_string())


if __name__ == "__main__":
    main()
