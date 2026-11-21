"""Step 7.4i: D15 + D22 (alpha_beta_EC subset): unbiased classification of the KO vs control trajectory genes by their
expression along the beta, alpha, EC2 and EC1 lineages in control and KO (the 7.3h design on the 7.4 trajectory; after
7.4g and 7.4h).

Decisions (user, 2026-11-21): repeat the D15 analysis on D15 + D22; genes = only the KO vs control trajectory results:
7.4g stop calls (lost / gained), 7.4g "on after stop, missing in KO", 7.4h window calls (higher / lower in KO, any
lineage, window or stage). CellRank programs (pseudotime kernel, 7.4d) are labels only. KO-line-wide genes (called in
>= 2/3 of the 7.4h tests) stay in and are flagged.
  curves    per lineage (7.4a, scVI latent): cells with curve weight >= W_MIN and a label on the lineage's path; NBIN
            equal-width pseudotime bins; per genotype and bin the curve-weighted mean log1p expression (SoupX counts,
            normalised to 10,000, log1p), 3-bin moving average; bins with < MIN_CELLS cells of a genotype are "not
            reached" and left out (KO cells barely reach D22 SC_beta / SC_alpha). Pseudotime runs from D15 into D22, so
            the late part of each curve is D22 (a different sequencing batch)
  features  the kept bins of the 8 curves, z-scored per gene across all of them
  clusters  k-means (n_init 50, seed 1234), k = K_MIN..K_MAX by silhouette; ARI vs Ward and vs 5 other seeds; optional
            extra k values as a comma-separated argument (files with the suffix _k<k>)
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine_<subset>/trajectory/ (reruns overwrite):
  gene_pattern_curves.csv, gene_patterns[_k<k>].csv, gene_pattern_clusters[_k<k>].csv, gene_pattern_k.csv,
  gene_pattern_split_k<best>_to_k<k>.csv, gene_patterns_heatmap[_k<k>].png, gene_patterns_curves[_k<k>].png
Run: sbatch analysis/scripts/7.4i_gene_patterns.sb [k,k,...]   (CELLRANK_ENV python; INT_SUBSET, default alpha_beta_EC)
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
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_samples, silhouette_score

OUT_DIR = Path(os.environ["OUT_DIR"])
SUBSET = os.environ.get("INT_SUBSET", "alpha_beta_EC")
INT = OUT_DIR / "day_analysis" / f"integration_D15_D22_endocrine_{SUBSET}"
TR = INT / "trajectory"
LINEAGES = {"beta": "SC_beta", "alpha": "SC_alpha", "EC2": "SC_EC2", "EC1": "SC_EC1"}
PROGS = ["beta_program", "alpha_program", "EC2_program", "EC1_program", "EC1_KO_program", "EC2_KO_program"]
GENOS = ["ctr", "KO"]
W_MIN, NBIN, MIN_CELLS = 0.5, 10, 15
K_MIN, K_MAX, SEED = 4, 20, 1234
COL = {("ctr", "beta"): "#2166ac", ("KO", "beta"): "#92c5de", ("ctr", "alpha"): "#b2182b",
       ("KO", "alpha"): "#f4a582", ("ctr", "EC2"): "#1b7837", ("KO", "EC2"): "#a6dba0",
       ("ctr", "EC1"): "#762a83", ("KO", "EC1"): "#c2a5cf"}


def gene_sources():
    src = {}

    def add(genes, tag):
        for g in genes:
            src.setdefault(g, set()).add(tag)
    f = pd.read_csv(TR / "cellrank_lost_in_KO.csv", usecols=["gene", "lineage", "call"])
    f = f[f.call != "-"]
    for (lin, call), d in f.groupby(["lineage", "call"]):
        add(d.gene, f"7.4g {lin} stop: {call}")
    o = pd.read_csv(TR / "cellrank_on_after_KO_stop.csv", usecols=["gene", "lineage", "call"])
    for (lin, call), d in o[o.call == "on after stop, missing in KO"].groupby(["lineage", "call"]):
        add(d.gene, f"7.4g {lin}: {call}")
    w = pd.read_csv(TR / "ko_vs_ctr_windows.csv.gz", usecols=["gene", "lineage", "window", "stage", "call"] + PROGS)
    c = w[w.call != "-"]
    for (lin, call), d in c.groupby(["lineage", "call"]):
        add(d.gene.unique(), f"7.4h {lin}: {call}")
    prog = w.drop_duplicates("gene").set_index("gene")[PROGS]
    n_win = w.groupby(["lineage", "window", "stage"]).ngroups
    nc = c.groupby(["gene", "call"]).size().unstack(fill_value=0)
    line_wide = nc.index[(nc.max(axis=1) >= 2 / 3 * n_win)]
    print(f"KO-line-wide genes (called in >= 2/3 of {n_win} tests): {len(line_wide)}")
    return src, prog, set(line_wide)


def main():
    src, prog, line_wide = gene_sources()
    a = ad.read_h5ad(INT / "merged_counts.h5ad")
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sl = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0)
    a = a[sl.index].copy()
    genes = [g for g in a.var_names if g in src]
    print(f"genes in any category: {len(src)}; in the expression matrix: {len(genes)}")
    tags = pd.Series({g: "; ".join(sorted(src[g])) for g in genes})
    cat_counts = pd.Series([t for g in genes for t in src[g]]).value_counts()
    print(cat_counts.to_string())
    X = a[:, genes].X
    X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    geno = a.obs.genotype.astype(str).to_numpy()
    ln = pd.read_csv(TR / "slingshot_scvi_lineages.csv")
    paths = {r.lineage: r.clusters.split(" -> ") for r in ln.itertuples()}
    end = {k: v[-1] for k, v in paths.items()}

    # ---- curves
    rows, feats, fnames, keepcol = [], [], [], []
    for lin, fate in LINEAGES.items():
        L = [k for k, v in end.items() if v == fate][0]
        pt, w = sl[f"pt_{L}"].to_numpy(), sl[f"w_{L}"].fillna(0).to_numpy()
        on = (w >= W_MIN) & ~np.isnan(pt) & sl.label.isin(paths[L]).to_numpy()
        edges = np.linspace(np.nanmin(pt[on]), np.nanmax(pt[on]), NBIN + 1)
        b = np.clip(np.digitize(pt, edges) - 1, 0, NBIN - 1)
        for g in GENOS:
            M = np.full((len(genes), NBIN), np.nan); n = np.zeros(NBIN, int)
            for i in range(NBIN):
                m = on & (geno == g) & (b == i)
                n[i] = m.sum()
                if n[i] >= MIN_CELLS:
                    ww = w[m]
                    M[:, i] = np.asarray(X[m].T @ ww).ravel() / ww.sum()
            S = pd.DataFrame(M.T).rolling(3, center=True, min_periods=1).mean().to_numpy().T
            S[:, n < MIN_CELLS] = np.nan
            for i in range(NBIN):
                rows.append(pd.DataFrame({"gene": genes, "lineage": lin, "genotype": g, "bin": i + 1,
                                          "pseudotime": round((edges[i] + edges[i + 1]) / 2, 1), "cells": int(n[i]),
                                          "mean": np.round(M[:, i], 4), "smoothed": np.round(S[:, i], 4)}))
            feats.append(S); fnames += [f"{lin}|{g}|{i + 1}" for i in range(NBIN)]
            keepcol += list(n >= MIN_CELLS)
            print(f"{lin} {g}: cells per bin {n.tolist()}")
    pd.concat(rows).to_csv(TR / "gene_pattern_curves.csv", index=False)
    F = np.hstack(feats); keepcol = np.array(keepcol)
    print("bins not reached (left out):", [f for f, k in zip(fnames, keepcol) if not k])
    Fk = F[:, keepcol]
    level = np.nanmax(Fk, axis=1)
    Z = (Fk - Fk.mean(1, keepdims=True)) / (Fk.std(1, keepdims=True) + 1e-9)
    ok = (level > 0) & np.isfinite(Z).all(1)
    print("genes with a usable profile:", int(ok.sum()), "of", len(genes))
    Z, gk = Z[ok], np.array(genes)[ok]

    # ---- choose k
    krows, fits = [], {}
    Zw = linkage(Z, method="ward")
    for k in range(K_MIN, K_MAX + 1):
        km = KMeans(n_clusters=k, n_init=50, random_state=SEED).fit(Z)
        fits[k] = km.labels_
        s = silhouette_score(Z, km.labels_, random_state=SEED, sample_size=min(len(Z), 5000))
        ari_w = adjusted_rand_score(km.labels_, fcluster(Zw, k, "maxclust"))
        ari_s = np.mean([adjusted_rand_score(km.labels_, KMeans(n_clusters=k, n_init=50, random_state=SEED + j)
                                             .fit(Z).labels_) for j in range(1, 6)])
        krows.append({"k": k, "silhouette": round(s, 4), "ARI_vs_ward": round(ari_w, 3), "ARI_vs_seeds": round(ari_s, 3)})
    kt = pd.DataFrame(krows); kt.to_csv(TR / "gene_pattern_k.csv", index=False)
    print(kt.to_string(index=False))
    best = int(kt.sort_values(["silhouette", "k"], ascending=[False, True]).k.iloc[0])
    ks = sorted({best} | ({int(x) for x in sys.argv[1].split(",")} if len(sys.argv) > 1 else set()))   # best always
    print(f"silhouette-best k = {best}; writing k = {ks}")
    ref = None
    for k in ks:
        sfx = "" if k == best else f"_k{k}"   # the silhouette-best k keeps the plain file names
        lab_ = fits[k]
        # order clusters by where the control curves peak (beta, alpha, EC order of columns)
        order = pd.Series(np.argmax(Z, axis=1)).groupby(lab_).median().sort_values().index
        remap = {old: i + 1 for i, old in enumerate(order)}
        cl = np.array([remap[x] for x in lab_])
        sil = silhouette_samples(Z, cl)
        print(f"\n===== k = {k}")

        gp = pd.DataFrame({"gene": gk, "cluster": cl, "silhouette": sil.round(3), "level_max": level[ok].round(3),
                           "line_wide": np.isin(gk, list(line_wide)), "sources": tags.reindex(gk).to_numpy()})
        gp = gp.merge(prog, left_on="gene", right_index=True, how="left")
        gp = gp.sort_values(["cluster", "silhouette"], ascending=[True, False], kind="stable")
        gp.to_csv(TR / f"gene_patterns{sfx}.csv", index=False)
        crow = []
        for c, d in gp.groupby("cluster"):
            r = {"cluster": c, "genes": len(d), "mean_silhouette": round(d.silhouette.mean(), 3),
                 "line_wide": int(d.line_wide.sum())}
            for p in PROGS:
                r[p] = int((d[p] == "yes").sum())
            st = pd.Series([t for s in d.sources for t in s.split("; ")]).value_counts()
            r["top_sources"] = "; ".join(f"{t} {v}" for t, v in st.head(6).items())
            r["top_genes"] = ", ".join(d.sort_values("level_max", ascending=False).gene.head(20))
            crow.append(r)
        ct = pd.DataFrame(crow); ct.to_csv(TR / f"gene_pattern_clusters{sfx}.csv", index=False)
        print(ct.drop(columns=["top_sources"]).to_string(index=False))
        for r in crow:
            print(f"\ncluster {r['cluster']} ({r['genes']}): {r['top_sources']}\n  {r['top_genes']}")

        # ---- figures
        Zfull = np.full((len(gk), len(fnames)), np.nan)
        Zfull[:, keepcol] = Z
        o = np.lexsort((-sil, cl))
        fig, ax = plt.subplots(figsize=(14, 11))
        im = ax.imshow(np.ma.masked_invalid(Zfull[o]), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                       interpolation="nearest")
        im.cmap.set_bad("#d9d9d9")
        for j in range(1, 2 * len(LINEAGES)):
            ax.axvline(j * NBIN - 0.5, c="k", lw=1.2 if j % 2 == 0 else 0.5)
        ax.set_xticks([NBIN * j + NBIN / 2 - 0.5 for j in range(2 * len(LINEAGES))])
        ax.set_xticklabels([f"{lin} {g}" for lin in LINEAGES for g in GENOS])
        bounds = np.cumsum(np.bincount(cl[o])[1:])
        for y in bounds[:-1]:
            ax.axhline(y - 0.5, c="k", lw=0.6)
        ax.set_yticks(bounds - np.bincount(cl[o])[1:] / 2); ax.set_yticklabels([f"C{c}" for c in range(1, k + 1)])
        ax.set_title(f"{len(gk)} trajectory genes, k-means k = {k} (z-score per gene across the 8 curves; grey = KO not "
                     "reached; pseudotime left to right within each curve)", fontsize=10)
        plt.colorbar(im, ax=ax, fraction=0.02); fig.tight_layout()
        fig.savefig(TR / f"gene_patterns_heatmap{sfx}.png", dpi=90); plt.close(fig)
        nc = 4; nr = int(np.ceil(k / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(5 * nc, 3.6 * nr), squeeze=False)
        x = np.arange(NBIN)
        for c in range(1, k + 1):
            ax = axs.flat[c - 1]; m = cl == c
            for j, (lin, g) in enumerate([(lin, g) for lin in LINEAGES for g in GENOS]):
                y = np.nanmean(Zfull[m, j * NBIN:(j + 1) * NBIN], axis=0)
                ax.plot(x + list(LINEAGES).index(lin) * (NBIN + 1), y, c=COL[(g, lin)],
                        lw=2, ls="-" if g == "ctr" else "--", label=f"{lin} {g}")
            ax.axhline(0, c="#ccc", lw=0.6)
            ax.set_title(f"C{c}: {m.sum()} genes (line-wide {int(gp[gp.cluster == c].line_wide.sum())})", fontsize=10)
            ax.set_xticks([(i + 0.5) * NBIN + i for i in range(len(LINEAGES))]); ax.set_xticklabels(list(LINEAGES))
        axs.flat[0].legend(fontsize=7, ncol=2)
        for ax in list(axs.flat)[k:]:
            ax.axis("off")
        fig.suptitle("Mean z-scored curves per cluster (solid = control, dashed = KO)", fontsize=12)
        fig.tight_layout(); fig.savefig(TR / f"gene_patterns_curves{sfx}.png", dpi=85); plt.close(fig)
        if k == best:
            ref = pd.Series(cl, index=gk)
    # how each finer k splits the silhouette-best clusters (needs the best k among ks or computed here)
    lab_b = fits[best]
    order = pd.Series(np.argmax(Z, axis=1)).groupby(lab_b).median().sort_values().index
    rb = {old: i + 1 for i, old in enumerate(order)}
    ref = pd.Series([rb[x] for x in lab_b], index=gk)
    for k in ks:
        if k == best:
            continue
        g = pd.read_csv(TR / f"gene_patterns_k{k}.csv").set_index("gene").cluster.reindex(gk)
        x = pd.crosstab(ref.rename(f"k{best} cluster"), g.rename(f"k{k} cluster"))
        x.to_csv(TR / f"gene_pattern_split_k{best}_to_k{k}.csv")
        print(f"\nk = {best} clusters (rows) split into k = {k} clusters (columns):\n", x.to_string())
    print("\nwrote", TR, "gene_patterns*.csv, gene_pattern_clusters*.csv, heatmaps, curves")


if __name__ == "__main__":
    main()
