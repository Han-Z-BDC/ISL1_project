"""Step 7.3h: unbiased classification of the D15 trajectory genes by their expression along the beta, alpha and EC
lineages in control and KO (after 7.3f and 7.3g).

Decisions (user, 2026-11-21): classify every gene identified in any category of the trajectory analysis by the shape of
its expression along the three lineages x two genotypes, to find patterns such as "on toward beta only", "stopped in
KO", "stopped in KO on every lineage" or "switched to the EC lineage in KO". Unbiased clustering, no rule-based classes;
binned curves now (do not wait for tradeSeq 7.3e; its smoothers can replace the bins later).
  genes     (user, 2026-11-21: only the KO vs control trajectory results) union of: 7.3f pseudotime-matched stop
            calls (lost / lost (one kernel) / gained; label matches not used), 7.3f "on after stop, missing in KO"
            (genes that only rise in control after the stop, with no KO difference, are left out), and 7.3g window
            calls (higher / lower in KO, any lineage or window). CellRank programs are not a source; their membership
            is kept as labels. The source categories are kept per gene. KO-line-wide genes (called in >= LINE_WIDE of
            the tested 7.3g windows) stay in and are flagged.
  curves    per lineage (Slingshot 7.3a PCs 1-30; beta = lineage ending in pre_beta2, alpha = pre_alpha, EC = pre_EC):
            cells with curve weight >= W_MIN, NBIN equal-width pseudotime bins over the lineage; per genotype and bin
            the curve-weighted mean log1p expression (7.1d X), then a 3-bin moving average within the lineage. A bin
            with < MIN_CELLS cells of a genotype is "not reached" (empty for every gene: KO cells never reach the late
            beta and alpha bins); those columns are left out of the clustering, not filled in.
  features  the kept bins of the 6 curves (beta ctr | beta KO | alpha ctr | alpha KO | EC ctr | EC KO), z-scored per
            gene across all of them together, so differences between lineages and between genotypes count, not the
            gene's overall level
  clusters  k-means (n_init 50, seed 1234) for k = K_MIN..K_MAX; k with the highest mean silhouette; agreement with
            Ward hierarchical clustering at the same k (adjusted Rand index) and with k-means from 5 other seeds
            reported; per-gene silhouette kept
Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
  gene_pattern_curves.csv     gene, lineage, genotype, bin, pseudotime, cells, mean, smoothed
  gene_patterns.csv           gene, cluster, silhouette, source categories, program flags, line-wide flag, level
  gene_pattern_clusters.csv   cluster: genes, mean silhouette, source / program / line-wide composition, top genes
  gene_pattern_k.csv          k, silhouette, ARI vs Ward, ARI vs other seeds
  gene_patterns_heatmap.png, gene_patterns_curves.png
  extra k   (user, 2026-11-21: try k 8-10) optional argument, comma-separated k values: each writes its own
            files with the suffix _k<k> (the silhouette-best k keeps the plain names) and a table of how the
            silhouette-best clusters split into it (gene_pattern_split_k<best>_to_k<k>.csv)
  pseudobulk (user, 2026-11-21) with GENESET=pseudobulk: genes = the 7.3l pseudobulk calls with curves on the three
            lineages (curve test on any lineage, or past-the-stop missing / kept on in KO; KO-state-only genes are a
            separate category, not clustered); sources = those calls; line-wide = curve call on all three lineages
            in the same direction; same curves and clustering; writes to trajectory/pseudobulk/gene_patterns/
Run: sbatch analysis/scripts/7.3h_gene_patterns.sb [k,k,...]   (CELLRANK_ENV python)
     GENESET=pseudobulk: sbatch --export=ALL,GENESET=pseudobulk analysis/scripts/7.3h_gene_patterns.sb 10
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
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_samples, silhouette_score

OUT_DIR = Path(os.environ["OUT_DIR"])
DD = OUT_DIR / "day_analysis" / "D15"
SUB, TR = DD / "subcluster_endocrine", DD / "trajectory"
GENESET = os.environ.get("GENESET", "")
OUTD = TR / "pseudobulk" / "gene_patterns" if GENESET == "pseudobulk" else TR
LINEAGES = {"beta": "pre_beta2", "alpha": "pre_alpha", "EC": "pre_EC"}
GENOS = ["ctr", "KO"]
W_MIN, NBIN, MIN_CELLS, LINE_WIDE = 0.5, 10, 15, 15
K_MIN, K_MAX, SEED = 4, 20, 1234
COL = {("ctr", "beta"): "#2166ac", ("KO", "beta"): "#92c5de", ("ctr", "alpha"): "#b2182b",
       ("KO", "alpha"): "#f4a582", ("ctr", "EC"): "#1b7837", ("KO", "EC"): "#a6dba0"}


def gene_sources_pseudobulk():
    g = pd.read_csv(TR / "pseudobulk" / "pb_genes.csv", keep_default_na=False)
    g = g[g.in_curves.astype(str) == "True"]
    src = {}
    for r in g.itertuples():
        tags = [f"pb curve {t}" for t in r.curve.split("; ") if t] + [f"pb past stop {t}" for t in r.after_stop.split("; ") if t]
        tags += [f"pb KO state {t}" for t in r.ko_state.split("; ") if t]
        src[r.gene] = set(tags)
    w = pd.read_csv(TR / "ko_vs_ctr_windows.csv.gz", usecols=["gene", "beta_program", "alpha_program", "EC_program",
                                                              "EC_KO_program"])
    prog = w.drop_duplicates("gene").set_index("gene")
    c = pd.read_csv(TR / "pseudobulk" / "pb_curve.csv")
    c = c[c.call.astype(str) == "True"]
    d = c.groupby("gene").agg(n=("lineage", "nunique"), dirs=("direction", lambda x: set(x)))
    line_wide = set(d.index[(d.n == 3) & d.dirs.apply(lambda x: len(x) == 1 and list(x)[0] in ("higher in KO", "lower in KO"))])
    print(f"pseudobulk gene set: {len(src)} genes; line-wide (curve call on all three lineages, same direction): {len(line_wide)}")
    return src, prog, line_wide


def gene_sources():
    if GENESET == "pseudobulk":
        return gene_sources_pseudobulk()
    src = {}

    def add(genes, tag):
        for g in genes:
            src.setdefault(g, set()).add(tag)
    f = pd.read_csv(TR / "cellrank_lost_in_KO.csv", usecols=["gene", "lineage", "match", "step", "call"])
    f = f[(f.match == "pseudotime") & (f.call != "-")]
    for (lin, call), d in f.groupby(["lineage", "call"]):
        add(d.gene, f"7.3f {lin} stop: {call}")
    o = pd.read_csv(TR / "cellrank_on_after_KO_stop.csv", usecols=["gene", "lineage", "call"])
    for (lin, call), d in o[o.call == "on after stop, missing in KO"].groupby(["lineage", "call"]):
        add(d.gene, f"7.3f {lin}: {call}")
    w = pd.read_csv(TR / "ko_vs_ctr_windows.csv.gz",
                    usecols=["gene", "lineage", "window", "call", "beta_program", "alpha_program", "EC_program",
                             "EC_KO_program"])
    c = w[w.call != "-"]
    for (lin, call), d in c.groupby(["lineage", "call"]):
        add(d.gene.unique(), f"7.3g {lin}: {call}")
    prog = w.drop_duplicates("gene").set_index("gene")[["beta_program", "alpha_program", "EC_program", "EC_KO_program"]]
    n_win = w.groupby(["lineage", "window"]).ngroups
    nc = c.groupby(["gene", "call"]).size().unstack(fill_value=0)
    line_wide = nc.index[(nc.max(axis=1) >= LINE_WIDE)]
    print(f"KO-line-wide genes (called in >= {LINE_WIDE} of {n_win} windows): {len(line_wide)}")
    return src, prog, set(line_wide)


def main():
    OUTD.mkdir(parents=True, exist_ok=True)
    src, prog, line_wide = gene_sources()
    a = ad.read_h5ad(SUB / "D15_endocrine.h5ad")
    lab = pd.read_csv(TR / "cellrank_velocity_fates.csv", index_col=0, keep_default_na=False)
    a = a[lab.index].copy()
    genes = [g for g in a.var_names if g in src]
    print(f"genes in any category: {len(src)}; in the expression matrix: {len(genes)}")
    tags = pd.Series({g: "; ".join(sorted(src[g])) for g in genes})
    cat_counts = pd.Series([t for g in genes for t in src[g]]).value_counts()
    print(cat_counts.to_string())
    X = a[:, genes].X
    X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    geno = a.obs.genotype.astype(str).to_numpy()
    sl = pd.read_csv(TR / "slingshot_pc30_cells.csv", index_col=0).reindex(a.obs_names)
    ln = pd.read_csv(TR / "slingshot_pc30_lineages.csv")
    end = {r.lineage: r.clusters.split(" -> ")[-1] for r in ln.itertuples()}

    # ---- curves
    rows, feats, fnames, keepcol = [], [], [], []
    for lin, fate in LINEAGES.items():
        L = [k for k, v in end.items() if v == fate][0]
        pt, w = sl[f"pt_{L}"].to_numpy(), sl[f"w_{L}"].fillna(0).to_numpy()
        on = (w >= W_MIN) & ~np.isnan(pt)
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
    pd.concat(rows).to_csv(OUTD / "gene_pattern_curves.csv", index=False)
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
    kt = pd.DataFrame(krows); kt.to_csv(OUTD / "gene_pattern_k.csv", index=False)
    print(kt.to_string(index=False))
    best = int(kt.sort_values(["silhouette", "k"], ascending=[False, True]).k.iloc[0])
    ks = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [best]
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
        gp.to_csv(OUTD / f"gene_patterns{sfx}.csv", index=False)
        crow = []
        for c, d in gp.groupby("cluster"):
            r = {"cluster": c, "genes": len(d), "mean_silhouette": round(d.silhouette.mean(), 3),
                 "line_wide": int(d.line_wide.sum())}
            for p in ["beta_program", "alpha_program", "EC_program", "EC_KO_program"]:
                r[p] = int((d[p] == "both").sum())
            st = pd.Series([t for s in d.sources for t in s.split("; ")]).value_counts()
            r["top_sources"] = "; ".join(f"{t} {v}" for t, v in st.head(6).items())
            r["top_genes"] = ", ".join(d.sort_values("level_max", ascending=False).gene.head(20))
            crow.append(r)
        ct = pd.DataFrame(crow); ct.to_csv(OUTD / f"gene_pattern_clusters{sfx}.csv", index=False)
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
        for j in range(1, 6):
            ax.axvline(j * NBIN - 0.5, c="k", lw=1.2 if j % 2 == 0 else 0.5)
        ax.set_xticks([NBIN * j + NBIN / 2 - 0.5 for j in range(6)])
        ax.set_xticklabels([f"{lin} {g}" for lin in LINEAGES for g in GENOS])
        bounds = np.cumsum(np.bincount(cl[o])[1:])
        for y in bounds[:-1]:
            ax.axhline(y - 0.5, c="k", lw=0.6)
        ax.set_yticks(bounds - np.bincount(cl[o])[1:] / 2); ax.set_yticklabels([f"C{c}" for c in range(1, k + 1)])
        ax.set_title(f"{len(gk)} trajectory genes, k-means k = {k} (z-score per gene across the 6 curves; grey = KO not "
                     "reached; pseudotime left to right within each curve)", fontsize=10)
        plt.colorbar(im, ax=ax, fraction=0.02); fig.tight_layout()
        fig.savefig(OUTD / f"gene_patterns_heatmap{sfx}.png", dpi=90); plt.close(fig)
        nc = 4; nr = int(np.ceil(k / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(5 * nc, 3.6 * nr), squeeze=False)
        x = np.arange(NBIN)
        for c in range(1, k + 1):
            ax = axs.flat[c - 1]; m = cl == c
            for j, (lin, g) in enumerate([(lin, g) for lin in LINEAGES for g in GENOS]):
                y = np.nanmean(Zfull[m, j * NBIN:(j + 1) * NBIN], axis=0)
                ax.plot(x + (0 if lin == "beta" else NBIN + 1 if lin == "alpha" else 2 * NBIN + 2), y, c=COL[(g, lin)],
                        lw=2, ls="-" if g == "ctr" else "--", label=f"{lin} {g}")
            ax.axhline(0, c="#ccc", lw=0.6)
            ax.set_title(f"C{c}: {m.sum()} genes (line-wide {int(gp[gp.cluster == c].line_wide.sum())})", fontsize=10)
            ax.set_xticks([NBIN / 2, 1.5 * NBIN + 1, 2.5 * NBIN + 2]); ax.set_xticklabels(["beta", "alpha", "EC"])
        axs.flat[0].legend(fontsize=7, ncol=2)
        for ax in list(axs.flat)[k:]:
            ax.axis("off")
        fig.suptitle("Mean z-scored curves per cluster (solid = control, dashed = KO)", fontsize=12)
        fig.tight_layout(); fig.savefig(OUTD / f"gene_patterns_curves{sfx}.png", dpi=85); plt.close(fig)
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
        g = pd.read_csv(OUTD / f"gene_patterns_k{k}.csv").set_index("gene").cluster.reindex(gk)
        x = pd.crosstab(ref.rename(f"k{best} cluster"), g.rename(f"k{k} cluster"))
        x.to_csv(OUTD / f"gene_pattern_split_k{best}_to_k{k}.csv")
        print(f"\nk = {best} clusters (rows) split into k = {k} clusters (columns):\n", x.to_string())
    print("\nwrote", OUTD, "gene_patterns*.csv, gene_pattern_clusters*.csv, heatmaps, curves")


if __name__ == "__main__":
    main()
