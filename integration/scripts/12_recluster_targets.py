"""Step 12 (re-clustering): k-means of the D15 ISL1 direct targets among the C1-C10 trajectory genes (user, 2026-11-24:
"only take all the D15 direct targets among C1-C10 and redo the k-means clustering").

Genes: the C1-C10 genes (scRNAseq 7.3h, k = 10) with an ISL1 site called at D15 (TSS +/- 1 kb or ABC enhancer;
  genes.tsv of step 12, bound_D15).
Features: exactly the 7.3h input -- the smoothed binned curves of gene_pattern_curves.csv (10 equal-width pseudotime
  bins on beta / alpha / EC x control / KO; bins with < 15 cells of a genotype left out), z-scored per gene across the
  six curves (population sd).
k-means as 7.3h: k 4-12, n_init 50, seed 1234; silhouette, agreement across 5 other seeds (ARI); new clusters numbered by
  the median column of the genes' maximum z (as 7.3h).
Outputs in <out dir>/recluster_D15_targets/: k.csv, clusters_k<k>.csv (gene, new cluster, old C1-C10, silhouette),
  crosstab_k<k>.csv (new x old), curves_k<k>.png (mean z per new cluster).
usage: python 12_recluster_targets.py <gene_pattern_curves.csv> <genes.tsv> <out dir>
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples, adjusted_rand_score

K_MIN, K_MAX, SEED, NBIN = 4, 12, 1234, 10
LIN, GENOS = ["beta", "alpha", "EC"], ["ctr", "KO"]
COL = {("ctr", "beta"): "#965a00", ("KO", "beta"): "#965a00", ("ctr", "alpha"): "#00b3c4",
       ("KO", "alpha"): "#00b3c4", ("ctr", "EC"): "#3c3cff", ("KO", "EC"): "#3c3cff"}


def main():
    cv_f, genes_f, od = sys.argv[1], sys.argv[2], Path(sys.argv[3]) / "recluster_D15_targets"
    od.mkdir(parents=True, exist_ok=True)
    g = pd.read_csv(genes_f, sep="\t", keep_default_na=False)
    g = g[(g.group != "background")]
    tgt = g[g.bound_D15.astype(str) == "True"]
    old = pd.Series(tgt.group.values, index=tgt.gene.values)
    print(f"C1-C10 genes: {len(g)}; D15 targets: {len(tgt)}")
    print(old.value_counts().reindex([f"C{i}" for i in range(1, 11)]).to_string())

    cv = pd.read_csv(cv_f, usecols=["gene", "lineage", "genotype", "bin", "smoothed"])
    cv = cv[cv.gene.isin(old.index)]
    cols = [(l, gt, b) for l in LIN for gt in GENOS for b in range(1, NBIN + 1)]
    F = cv.pivot_table(index="gene", columns=["lineage", "genotype", "bin"], values="smoothed", dropna=False)
    F = F.reindex(columns=pd.MultiIndex.from_tuples(cols))
    keep = F.notna().any(axis=0).to_numpy()            # bins with >= 15 cells (empty for every gene otherwise)
    Fk = F.loc[:, keep].to_numpy()
    Z = (Fk - Fk.mean(1, keepdims=True)) / (Fk.std(1, keepdims=True) + 1e-9)
    ok = np.isfinite(Z).all(1) & (np.nanmax(Fk, 1) > 0)
    Z, gk = Z[ok], F.index.to_numpy()[ok]
    print(f"genes with a usable profile: {len(gk)}; feature columns: {keep.sum()}")

    rows, fits = [], {}
    for k in range(K_MIN, K_MAX + 1):
        lab = KMeans(n_clusters=k, n_init=50, random_state=SEED).fit(Z).labels_
        fits[k] = lab
        s = silhouette_score(Z, lab, random_state=SEED)
        ari = np.mean([adjusted_rand_score(lab, KMeans(n_clusters=k, n_init=50, random_state=SEED + j).fit(Z).labels_)
                       for j in range(1, 6)])
        ari_old = adjusted_rand_score(old.reindex(gk).values, lab)
        rows.append({"k": k, "silhouette": round(s, 4), "ARI_vs_seeds": round(ari, 3), "ARI_vs_C1_C10": round(ari_old, 3)})
    kt = pd.DataFrame(rows); kt.to_csv(od / "k.csv", index=False); print(kt.to_string(index=False))
    best = int(kt.sort_values(["silhouette", "k"], ascending=[False, True]).k.iloc[0])
    print(f"silhouette-best k = {best}")

    Zfull = np.full((len(gk), len(cols)), np.nan); Zfull[:, keep] = Z
    for k in range(K_MIN, K_MAX + 1):
        lab = fits[k]
        order = pd.Series(np.argmax(Z, axis=1)).groupby(lab).median().sort_values().index
        remap = {o: i + 1 for i, o in enumerate(order)}
        cl = np.array([remap[x] for x in lab])
        sil = silhouette_samples(Z, cl)
        out = pd.DataFrame({"gene": gk, "cluster": cl, "old_cluster": old.reindex(gk).values, "silhouette": sil.round(3)})
        out = out.sort_values(["cluster", "silhouette"], ascending=[True, False], kind="stable")
        out.to_csv(od / f"clusters_k{k}.csv", index=False)
        x = pd.crosstab(out.cluster.rename("new"), out.old_cluster.rename("old"))
        x = x.reindex(columns=[f"C{i}" for i in range(1, 11)], fill_value=0)
        x.to_csv(od / f"crosstab_k{k}.csv")
        if k in (best, 10) or k == K_MIN:
            print(f"\n===== k = {k}: new clusters (rows) x old C1-C10 (columns)\n{x.to_string()}")
        nc = 4; nr = int(np.ceil(k / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 2.8 * nr), squeeze=False)
        for c in range(1, k + 1):
            ax = axs.flat[c - 1]; m = cl == c
            for j, (l, gt) in enumerate([(l, gt) for l in LIN for gt in GENOS]):
                y = np.nanmean(Zfull[m, j * NBIN:(j + 1) * NBIN], axis=0)
                off = LIN.index(l) * (NBIN + 1)
                ax.plot(np.arange(NBIN) + off, y, c=COL[(gt, l)], lw=1.8, ls="-" if gt == "ctr" else "--")
            ax.axhline(0, c="#ccc", lw=0.6)
            top = x.loc[c].sort_values(ascending=False)
            comp = ", ".join(f"{o} {v}" for o, v in top[top > 0].head(3).items())
            ax.set_title(f"T{c}: {m.sum()} genes ({comp})", fontsize=8)
            ax.set_xticks([NBIN / 2, 1.5 * NBIN + 1, 2.5 * NBIN + 2]); ax.set_xticklabels(LIN)
        for ax in list(axs.flat)[k:]:
            ax.axis("off")
        fig.suptitle(f"D15 ISL1 targets among C1-C10 ({len(gk)} genes), k = {k}: mean z per new cluster "
                     "(solid control, dashed KO; title: largest old clusters)", fontsize=9)
        fig.tight_layout(); fig.savefig(od / f"curves_k{k}.png", dpi=110); plt.close(fig)
    print("wrote", od)


if __name__ == "__main__":
    main()
