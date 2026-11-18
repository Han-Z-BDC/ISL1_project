"""Step 7.1b: cluster one day (after 7.1_day_merged), Leiden resolution sweep.

Decisions (user): cluster day by day, starting with D15, on the per-day merge from 7.1 (no batch correction;
damaged cells removed) (2026-11-17). The graph is built from SCT B (2026-11-18): SCTransform v2 on the merged
day, S-G2M regressed, 3,000 variable features without MT- genes (7.1c, sct/SCT_B_pca.csv.gz); 30 PCs,
15 neighbours, seed 1234. The UMAP is 7.1c's SCT B UMAP (sct/sct_umaps.csv.gz), made from the same graph.
Leiden: igraph, n_iterations 2, undirected, seed 1234, at RESOLUTIONS (6.1 / earlier 7.1 convention).
Markers: Wilcoxon on the unregressed log1p data (X), each cluster vs the rest.

Writes ${OUT_DIR}/day_analysis/{day}/clustering/ (reruns overwrite):
  clusters.csv               barcode -> cluster at every resolution;  n_clusters.csv
  lineage.csv                each cluster's main parent at the next-lower resolution (+ % of its cells)
  composition_r*.csv         per cluster: cells per sample, % of all ctr / KO cells, KO_over_ctr,
                             % S and % G2M by genotype, median QC
  markers_r*.csv             top 50 Wilcoxon markers per cluster (logFC, % in / out, adj p)
  umap_clusters.png          all resolutions;  umap_r*_large.png  one resolution, labelled, + genotype split
  composition_r*.png         % of each sample's cells per cluster
  dotplot_genes_r*.png       fixed gene list by cluster;  dotplot_markers_r*.png  top 3 markers per cluster
  feature_genes.png          UMAP feature plots of the fixed gene list (log1p, capped at the 99th percentile of
                             expressing cells);  feature_markers_r*.png  top 2 markers of each cluster (user, 2026-11-18)
  {day}_clusters.html
and adds leiden_r* to obs, obsm X_pca_SCT_B / X_umap_SCT_B and the SCT_B graph to {day}_merged.h5ad
(overwritten in place; X_pca / X_umap from the 7.1 log1p route are kept).
Run: sbatch analysis/scripts/7.1b_day_cluster.sb D15
"""
import base64
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
DAY = sys.argv[1]
DAY_DIR = OUT_DIR / "day_analysis" / DAY
RES_DIR = DAY_DIR / "clustering"
SEED = 1234
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
N_PCS, N_NEIGHBORS = 30, 15
GENES = ["PDX1", "PTF1A", "NKX6-1", "SOX9", "ONECUT1", "ARX", "GCG", "PAX4", "PAX6", "INS", "SLC18A1", "FEV",
         "LMX1A", "SST", "GHRL", "NKX2-2", "NEUROG3", "CHGA", "RFX6", "ISL1", "AFP", "COL1A2", "MKI67"]
SAMPLE_ORDER = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
SAMPLE_PAL = {"ctr_rep1": "#2f5d9e", "ctr_rep2": "#7fa7dd", "KO_rep1": "#b23b2a", "KO_rep2": "#eb9a84"}
GREY = "#e3e3e3"


def key(r):
    return f"leiden_r{str(r).replace('.', 'p')}"


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def palette(cats):
    cmap = plt.get_cmap("tab20")
    extra = plt.get_cmap("tab20b")
    return {c: matplotlib.colors.to_hex((cmap if i < 20 else extra)(i % 20)) for i, c in enumerate(cats)}


FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", ["#dcdcdc", "#fdbb84", "#e34a33", "#7f0000"])


def feature_grid(a, u, genes, path, title, labels=None, ncol=6):
    """UMAP coloured by log1p expression, one panel per gene; expressing cells drawn on top."""
    nrow = int(np.ceil(len(genes) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.7 * ncol, 3.5 * nrow), squeeze=False)
    for i, (ax, g) in enumerate(zip(axs.flat, genes)):
        x = a[:, g].X
        v = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.ravel(x)
        vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1
        o = np.argsort(v, kind="stable")
        sc_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=3, cmap=FEAT_CMAP, vmin=0, vmax=vmax, linewidths=0, rasterized=True)
        plt.colorbar(sc_, ax=ax, fraction=0.04, pad=0.01)
        ax.set_title(g if labels is None else f"{g}  ({labels[i]})", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in list(axs.flat)[len(genes):]:
        ax.axis("off")
    fig.suptitle(title); fig.tight_layout(); fig.savefig(path, dpi=90); plt.close(fig)


def draw_clusters(ax, u, cl, title, mask=None, labels=True):
    cats = sorted(cl.unique(), key=int)
    pal = palette(cats)
    if mask is None:
        mask = np.ones(len(cl), bool)
    ax.scatter(u[~mask, 0], u[~mask, 1], c=GREY, s=2.5, linewidths=0, rasterized=True)
    idx = np.random.default_rng(SEED).permutation(np.where(mask)[0])
    ax.scatter(u[idx, 0], u[idx, 1], c=cl.map(pal).to_numpy()[idx], s=3.5, linewidths=0, rasterized=True)
    if labels:
        for c in cats:
            xy = np.median(u[(cl == c).to_numpy()], axis=0)
            ax.text(xy[0], xy[1], c, fontsize=10, weight="bold", ha="center", va="center")
    ax.set_title(title, fontsize=11); ax.set_xticks([]); ax.set_yticks([])


def main():
    RES_DIR.mkdir(parents=True, exist_ok=True)
    for p in RES_DIR.rglob("*"):
        if p.is_file():
            p.unlink()
    h5 = DAY_DIR / f"{DAY}_merged.h5ad"
    a = ad.read_h5ad(h5)
    a.obs = a.obs[[c for c in a.obs.columns if not c.startswith("leiden_r")]]
    pca = pd.read_csv(DAY_DIR / "sct" / "SCT_B_pca.csv.gz", index_col=0)
    um = pd.read_csv(DAY_DIR / "sct" / "sct_umaps.csv.gz", index_col=0)
    assert (pca.index == a.obs_names).all() and (um.index == a.obs_names).all()
    a.obsm["X_pca_SCT_B"] = pca.to_numpy()
    a.obsm["X_umap_SCT_B"] = um[["SCT_B_UMAP1", "SCT_B_UMAP2"]].to_numpy()
    sc.pp.neighbors(a, use_rep="X_pca_SCT_B", n_pcs=N_PCS, n_neighbors=N_NEIGHBORS, random_state=SEED, key_added="SCT_B")
    u = a.obsm["X_umap_SCT_B"]
    for r in RESOLUTIONS:
        sc.tl.leiden(a, resolution=r, key_added=key(r), random_state=SEED, flavor="igraph",
                     n_iterations=2, directed=False, neighbors_key="SCT_B")
    n_cl = {r: a.obs[key(r)].nunique() for r in RESOLUTIONS}
    print(f"{DAY}: {a.n_obs} cells; clusters {n_cl}")
    a.obs[[key(r) for r in RESOLUTIONS]].to_csv(RES_DIR / "clusters.csv")
    pd.Series(n_cl, name="n_clusters").rename_axis("resolution").to_csv(RES_DIR / "n_clusters.csv")

    # lineage: main parent of each cluster at the next-lower resolution
    rows = []
    for lo, hi in zip(RESOLUTIONS[:-1], RESOLUTIONS[1:]):
        ct = pd.crosstab(a.obs[key(hi)], a.obs[key(lo)])
        for c in ct.index:
            par = ct.loc[c].idxmax()
            rows.append({"resolution": hi, "cluster": int(c), "parent_resolution": lo, "parent": int(par),
                         "pct_cells_from_parent": round(100 * ct.loc[c, par] / ct.loc[c].sum(), 1)})
    pd.DataFrame(rows).to_csv(RES_DIR / "lineage.csv", index=False)

    fig, axs = plt.subplots(2, 3, figsize=(20, 12))
    for ax, r in zip(axs.flat, RESOLUTIONS):
        draw_clusters(ax, u, a.obs[key(r)].astype(str), f"resolution {r}: {n_cl[r]} clusters")
    fig.suptitle(f"{DAY}: Leiden clusters on SCT B (no batch correction)")
    fig.tight_layout(); fig.savefig(RES_DIR / "umap_clusters.png", dpi=110); plt.close(fig)
    feature_grid(a, u, [g for g in GENES if g in a.var_names], RES_DIR / "feature_genes.png",
                 f"{DAY}: lineage marker genes (log1p)")

    geno = a.obs["genotype"].astype(str).to_numpy()
    n_geno = pd.Series(geno).value_counts()
    genes = [g for g in GENES if g in a.var_names]
    for r in RESOLUTIONS:
        k = key(r); cl = a.obs[k].astype(str)
        # large UMAP + genotype split
        fig, axs = plt.subplots(1, 3, figsize=(25, 8))
        draw_clusters(axs[0], u, cl, f"resolution {r}: {n_cl[r]} clusters")
        for ax, gt, lab in zip(axs[1:], ("ctr", "KO"), ("Control", "ISL1 KO")):
            draw_clusters(ax, u, cl, f"{lab} cells (other genotype grey)", mask=(geno == gt))
        fig.tight_layout(); fig.savefig(RES_DIR / f"umap_r{r}_large.png", dpi=110); plt.close(fig)

        # composition
        ct = pd.crosstab(a.obs[k], a.obs["sample"]).reindex(columns=SAMPLE_ORDER, fill_value=0)
        ct.index = ct.index.astype(int); ct = ct.sort_index()
        comp = ct.copy()
        comp["cells"] = ct.sum(axis=1)
        comp["pct_of_ctr_cells"] = (100 * ct[["ctr_rep1", "ctr_rep2"]].sum(axis=1) / n_geno["ctr"]).round(2)
        comp["pct_of_KO_cells"] = (100 * ct[["KO_rep1", "KO_rep2"]].sum(axis=1) / n_geno["KO"]).round(2)
        comp["KO_over_ctr"] = (comp["pct_of_KO_cells"] / comp["pct_of_ctr_cells"]).round(2)
        obs = a.obs.assign(cl=a.obs[k].astype(int))
        for gt in ("ctr", "KO"):
            sub = obs[obs["genotype"] == gt]
            ph = pd.crosstab(sub["cl"], sub["phase"]).reindex(columns=["G1", "S", "G2M"], fill_value=0)
            comp[f"{gt}_pct_S"] = (100 * ph["S"] / ph.sum(axis=1)).reindex(comp.index).round(1)
            comp[f"{gt}_pct_G2M"] = (100 * ph["G2M"] / ph.sum(axis=1)).reindex(comp.index).round(1)
        med = obs.groupby("cl")[["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
                                 "nuclear_fraction"]].median()
        med.columns = ["median_counts", "median_genes", "median_pct_mt", "median_pct_ribo", "median_NF"]
        comp = comp.join(med.round(2))
        comp.to_csv(RES_DIR / f"composition_r{r}.csv")
        pct = 100 * ct / ct.sum(axis=0)
        fig, ax = plt.subplots(figsize=(max(10, 0.6 * len(pct)), 4))
        x = np.arange(len(pct)); w = 0.2
        for i, s in enumerate(SAMPLE_ORDER):
            ax.bar(x + (i - 1.5) * w, pct[s], width=w, color=SAMPLE_PAL[s], label=s)
        ax.set_xticks(x); ax.set_xticklabels(pct.index); ax.set_ylabel("% of sample's cells"); ax.set_xlabel("cluster")
        ax.legend(fontsize=8, frameon=False, ncol=4); ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"{DAY} resolution {r}: % of each sample's cells per cluster", fontsize=10)
        fig.tight_layout(); fig.savefig(RES_DIR / f"composition_r{r}.png", dpi=110); plt.close(fig)

        # markers (unregressed log1p)
        a.obs[k] = pd.Categorical(a.obs[k], sorted(a.obs[k].unique(), key=int))
        sc.tl.rank_genes_groups(a, k, method="wilcoxon", n_genes=50, use_raw=False, pts=True)
        mk = sc.get.rank_genes_groups_df(a, None)
        mk = mk.groupby("group", observed=True).head(50)
        mk.to_csv(RES_DIR / f"markers_r{r}.csv", index=False)
        dp = sc.pl.dotplot(a, genes, groupby=k, show=False, return_fig=True)
        dp.savefig(RES_DIR / f"dotplot_genes_r{r}.png", dpi=110); plt.close("all")
        top = []
        for grp, g in mk[mk["logfoldchanges"] > 0].groupby("group", observed=True, sort=False):
            for nm in g["names"].head(3):
                if nm not in top:
                    top.append(nm)
        dp = sc.pl.dotplot(a, top, groupby=k, show=False, return_fig=True, standard_scale="var")
        dp.savefig(RES_DIR / f"dotplot_markers_r{r}.png", dpi=110); plt.close("all")
        a.uns.pop("rank_genes_groups", None)
        fm, fl = [], []
        for grp, g in mk[mk["logfoldchanges"] > 0].groupby("group", observed=True, sort=False):
            for nm in g["names"].head(2):
                if nm not in fm:
                    fm.append(nm); fl.append(f"cl {grp}")
        feature_grid(a, u, fm, RES_DIR / f"feature_markers_r{r}.png",
                     f"{DAY} resolution {r}: top 2 markers of each cluster (log1p)", labels=fl)

    a.write_h5ad(h5, compression="gzip")

    # ---- HTML
    def table(df):
        return df.to_html(border=0, classes="t")
    secs = ""
    for r in RESOLUTIONS:
        comp = pd.read_csv(RES_DIR / f"composition_r{r}.csv", index_col=0)
        mk = pd.read_csv(RES_DIR / f"markers_r{r}.csv")
        top10 = (mk[mk["logfoldchanges"] > 0].groupby("group").head(10).groupby("group")["names"]
                 .apply(lambda v: ", ".join(v)).rename("top 10 markers (up)"))
        top10.index = top10.index.astype(int)
        secs += f"""<h2 id="r{r}">Resolution {r}: {n_cl[r]} clusters</h2>
<img src="data:image/png;base64,{b64(RES_DIR / f'umap_r{r}_large.png')}" alt="clusters at resolution {r}, with genotype split">
<img src="data:image/png;base64,{b64(RES_DIR / f'composition_r{r}.png')}" alt="composition at resolution {r}">
<div class=tw>{table(comp)}</div><div class=tw>{table(top10.to_frame())}</div>
<img src="data:image/png;base64,{b64(RES_DIR / f'dotplot_genes_r{r}.png')}" alt="gene dotplot at resolution {r}">
<img src="data:image/png;base64,{b64(RES_DIR / f'dotplot_markers_r{r}.png')}" alt="top marker dotplot at resolution {r}">
<h3>Feature plots: top 2 markers of each cluster</h3>
<img src="data:image/png;base64,{b64(RES_DIR / f'feature_markers_r{r}.png')}" alt="UMAP feature plots of the top markers at resolution {r}">"""
    nav = " · ".join(f'<a href="#r{r}">{r}</a>' for r in RESOLUTIONS)
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{DAY} Clusters</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1240px;margin:0 auto;padding:24px 16px 60px}} a{{color:inherit}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:34px 0 8px}} h3{{font-size:16px;margin:20px 0 6px}}
.sub{{color:var(--mut)}} .tw{{overflow-x:auto;margin:8px 0}}
table.t{{border-collapse:collapse;font-size:12.5px}} .t th,.t td{{padding:3px 7px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
.t td:last-child{{text-align:left}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>{DAY}: Leiden clusters</h1>
<p class=sub>{a.n_obs:,} cells; per-day merge without batch correction (7.1). Graph and UMAP from SCT B: SCTransform v2
on the merged day, S&minus;G2M regressed, 3,000 variable features without MT- genes; 30 PCs, 15 neighbours.
Leiden (igraph, 2 iterations, seed {SEED}). Markers: Wilcoxon, cluster vs rest, unregressed log1p. KO_over_ctr = % of
KO cells in the cluster / % of ctr cells. Resolutions: {nav}</p>
<img src="data:image/png;base64,{b64(RES_DIR / 'umap_clusters.png')}" alt="clusters at all resolutions">
<h2>Feature plots: lineage marker genes</h2>
<p class=sub>log1p expression; colour capped at the 99th percentile of expressing cells; expressing cells drawn on top.</p>
<img src="data:image/png;base64,{b64(RES_DIR / 'feature_genes.png')}" alt="UMAP feature plots of lineage marker genes">
{secs}
<p class=sub>Full marker tables (top 50 per cluster) and lineage.csv (parent cluster at the next-lower resolution) are in this folder.</p>
</main></body></html>"""
    (RES_DIR / f"{DAY}_clusters.html").write_text(html)
    print("wrote", RES_DIR)


if __name__ == "__main__":
    main()
