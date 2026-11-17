"""Step 7.1: one time point on its own -- merge its 4 samples, no batch correction.

Decision (user, 2026-11-17): look at the data day by day, starting with D15; merge only, no batch
correction of any kind (all 4 samples of a day share one differentiation batch and OCM channel).

Inputs (current QC: 4 MAD + min 2,000 genes + 10% MT, SoupX, scDblFinder only):
  ${OUT_DIR}/detect_doublets/h5ad/{sample}.doublets_removed.h5ad for the day's 4 samples
  (X = SoupX-corrected counts; layers counts / soupx_counts).
Processing (pipeline 5.3 log path + 6.1 graph conventions, on this day's cells only):
  concat -> normalize_total (1e4) -> log1p -> 2,000 HVGs (flavor seurat) -> PCA (50 comps)
  -> neighbours on 30 PCs (n_neighbors 15) -> UMAP -> Leiden (igraph, n_iterations 2) at RESOLUTIONS.
  Seed 1234 throughout.
Writes ${OUT_DIR}/day_analysis/{day}/ (reruns overwrite):
  {day}_merged.h5ad        cells x genes; X = log1p-normalized, layers counts / soupx_counts,
                           obsm X_pca + X_umap, obs leiden_r* for every resolution
  per_sample.csv           cells and QC medians per sample
  clusters.csv             barcode -> cluster at every resolution;  n_clusters.csv
  composition_r*.csv       per cluster: cells per sample, % of ctr and of KO cells
  markers_r*.csv           top 10 Wilcoxon markers per cluster
  umap_samples.png, umap_qc.png, umap_clusters.png, composition.png, dotplot_r*.png,
  {day}_summary.html
Run: sbatch analysis/scripts/7.1_day_merged.sb D15
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

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = Path(os.environ["OUT_DIR"])
DAY = sys.argv[1]
RES_DIR = OUT_DIR / "day_analysis" / DAY
SEED = 1234
N_HVG, N_PCS_STORE, N_PCS, N_NEIGHBORS = 2000, 50, 30, 15
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
GENES = ["PDX1", "PTF1A", "NKX6-1", "ARX", "GCG", "PAX4", "PAX6", "INS", "SLC18A1", "FEV", "SST",
         "NKX2-2", "RFX6", "ISL1"]
SAMPLE_ORDER = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
SAMPLE_PAL = {"ctr_rep1": "#2f5d9e", "ctr_rep2": "#7fa7dd", "KO_rep1": "#b23b2a", "KO_rep2": "#eb9a84"}
GENO_PAL = {"ctr": "#3b6fb6", "KO": "#d0573f"}


def key(res):
    return f"leiden_r{str(res).replace('.', 'p')}"


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def load_day():
    rename = pd.read_csv(REPO / "study_config" / "rename_manifest.csv")
    rows = rename[rename["new_name"].str.contains(f"_{DAY}-")]
    assert len(rows) == 4, rows
    parts = []
    for old, new in zip(rows["old_name"], rows["new_name"]):
        a = ad.read_h5ad(OUT_DIR / "detect_doublets" / "h5ad" / f"{old}.doublets_removed.h5ad")
        tag = new.split("-", 1)[1]                      # ctr_rep1 ...
        a.obs["sample"] = tag
        a.obs["genotype"] = tag.split("_")[0]
        a.obs["replicate"] = tag.split("_")[1]
        a.obs_names = [f"{new}#{b}" for b in a.obs_names]
        keep_obs = ["sample", "genotype", "replicate", "total_counts", "n_genes_by_counts", "pct_counts_mt",
                    "pct_counts_ribo", "scDblFinder_doublet_score"]
        a.obs = a.obs[keep_obs]
        a.obsm, a.uns, a.obsp = {}, {}, {}
        parts.append(a)
    var = parts[0].var[["gene_ids", "feature_types"]].copy()
    for p in parts:
        assert (p.var_names == var.index).all()
    m = ad.concat(parts, join="inner", merge=None)
    m.var = var
    m.obs["sample"] = pd.Categorical(m.obs["sample"], SAMPLE_ORDER)
    m.obs["genotype"] = pd.Categorical(m.obs["genotype"], ["ctr", "KO"])
    return m


def scatter(ax, u, values, title, pal=None, cmap="viridis", vmax=None, order=None):
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    if pal is not None:
        idx = np.random.default_rng(SEED).permutation(len(values))   # no category drawn last on top
        cols = pd.Series(values).astype(str).map(pal).to_numpy()
        ax.scatter(u[idx, 0], u[idx, 1], c=cols[idx], s=1.5, linewidths=0, rasterized=True)
        for k, c in pal.items():
            ax.scatter([], [], c=c, s=20, label=k)
        ax.legend(fontsize=7, frameon=False, markerscale=1, loc="best")
    else:
        v = np.asarray(values, float)
        o = np.argsort(v)
        sc_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=1.5, cmap=cmap, vmax=vmax, linewidths=0, rasterized=True)
        plt.colorbar(sc_, ax=ax, fraction=0.04)


def main():
    RES_DIR.mkdir(parents=True, exist_ok=True)
    for p in RES_DIR.iterdir():                       # rerun replaces the day's outputs
        if p.is_file():
            p.unlink()
        else:
            for q in p.rglob("*"):
                if q.is_file():
                    q.unlink()
    a = load_day()
    print(f"{DAY}: {a.n_obs} cells, {a.n_vars} genes")

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat")
    hv = a[:, a.var["highly_variable"]].copy()
    sc.tl.pca(hv, n_comps=N_PCS_STORE, random_state=SEED)
    a.obsm["X_pca"] = hv.obsm["X_pca"]
    a.uns["pca_variance_ratio"] = hv.uns["pca"]["variance_ratio"]
    sc.pp.neighbors(a, use_rep="X_pca", n_pcs=N_PCS, n_neighbors=N_NEIGHBORS, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    for r in RESOLUTIONS:
        sc.tl.leiden(a, resolution=r, key_added=key(r), random_state=SEED, flavor="igraph",
                     n_iterations=2, directed=False)
    u = a.obsm["X_umap"]

    # ---- per-sample table
    ps = a.obs.groupby("sample", observed=True).agg(
        cells=("total_counts", "size"), median_counts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"), median_pct_mt=("pct_counts_mt", "median"),
        median_pct_ribo=("pct_counts_ribo", "median")).round(1)
    ps.to_csv(RES_DIR / "per_sample.csv")

    # ---- UMAPs: samples / genotype
    fig, axs = plt.subplots(2, 3, figsize=(16, 10))
    scatter(axs[0, 0], u, a.obs["sample"], "sample", pal=SAMPLE_PAL)
    scatter(axs[0, 1], u, a.obs["genotype"], "genotype", pal=GENO_PAL)
    for ax, g in zip([axs[0, 2], axs[1, 0]], ["ctr", "KO"]):
        ax.scatter(u[:, 0], u[:, 1], c="#dddddd", s=1, linewidths=0, rasterized=True)
        m = (a.obs["genotype"] == g).to_numpy()
        sub = a.obs.loc[m, "sample"]
        scatter(ax, u[m], sub, f"{g} only", pal={k: v for k, v in SAMPLE_PAL.items() if k.startswith(g)})
    scatter(axs[1, 1], u, a.obs["replicate"], "replicate", pal={"rep1": "#555555", "rep2": "#e0a030"})
    axs[1, 2].axis("off")
    fig.suptitle(f"{DAY}: 4 samples merged, no batch correction ({a.n_obs:,} cells)")
    fig.tight_layout(); fig.savefig(RES_DIR / "umap_samples.png", dpi=110); plt.close(fig)

    # ---- UMAPs: QC
    fig, axs = plt.subplots(2, 3, figsize=(16, 10))
    for ax, (col, lab, vmax) in zip(axs.flat, [("total_counts", "log10 counts", None), ("n_genes_by_counts", "log10 genes", None),
                                             ("pct_counts_mt", "% MT", 10), ("pct_counts_ribo", "% ribo", None),
                                             ("scDblFinder_doublet_score", "scDblFinder score", None)]):
        v = a.obs[col].to_numpy(float)
        if lab.startswith("log10"):
            v = np.log10(v)
        scatter(ax, u, v, lab, vmax=vmax)
    axs[1, 2].axis("off")
    fig.tight_layout(); fig.savefig(RES_DIR / "umap_qc.png", dpi=110); plt.close(fig)

    # ---- UMAPs: clusters at every resolution
    n_cl = {r: a.obs[key(r)].nunique() for r in RESOLUTIONS}
    pd.Series(n_cl, name="n_clusters").rename_axis("resolution").to_csv(RES_DIR / "n_clusters.csv")
    a.obs[[key(r) for r in RESOLUTIONS]].to_csv(RES_DIR / "clusters.csv")
    fig, axs = plt.subplots(2, 3, figsize=(18, 11))
    for ax, r in zip(axs.flat, RESOLUTIONS):
        cl = a.obs[key(r)].astype(str)
        cats = sorted(cl.unique(), key=int)
        cmap = plt.get_cmap("tab20" if len(cats) > 10 else "tab10")
        pal = {c: matplotlib.colors.to_hex(cmap(i % cmap.N)) for i, c in enumerate(cats)}
        idx = np.random.default_rng(SEED).permutation(len(cl))
        ax.scatter(u[idx, 0], u[idx, 1], c=cl.map(pal).to_numpy()[idx], s=1.5, linewidths=0, rasterized=True)
        for c in cats:
            xy = np.median(u[(cl == c).to_numpy()], axis=0)
            ax.text(xy[0], xy[1], c, fontsize=9, weight="bold", ha="center", va="center")
        ax.set_title(f"resolution {r}: {n_cl[r]} clusters", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(RES_DIR / "umap_clusters.png", dpi=110); plt.close(fig)

    # ---- composition + markers + dotplots per resolution
    n_geno = a.obs["genotype"].value_counts()
    fig, axs = plt.subplots(len(RESOLUTIONS), 1, figsize=(14, 3.2 * len(RESOLUTIONS)))
    genes = [g for g in GENES if g in a.var_names]
    for ax, r in zip(axs, RESOLUTIONS):
        k = key(r)
        ct = pd.crosstab(a.obs[k], a.obs["sample"]).reindex(columns=SAMPLE_ORDER, fill_value=0)
        ct.index = ct.index.astype(int); ct = ct.sort_index()
        comp = ct.copy()
        comp["cells"] = ct.sum(axis=1)
        comp["pct_of_ctr_cells"] = (100 * ct[["ctr_rep1", "ctr_rep2"]].sum(axis=1) / n_geno["ctr"]).round(2)
        comp["pct_of_KO_cells"] = (100 * ct[["KO_rep1", "KO_rep2"]].sum(axis=1) / n_geno["KO"]).round(2)
        comp["KO_over_ctr"] = (comp["pct_of_KO_cells"] / comp["pct_of_ctr_cells"]).round(2)
        comp.to_csv(RES_DIR / f"composition_r{r}.csv")
        # bars: % of each sample's cells in each cluster (replicates side by side)
        pct = 100 * ct / ct.sum(axis=0)
        x = np.arange(len(pct)); w = 0.2
        for i, s in enumerate(SAMPLE_ORDER):
            ax.bar(x + (i - 1.5) * w, pct[s], width=w, color=SAMPLE_PAL[s], label=s)
        ax.set_xticks(x); ax.set_xticklabels(pct.index); ax.set_ylabel("% of sample's cells")
        ax.set_title(f"resolution {r}", fontsize=10); ax.spines[["top", "right"]].set_visible(False)
        axs[0].legend(fontsize=8, frameon=False, ncol=4)

        sc.tl.rank_genes_groups(a, k, method="wilcoxon", n_genes=10, use_raw=False)
        mk = sc.get.rank_genes_groups_df(a, None)
        mk = mk.groupby("group", observed=True).head(10)
        mk.to_csv(RES_DIR / f"markers_r{r}.csv", index=False)
        dp = sc.pl.dotplot(a, genes, groupby=k, show=False, return_fig=True, standard_scale=None)
        dp.savefig(RES_DIR / f"dotplot_r{r}.png", dpi=110); plt.close("all")
    fig.tight_layout(); fig.savefig(RES_DIR / "composition.png", dpi=110); plt.close(fig)
    a.uns.pop("rank_genes_groups", None)

    a.write_h5ad(RES_DIR / f"{DAY}_merged.h5ad", compression="gzip")

    # ---- HTML
    def table(df, idx=True):
        return df.to_html(index=idx, border=0, classes="t")
    comp_tabs = "".join(
        f"<h3>Resolution {r} ({n_cl[r]} clusters)</h3><div class=tw>{table(pd.read_csv(RES_DIR / f'composition_r{r}.csv', index_col=0))}</div>"
        f"<img src='data:image/png;base64,{b64(RES_DIR / f'dotplot_r{r}.png')}' alt='dotplot resolution {r}'>"
        for r in RESOLUTIONS)
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{DAY} Merged Summary</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1180px;margin:0 auto;padding:24px 16px 60px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:30px 0 8px}} h3{{font-size:15px;margin:22px 0 6px}}
.sub{{color:var(--mut)}} .tw{{overflow-x:auto}}
table.t{{border-collapse:collapse;font-size:13px}} .t th,.t td{{padding:3px 8px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>{DAY}: four samples merged, no batch correction</h1>
<p class=sub>{a.n_obs:,} cells after QC (4 MAD + at least 2,000 genes + %MT &le; 10, SoupX, scDblFinder).
normalize_total (10,000) &rarr; log1p &rarr; 2,000 HVGs &rarr; PCA &rarr; 30 PCs, 15 neighbours &rarr; UMAP, Leiden. Seed 1234.</p>
<h2>Cells per sample</h2><div class=tw>{table(ps)}</div>
<h2>Samples on the UMAP</h2><img src="data:image/png;base64,{b64(RES_DIR / 'umap_samples.png')}" alt="UMAP by sample and genotype">
<h2>QC metrics on the UMAP</h2><img src="data:image/png;base64,{b64(RES_DIR / 'umap_qc.png')}" alt="UMAP coloured by QC metrics">
<h2>Clusters</h2><img src="data:image/png;base64,{b64(RES_DIR / 'umap_clusters.png')}" alt="Leiden clusters at six resolutions">
<h2>Cluster composition</h2><p class=sub>Bars: % of each sample's cells in each cluster. Tables: cells per sample;
% of all ctr / KO cells in the cluster; KO_over_ctr = ratio of those two.</p>
<img src="data:image/png;base64,{b64(RES_DIR / 'composition.png')}" alt="cluster composition per sample">
{comp_tabs}
<p class=sub>Top 10 Wilcoxon markers per cluster: markers_r*.csv in this folder.</p>
</main></body></html>"""
    (RES_DIR / f"{DAY}_summary.html").write_text(html)
    print("clusters:", n_cl)
    print("wrote", RES_DIR)


if __name__ == "__main__":
    main()
