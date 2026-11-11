"""Step 7.1: one time point on its own, no batch correction -- summary + genes split by genotype.

Decision (user, 2026-11-09): analyse each day separately, no integration (5.6 showed replicates
already mix, and per-sample correction pulls ctr and KO together). Starts with D15.
No clustering here (not requested yet).

Inputs:
  - ${OUT_DIR}/integration_per_day/{day}.h5ad (5.6): obs and the three UMAPs from 5.6, used as is:
    X_umap_none (log-path PCA on this day's cells only), X_umap_harmony (Harmony on sample_id),
    X_umap_scvi (scVI, batch_key=sample_id). Same genes and colour scales on each, for comparison.
  - ${OUT_DIR}/integration/merged_rna_raw.h5ad (5.2): soupx_counts for those cells ->
    normalize_total(1e4) -> log1p (same as 6.1).

Writes ${OUT_DIR}/day_analysis/{day}/:
  per_sample.csv, genes_by_genotype.csv, umap_genotype.png (rows = method),
  genes_by_genotype_{method}.png, dotplot_genes.png, {day}_summary.html
Clustering (user, 2026-11-09): for each method, neighbours on that method's embedding (X_none /
X_harmony / X_scvi; scanpy defaults, the same graph 5.6 built its UMAP from) -> Leiden (igraph,
n_iterations=2) at RESOLUTIONS. Writes clusters_{day}.csv (barcode x method_resolution),
n_clusters.csv, cluster_ari.csv (agreement between methods), cluster_composition.csv,
clusters_{method}.png (UMAPs + % of ctr / KO cells per cluster) and
dotplot_{method}_r{res}.png (requested genes by cluster).
Run: sbatch analysis/scripts/7.1_day_summary.sb D15 --cluster   (D22: without --cluster)
"""
import argparse
import base64
import os
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
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

OUT_DIR = Path(os.environ["OUT_DIR"])
GENES = ["PDX1", "PTF1A", "NKX6-1", "ARX", "GCG", "PAX4", "PAX6", "INS", "SLC18A1", "FEV", "SST",
         "NKX2-2", "RFX6"]
EXTRA = ["ISL1"]                     # table / dot plot only: check of the KO
GENO_PAL = {"ctr": "#3b6fb6", "KO": "#d0573f"}
SAMPLE_ORDER = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
SEED = 1234
METHODS = {"none": "no integration", "harmony": "Harmony (per sample)", "scvi": "scVI (per sample)"}
SAMPLE_PAL = {"ctr_rep1": "#2f5d9e", "ctr_rep2": "#7fa7dd", "KO_rep1": "#b23b2a", "KO_rep2": "#eb9a84"}


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("day")
    ap.add_argument("--cluster", action="store_true",
                    help="also run the Leiden sweep on each method (D15: requested 2026-11-09)")
    args = ap.parse_args()
    day = args.day
    out = OUT_DIR / "day_analysis" / day
    out.mkdir(parents=True, exist_ok=True)

    d = sc.read_h5ad(OUT_DIR / "integration_per_day" / f"{day}.h5ad")
    raw = ad.read_h5ad(OUT_DIR / "integration" / "merged_rna_raw.h5ad", backed="r")
    lay = "soupx_counts" if "soupx_counts" in raw.layers else "counts"
    idx = raw.obs_names.get_indexer(d.obs_names)
    assert (idx >= 0).all(), "cells missing from merged_rna_raw"
    sub = raw[np.sort(idx)].to_memory()
    sub = sub[d.obs_names]
    a = ad.AnnData(X=sub.layers[lay].copy(), obs=d.obs.copy(), var=sub.var[[]].copy())
    del sub, raw
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    umaps = {m: np.asarray(d.obsm[f"X_umap_{m}"]) for m in METHODS}
    print(f"{day}: {a.shape}, counts layer {lay}", flush=True)
    geno = a.obs["genotype"].astype(str).to_numpy()
    samp = a.obs["sample_short"].astype(str).to_numpy()

    # ---- tables ------------------------------------------------------------------------------
    per_sample = a.obs.groupby("sample_short", observed=True).agg(
        cells=("total_counts", "size"), median_UMIs=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"), median_pct_mt=("pct_counts_mt", "median")
    ).reindex(SAMPLE_ORDER).round(1)
    per_sample.to_csv(out / "per_sample.csv")

    genes = [g for g in GENES if g in a.var_names]
    missing = [g for g in GENES if g not in a.var_names]
    tg = genes + [g for g in EXTRA if g in a.var_names]
    X = a[:, tg].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    rows = []
    for j, g in enumerate(tg):
        r = {"gene": g}
        for grp, msk in [("ctr", geno == "ctr"), ("KO", geno == "KO")] + [(s, samp == s) for s in SAMPLE_ORDER]:
            r[f"pct_expr_{grp}"] = round(100 * float((X[msk, j] > 0).mean()), 1)
        for grp in ("ctr", "KO"):
            r[f"mean_log_{grp}"] = round(float(X[geno == grp, j].mean()), 3)
        rows.append(r)
    gtab = pd.DataFrame(rows)
    gtab.to_csv(out / "genes_by_genotype.csv", index=False)
    print(per_sample.to_string(), gtab.to_string(), sep="\n\n", flush=True)

    # ---- genotype UMAP: rows = method -------------------------------------------------------
    order = np.random.default_rng(0).permutation(a.n_obs)
    fig, axes = plt.subplots(len(METHODS), 3, figsize=(17, 5.4 * len(METHODS)))
    for i, (meth, mname) in enumerate(METHODS.items()):
        u = umaps[meth]
        ax = axes[i, 0]
        for g, c in GENO_PAL.items():
            ax.scatter([], [], c=c, s=30, label=f"{g} (n={(geno == g).sum():,})")
        ax.scatter(u[order, 0], u[order, 1], s=1.5, c=[GENO_PAL[x] for x in geno[order]], linewidths=0)
        ax.legend(frameon=False, fontsize=9); ax.set_title(f"{mname}: genotype", fontsize=11)
        for ax, g in zip(axes[i, 1:], ("ctr", "KO")):
            m = geno == g
            ax.scatter(u[~m, 0], u[~m, 1], s=1.5, c="#e2e2e2", linewidths=0)
            for s_ in [x for x in SAMPLE_ORDER if x.startswith(g)]:
                mm = samp == s_
                ax.scatter(u[mm, 0], u[mm, 1], s=1.5, c=SAMPLE_PAL[s_], linewidths=0)
                ax.scatter([], [], c=SAMPLE_PAL[s_], s=30, label=f"{s_} (n={mm.sum():,})")
            ax.legend(frameon=False, fontsize=9)
            ax.set_title(f"{mname}: {g} only, by replicate (other genotype grey)", fontsize=11)
        for ax in axes[i]:
            clean(ax)
    fig.tight_layout()
    fig.savefig(out / "umap_genotype.png", dpi=80, bbox_inches="tight"); plt.close(fig)

    # ---- genes: per method, per gene a ctr and a KO panel, same colour scale, other genotype grey
    vmaxes = []
    for j in range(len(genes)):
        pos = X[X[:, j] > 0, j]
        vmaxes.append(float(np.percentile(pos, 99)) if pos.size else 1.0)
    ncol_pairs = 2
    nrow = int(np.ceil(len(genes) / ncol_pairs))
    for meth, mname in METHODS.items():
        u = umaps[meth]
        fig, axes = plt.subplots(nrow, 2 * ncol_pairs, figsize=(5.0 * 2 * ncol_pairs, 4.4 * nrow),
                                 sharex=True, sharey=True)
        axes = np.atleast_2d(axes)
        for j, g in enumerate(genes):
            x = X[:, j]
            rr, cc = divmod(j, ncol_pairs)
            for t, grp in enumerate(("ctr", "KO")):
                ax = axes[rr, 2 * cc + t]
                m = geno == grp
                ax.scatter(u[~m, 0], u[~m, 1], s=1, c="#efefef", linewidths=0)
                ii = np.where(m)[0]
                ii = ii[np.argsort(x[ii], kind="stable")]
                sca = ax.scatter(u[ii, 0], u[ii, 1], s=1.6, c=x[ii], cmap="viridis", vmin=0, vmax=vmaxes[j],
                                 linewidths=0)
                ax.set_title(f"{g}  |  {grp}  ({100 * (x[m] > 0).mean():.1f}% expr)", fontsize=11,
                             color=GENO_PAL[grp])
                clean(ax)
                if t == 1:
                    fig.colorbar(sca, ax=ax, fraction=0.04, pad=0.01, label="log1p norm")
        for j in range(len(genes), nrow * ncol_pairs):
            rr, cc = divmod(j, ncol_pairs)
            axes[rr, 2 * cc].axis("off"); axes[rr, 2 * cc + 1].axis("off")
        fig.suptitle(f"{day}, {mname} UMAP: expression by genotype (colour scale shared within each gene, "
                     f"capped at the 99th percentile of expressing cells)", fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.savefig(out / f"genes_by_genotype_{meth}.png", dpi=70, bbox_inches="tight"); plt.close(fig)

    a.obs["sample_short"] = pd.Categorical(a.obs["sample_short"].astype(str), categories=SAMPLE_ORDER)
    dp = sc.pl.dotplot(a, tg, groupby="sample_short", show=False, return_fig=True, title=f"{day}: by sample")
    dp.savefig(out / "dotplot_genes.png", dpi=110, bbox_inches="tight")
    plt.close("all")

    if args.cluster:
        # ---- clustering on each method's embedding ------------------------------------------------
        def rk(meth, r):
            return f"{meth}_r{str(r).replace('.', 'p')}"

        labels, ncl, comp_rows = {}, [], []
        for meth in METHODS:
            g_ = ad.AnnData(X=np.empty((a.n_obs, 0), dtype=np.float32), obs=pd.DataFrame(index=a.obs_names))
            g_.obsm["X_emb"] = np.asarray(d.obsm[f"X_{meth}"])
            sc.pp.neighbors(g_, use_rep="X_emb", random_state=SEED)
            for r in RESOLUTIONS:
                sc.tl.leiden(g_, resolution=r, key_added="l", random_state=SEED, flavor="igraph",
                             n_iterations=2, directed=False)
                lab = g_.obs["l"].astype(str).to_numpy()
                labels[rk(meth, r)] = lab
                ncl.append({"method": meth, "resolution": r, "n_clusters": len(set(lab))})
                ct = pd.crosstab(lab, samp).reindex(columns=SAMPLE_ORDER, fill_value=0)
                for c in sorted(ct.index, key=int):
                    row = {"method": meth, "resolution": r, "cluster": int(c), "cells": int(ct.loc[c].sum())}
                    for s_ in SAMPLE_ORDER:
                        row[f"%_{s_}"] = round(100 * ct.loc[c, s_] / ct[s_].sum(), 2)
                    for g in ("ctr", "KO"):
                        cols = [x for x in SAMPLE_ORDER if x.startswith(g)]
                        row[f"%_of_{g}"] = round(100 * ct.loc[c, cols].sum() / ct[cols].sum().sum(), 2)
                    row["KO_share_norm"] = round(100 * row["%_of_KO"] / max(row["%_of_KO"] + row["%_of_ctr"], 1e-9), 1)
                    comp_rows.append(row)
                print(f"  {meth} r={r}: {len(set(lab))} clusters", flush=True)
        pd.DataFrame({"barcode": a.obs_names, **labels}).to_csv(out / f"clusters_{day}.csv", index=False)
        ncl = pd.DataFrame(ncl)
        ncl_w = ncl.pivot(index="method", columns="resolution", values="n_clusters").reindex(list(METHODS))
        ncl_w.index = [METHODS[m] for m in ncl_w.index]
        ncl.to_csv(out / "n_clusters.csv", index=False)
        comp = pd.DataFrame(comp_rows)
        comp.to_csv(out / "cluster_composition.csv", index=False)
        pairs = [("none", "harmony"), ("none", "scvi"), ("harmony", "scvi")]
        ari = pd.DataFrame({f"{METHODS[x]} vs {METHODS[y]}": [round(adjusted_rand_score(labels[rk(x, r)], labels[rk(y, r)]), 2)
                                                             for r in RESOLUTIONS] for x, y in pairs},
                           index=pd.Index(RESOLUTIONS, name="resolution")).T
        ari.to_csv(out / "cluster_ari.csv")
        print(ncl_w.to_string(), ari.to_string(), sep="\n\n", flush=True)

        cmap = plt.get_cmap("tab20")
        def pal(i):
            return cmap(i % 20) if i < 20 else plt.get_cmap("tab20b")(i % 20)

        for meth, mname in METHODS.items():
            u = umaps[meth]
            nR = len(RESOLUTIONS)
            fig, axes = plt.subplots(2, nR, figsize=(5.0 * nR, 9.6), gridspec_kw={"height_ratios": [1.15, 0.85]})
            for t, r in enumerate(RESOLUTIONS):
                lab = labels[rk(meth, r)]
                cats = sorted(set(lab), key=int)
                ax = axes[0, t]
                for i, c in enumerate(cats):
                    m = lab == c
                    ax.scatter(u[m, 0], u[m, 1], s=1.2, color=pal(i), linewidths=0)
                for c in cats:
                    m = lab == c
                    ax.text(np.median(u[m, 0]), np.median(u[m, 1]), c, fontsize=8, weight="bold", ha="center",
                            va="center", bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
                ax.set_title(f"r={r}  ({len(cats)} clusters)", fontsize=12)
                clean(ax)
                ax = axes[1, t]
                cc = comp[(comp.method == meth) & (comp.resolution == r)].sort_values("cluster")
                x = np.arange(len(cc))
                ax.bar(x - 0.2, cc["%_of_ctr"], width=0.4, color=GENO_PAL["ctr"], label="% of ctr cells")
                ax.bar(x + 0.2, cc["%_of_KO"], width=0.4, color=GENO_PAL["KO"], label="% of KO cells")
                for s_, off in (("ctr_rep1", -0.2), ("ctr_rep2", -0.2), ("KO_rep1", 0.2), ("KO_rep2", 0.2)):
                    ax.scatter(x + off, cc[f"%_{s_}"], s=6, color="black", zorder=3, linewidths=0)
                ax.set_xticks(x); ax.set_xticklabels(cc["cluster"], fontsize=7 if len(cc) > 15 else 9)
                ax.spines[["top", "right"]].set_visible(False)
                if t == 0:
                    ax.set_ylabel("% of genotype's cells (dots = replicates)")
                    ax.legend(frameon=False, fontsize=8)
            fig.suptitle(f"{day}, {mname}: Leiden clusters on this method's embedding", fontsize=14)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            fig.savefig(out / f"clusters_{meth}.png", dpi=70, bbox_inches="tight"); plt.close(fig)

            for r in RESOLUTIONS:
                lab = labels[rk(meth, r)]
                cats = sorted(set(lab), key=int)
                a.obs["_cl"] = pd.Categorical(lab, categories=cats)
                dp = sc.pl.dotplot(a, tg, groupby="_cl", show=False, return_fig=True,
                                   title=f"{mname}, r={r}")
                dp.savefig(out / f"dotplot_{meth}_r{r}.png", dpi=90, bbox_inches="tight")
                plt.close("all")
        a.obs.drop(columns="_cl", inplace=True)

    # ---- html --------------------------------------------------------------------------------
    def table(df, index=True):
        return df.to_html(index=index, border=0, classes="t", na_rep="")

    cluster_html, cluster_link = "", ""
    if args.cluster:
        def comp_tab(meth, r):
            cc = comp[(comp.method == meth) & (comp.resolution == r)].drop(columns=["method", "resolution"])
            return table(cc.set_index("cluster"))

        cluster_sections = "".join(
            f'<h3>{mn}</h3><img src="data:image/png;base64,{b64(out / f"clusters_{m}.png")}" alt="clusters, {mn}">'
            + "".join(f'<details><summary>{mn}, r={r}: composition table and gene dot plot</summary>'
                      f'<div class=tw>{comp_tab(m, r)}</div>'
                      f'<img src="data:image/png;base64,{b64(out / f"dotplot_{m}_r{r}.png")}" alt="dot plot {mn} r{r}">'
                      f'</details>' for r in RESOLUTIONS)
            for m, mn in METHODS.items())
        cluster_html = f"""<h2 id="clusters">Clustering on each method</h2>
<p class=sub>Leiden (igraph) on each method's own neighbour graph (the graph its UMAP was built from). Bars: % of
all ctr (blue) or KO (red) cells that fall in each cluster; black dots = the two replicates. A cluster with a red
bar and no blue bar is KO-specific, and vice versa. Expand a resolution for the full table (KO_share_norm =
KO share after equalising genotype sizes; 50 = no enrichment) and a dot plot of the requested genes.</p>
<h3>Number of clusters</h3><div class=tw>{table(ncl_w)}</div>
<h3>Agreement between methods (adjusted Rand index; 1 = identical clusters)</h3><div class=tw>{table(ari)}</div>
{cluster_sections}"""
        cluster_link = '<a href="#clusters">Clustering</a> · '
    gene_sections = "".join(
        f'<h2 id="genes_{m}">Genes by genotype: {mn} UMAP</h2>'
        f'<img src="data:image/png;base64,{b64(out / f"genes_by_genotype_{m}.png")}" alt="genes, {mn}">'
        for m, mn in METHODS.items())
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{day} Summary</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1300px;margin:0 auto;padding:24px 16px 60px}}
h3{{font-size:16px;margin:22px 0 6px}} details{{margin:6px 0}} summary{{cursor:pointer;color:var(--mut)}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:30px 0 8px}}
.sub{{color:var(--mut);margin:0 0 16px}} .tw{{overflow-x:auto}}
table.t{{border-collapse:collapse;font-size:13px}}
table.t th,table.t td{{padding:3px 9px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px}}
</style></head><body><main>
<h1>{day}: summary and marker genes on each single-day UMAP</h1>
<p class=sub>{a.n_obs:,} cells, 4 samples · UMAPs from 5.6, this day's cells only: no integration (log PCA),
Harmony (sample_id), scVI (sample_id) · expression values are the same in every method; only the layout differs ·
expression = SoupX counts, normalize_total 1e4, log1p · out_all/day_analysis/{day}</p>
<h2>Samples</h2><div class=tw>{table(per_sample)}</div>
<p class=sub>{cluster_link}Jump to genes on: {" · ".join(f'<a href="#genes_{m}">{mn}</a>' for m, mn in METHODS.items())}</p>
<h2>UMAP by genotype (rows: no integration, Harmony, scVI)</h2><img src="data:image/png;base64,{b64(out / 'umap_genotype.png')}" alt="genotype UMAP">
<h2>Requested genes by genotype (table)</h2>
{"<p class=sub>Not found: " + ", ".join(missing) + "</p>" if missing else ""}
<p class=sub>pct_expr = % of cells with any count; mean_log = mean log1p-normalized expression over all cells
of that genotype. ISL1 added as a check of the KO.</p>
<div class=tw>{table(gtab, index=False)}</div>
{gene_sections}
<h2>Dot plot by sample</h2>
<img src="data:image/png;base64,{b64(out / 'dotplot_genes.png')}" alt="dot plot">
{cluster_html}
</main></body></html>"""
    (out / f"{day}_summary.html").write_text(html)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
