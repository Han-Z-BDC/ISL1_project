"""Step 7.1d: subcluster one compartment of one day (cells taken from 7.1b clusters), with 2D and 3D UMAPs.

Request (user, 2026-11-18): D15, resolution 0.2 kept for the overview; the endocrine compartment (7.1b clusters
0, 1, 5 at r = 0.2) taken out and subclustered, with a 3D UMAP as well.

Two stages, run in order by 7.1d_subcluster.sb:
  prep     write the subset's barcodes (cells.txt) from 7.1b clusters.csv
  (R)      7.1c_day_sct.R in subset mode: SCT B (SCTransform v2, glmGamPoi, S - G2M regressed, 3,000 variable
           features without MT- genes) refit on the subset's cells only, so the variable genes and PCs are the
           compartment's own -> SCT_B_pca.csv.gz in the output folder
  cluster  30 PCs, 15 neighbours (as 7.1b); UMAP in 2D and 3D (same graph, seed 1234); Leiden (igraph,
           2 iterations, seed 1234) at RESOLUTIONS; composition, markers (Wilcoxon, subcluster vs rest of the
           subset, unregressed log1p), dot plots, interactive 3D UMAP (plotly) and a report.

Writes ${OUT_DIR}/day_analysis/{day}/subcluster_{name}/ (reruns overwrite):
  cells.txt, SCT_B_pca.csv.gz, SCT_B_variable_features.csv, SCT_B_pca_variance.csv   (prep + R)
  subclusters.csv            barcode -> parent cluster, subcluster at every resolution, UMAP 2D + 3D coordinates
  n_clusters.csv, lineage.csv, composition_r*.csv, markers_r*.csv
  umap_parent.png            subset UMAP coloured by the 7.1b parent cluster, sample and genotype
  umap_subclusters.png       all resolutions;  umap_r*_large.png  one resolution + genotype split
  umap_genes.png             endocrine genes on the 2D UMAP
  composition_r*.png, dotplot_genes_r*.png, dotplot_markers_r*.png
  umap3d.html                interactive 3D UMAP (colour by parent / genotype / sample / each resolution / genes)
  {day}_{name}_subclusters.html
Run: sbatch analysis/scripts/7.1d_subcluster.sb D15 endocrine 0.2 0,1,5
"""
import base64
import importlib.util
import json
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
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
OUT_DIR = Path(os.environ["OUT_DIR"])
SEED = 1234
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
N_PCS, N_NEIGHBORS = 30, 15
GENES = ["CHGA", "NEUROG3", "NKX2-2", "NEUROD1", "INSM1", "PAX4", "ARX", "IRX2", "PAX6", "ISL1", "MAFB", "MAFA",
         "INS", "GCG", "SST", "GHRL", "PPY", "PDX1", "NKX6-1", "FEV", "LMX1A", "TPH1", "SLC18A1", "DDC",
         "PHOX2A", "SCGN", "FOXP2", "TOX", "DCC", "MKI67"]
SAMPLE_ORDER = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
SAMPLE_PAL = {"ctr_rep1": "#2f5d9e", "ctr_rep2": "#7fa7dd", "KO_rep1": "#b23b2a", "KO_rep2": "#eb9a84"}
GENO_PAL = {"ctr": "#3b6fb6", "KO": "#d0573f"}
GREY = "#e3e3e3"

# drawing helpers shared with 7.1b (palette, draw_clusters, b64)
_spec = importlib.util.spec_from_file_location("s71b", HERE / "7.1b_day_cluster.py")
s71b = importlib.util.module_from_spec(_spec)
_argv = sys.argv; sys.argv = [_argv[0], "_"]; _spec.loader.exec_module(s71b); sys.argv = _argv
palette, draw_clusters, b64 = s71b.palette, s71b.draw_clusters, s71b.b64


def key(r):
    return f"sub_r{str(r).replace('.', 'p')}"


def paths(day, name):
    dd = OUT_DIR / "day_analysis" / day
    return dd, dd / f"subcluster_{name}"


def prep(day, name, res, clusters):
    dd, out = paths(day, name)
    out.mkdir(parents=True, exist_ok=True)
    for p in out.iterdir():
        if p.is_file():
            p.unlink()
    col = f"leiden_r{str(res).replace('.', 'p')}"
    cl = pd.read_csv(dd / "clustering" / "clusters.csv", index_col=0)[col]
    cells = cl.index[cl.isin(clusters)]
    (out / "cells.txt").write_text("\n".join(cells) + "\n")
    (out / "parent.json").write_text(json.dumps({"resolution": res, "clusters": clusters, "column": col}))
    print(f"{day} {name}: {len(cells)} cells from {col} clusters {clusters}")
    print(cl[cl.isin(clusters)].value_counts().sort_index().to_string())


def umap3d_html(a, out, title):
    """Self-contained plotly page; colour switch by dropdown. Loads plotly.js from cdnjs."""
    u = a.obsm["X_umap_3d"]
    geno = a.obs["genotype"].astype(str)
    layers = {"7.1b parent cluster": a.obs["parent"].astype(str), "genotype": geno,
              "sample": a.obs["sample"].astype(str)}
    for r in RESOLUTIONS:
        layers[f"subclusters r{r}"] = a.obs[key(r)].astype(str)
    cat_traces, buttons, vis_all = [], [], []
    for lab, s in layers.items():
        cats = sorted(s.unique(), key=lambda v: (not v.isdigit(), int(v) if v.isdigit() else v))
        pal = GENO_PAL if lab == "genotype" else SAMPLE_PAL if lab == "sample" else palette(cats)
        tr = []
        for c in cats:
            m = (s == c).to_numpy()
            tr.append({"type": "scatter3d", "mode": "markers", "name": f"{c} ({m.sum()})",
                       "x": u[m, 0].round(3).tolist(), "y": u[m, 1].round(3).tolist(), "z": u[m, 2].round(3).tolist(),
                       "marker": {"size": 1.6, "color": pal[c], "opacity": 0.85}, "hoverinfo": "name",
                       "visible": False})
        vis_all.append(len(tr)); cat_traces.extend(tr)
    # genes: one continuous trace each (log1p)
    genes = [g for g in GENES if g in a.var_names]
    X = a[:, genes].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    for i, g in enumerate(genes):
        cat_traces.append({"type": "scatter3d", "mode": "markers", "name": g,
                           "x": u[:, 0].round(3).tolist(), "y": u[:, 1].round(3).tolist(), "z": u[:, 2].round(3).tolist(),
                           "marker": {"size": 1.6, "color": X[:, i].round(2).tolist(), "colorscale": "Viridis",
                                      "colorbar": {"title": g, "thickness": 12}, "opacity": 0.85},
                           "hoverinfo": "skip", "visible": False})
        vis_all.append(1)
    names = list(layers) + [f"gene: {g}" for g in genes]
    start = 0
    for lab, n in zip(names, vis_all):
        vis = [False] * len(cat_traces)
        vis[start:start + n] = [True] * n
        buttons.append({"label": lab, "method": "update", "args": [{"visible": vis}]})
        start += n
    first = names.index("subclusters r0.6") if "subclusters r0.6" in names else 0
    s0 = sum(vis_all[:first])
    for t in cat_traces[s0:s0 + vis_all[first]]:
        t["visible"] = True
    layout = {"title": title, "height": 780, "legend": {"itemsizing": "constant"},
              "scene": {"xaxis": {"title": "UMAP1"}, "yaxis": {"title": "UMAP2"}, "zaxis": {"title": "UMAP3"}},
              "updatemenus": [{"buttons": buttons, "active": first, "x": 0, "y": 1.08, "xanchor": "left"}],
              "margin": {"l": 0, "r": 0, "t": 60, "b": 0}}
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js"></script>
<style>body{{margin:0;font:14px system-ui,sans-serif;background:#fff}} #p{{width:100%}}</style></head>
<body><div id="p"></div><script>
Plotly.newPlot("p", {json.dumps(cat_traces)}, {json.dumps(layout)}, {{responsive: true}});
</script></body></html>"""
    (out / "umap3d.html").write_text(html)


def cluster(day, name):
    dd, out = paths(day, name)
    par = json.loads((out / "parent.json").read_text())
    cells = (out / "cells.txt").read_text().split()
    a = ad.read_h5ad(dd / f"{day}_merged.h5ad")
    a = a[cells].copy()
    a.obs = a.obs[[c for c in a.obs.columns if not c.startswith("leiden_r")]]
    a.obs["parent"] = pd.read_csv(dd / "clustering" / "clusters.csv", index_col=0)[par["column"]] \
        .reindex(a.obs_names).astype(str).to_numpy()
    pca = pd.read_csv(out / "SCT_B_pca.csv.gz", index_col=0)
    assert (pca.index == a.obs_names).all()
    a.obsm["X_pca_sub"] = pca.to_numpy()
    sc.pp.neighbors(a, use_rep="X_pca_sub", n_pcs=N_PCS, n_neighbors=N_NEIGHBORS, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    a.obsm["X_umap_2d"] = a.obsm["X_umap"].copy()
    sc.tl.umap(a, random_state=SEED, n_components=3)
    a.obsm["X_umap_3d"] = a.obsm["X_umap"].copy()
    a.obsm["X_umap"] = a.obsm["X_umap_2d"]
    for r in RESOLUTIONS:
        sc.tl.leiden(a, resolution=r, key_added=key(r), random_state=SEED, flavor="igraph",
                     n_iterations=2, directed=False)
    n_cl = {r: a.obs[key(r)].nunique() for r in RESOLUTIONS}
    print(f"{day} {name}: {a.n_obs} cells; subclusters {n_cl}")
    u, u3 = a.obsm["X_umap_2d"], a.obsm["X_umap_3d"]
    tab = a.obs[["parent"] + [key(r) for r in RESOLUTIONS]].copy()
    tab[["UMAP1", "UMAP2"]] = u; tab[["UMAP3D_1", "UMAP3D_2", "UMAP3D_3"]] = u3
    tab.to_csv(out / "subclusters.csv")
    pd.Series(n_cl, name="n_clusters").rename_axis("resolution").to_csv(out / "n_clusters.csv")
    rows = []
    for lo, hi in zip(RESOLUTIONS[:-1], RESOLUTIONS[1:]):
        ct = pd.crosstab(a.obs[key(hi)], a.obs[key(lo)])
        for c in ct.index:
            p = ct.loc[c].idxmax()
            rows.append({"resolution": hi, "cluster": int(c), "parent_resolution": lo, "parent": int(p),
                         "pct_cells_from_parent": round(100 * ct.loc[c, p] / ct.loc[c].sum(), 1)})
    pd.DataFrame(rows).to_csv(out / "lineage.csv", index=False)

    geno = a.obs["genotype"].astype(str).to_numpy()
    n_geno = pd.Series(geno).value_counts()
    # parent / sample / genotype on the subset UMAP
    fig, axs = plt.subplots(1, 3, figsize=(24, 7.5))
    draw_clusters(axs[0], u, a.obs["parent"].astype(str), f"7.1b parent cluster (r = {par['resolution']})")
    for ax, (lab, s, pal) in zip(axs[1:], [("sample", a.obs["sample"].astype(str), SAMPLE_PAL),
                                          ("genotype", a.obs["genotype"].astype(str), GENO_PAL)]):
        idx = np.random.default_rng(SEED).permutation(len(s))
        ax.scatter(u[idx, 0], u[idx, 1], c=s.map(pal).to_numpy()[idx], s=1.5, linewidths=0, rasterized=True)
        for k_, c_ in pal.items():
            ax.scatter([], [], c=c_, s=20, label=k_)
        ax.legend(frameon=False, markerscale=1.5); ax.set_title(lab); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(out / "umap_parent.png", dpi=110); plt.close(fig)

    fig, axs = plt.subplots(2, 3, figsize=(20, 12))
    for ax, r in zip(axs.flat, RESOLUTIONS):
        draw_clusters(ax, u, a.obs[key(r)].astype(str), f"resolution {r}: {n_cl[r]} subclusters")
    fig.suptitle(f"{day} {name}: subclusters (SCT B refit on the subset)")
    fig.tight_layout(); fig.savefig(out / "umap_subclusters.png", dpi=110); plt.close(fig)

    genes = [g for g in GENES if g in a.var_names]
    sc.pl.umap(a, color=genes, ncols=6, cmap="viridis", s=6, frameon=False, show=False, vmax="p99")
    plt.gcf().savefig(out / "umap_genes.png", dpi=90); plt.close("all")

    for r in RESOLUTIONS:
        k = key(r); cl = a.obs[k].astype(str)
        fig, axs = plt.subplots(1, 3, figsize=(25, 8))
        draw_clusters(axs[0], u, cl, f"resolution {r}: {n_cl[r]} subclusters")
        for ax, gt, lab in zip(axs[1:], ("ctr", "KO"), ("Control", "ISL1 KO")):
            draw_clusters(ax, u, cl, f"{lab} cells (other genotype grey)", mask=(geno == gt))
        fig.tight_layout(); fig.savefig(out / f"umap_r{r}_large.png", dpi=110); plt.close(fig)

        ct = pd.crosstab(a.obs[k], a.obs["sample"]).reindex(columns=SAMPLE_ORDER, fill_value=0)
        ct.index = ct.index.astype(int); ct = ct.sort_index()
        comp = ct.copy()
        comp["cells"] = ct.sum(axis=1)
        comp["pct_of_ctr_subset"] = (100 * ct[["ctr_rep1", "ctr_rep2"]].sum(axis=1) / n_geno["ctr"]).round(2)
        comp["pct_of_KO_subset"] = (100 * ct[["KO_rep1", "KO_rep2"]].sum(axis=1) / n_geno["KO"]).round(2)
        comp["KO_over_ctr"] = (comp["pct_of_KO_subset"] / comp["pct_of_ctr_subset"]).round(2)
        pc = pd.crosstab(a.obs[k], a.obs["parent"]); pc.index = pc.index.astype(int)
        comp = comp.join((100 * pc.div(pc.sum(axis=1), axis=0)).round(1).add_prefix("pct_from_parent_"))
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
        comp.to_csv(out / f"composition_r{r}.csv")
        pct = 100 * ct / ct.sum(axis=0)
        fig, ax = plt.subplots(figsize=(max(10, 0.6 * len(pct)), 4))
        x = np.arange(len(pct)); w = 0.2
        for i, s in enumerate(SAMPLE_ORDER):
            ax.bar(x + (i - 1.5) * w, pct[s], width=w, color=SAMPLE_PAL[s], label=s)
        ax.set_xticks(x); ax.set_xticklabels(pct.index); ax.set_ylabel(f"% of sample's {name} cells")
        ax.set_xlabel("subcluster"); ax.legend(fontsize=8, frameon=False, ncol=4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"{day} {name} resolution {r}: % of each sample's {name} cells per subcluster", fontsize=10)
        fig.tight_layout(); fig.savefig(out / f"composition_r{r}.png", dpi=110); plt.close(fig)

        a.obs[k] = pd.Categorical(a.obs[k], sorted(a.obs[k].unique(), key=int))
        sc.tl.rank_genes_groups(a, k, method="wilcoxon", n_genes=50, use_raw=False, pts=True)
        mk = sc.get.rank_genes_groups_df(a, None).groupby("group", observed=True).head(50)
        mk.to_csv(out / f"markers_r{r}.csv", index=False)
        dp = sc.pl.dotplot(a, genes, groupby=k, show=False, return_fig=True)
        dp.savefig(out / f"dotplot_genes_r{r}.png", dpi=110); plt.close("all")
        top = []
        for _, g in mk[mk["logfoldchanges"] > 0].groupby("group", observed=True, sort=False):
            for nm in g["names"].head(3):
                if nm not in top:
                    top.append(nm)
        dp = sc.pl.dotplot(a, top, groupby=k, show=False, return_fig=True, standard_scale="var")
        dp.savefig(out / f"dotplot_markers_r{r}.png", dpi=110); plt.close("all")
        a.uns.pop("rank_genes_groups", None)

    umap3d_html(a, out, f"{day} {name}: 3D UMAP")
    a.write_h5ad(out / f"{day}_{name}.h5ad", compression="gzip")

    secs = ""
    for r in RESOLUTIONS:
        comp = pd.read_csv(out / f"composition_r{r}.csv", index_col=0)
        mk = pd.read_csv(out / f"markers_r{r}.csv")
        top10 = (mk[mk["logfoldchanges"] > 0].groupby("group").head(10).groupby("group")["names"]
                 .apply(lambda v: ", ".join(v)).rename("top 10 markers (up)"))
        top10.index = top10.index.astype(int)
        secs += f"""<h2 id="r{r}">Resolution {r}: {n_cl[r]} subclusters</h2>
<img src="data:image/png;base64,{b64(out / f'umap_r{r}_large.png')}" alt="subclusters at resolution {r}, with genotype split">
<img src="data:image/png;base64,{b64(out / f'composition_r{r}.png')}" alt="composition at resolution {r}">
<div class=tw>{comp.to_html(border=0, classes="t")}</div><div class=tw>{top10.to_frame().to_html(border=0, classes="t")}</div>
<img src="data:image/png;base64,{b64(out / f'dotplot_genes_r{r}.png')}" alt="gene dotplot at resolution {r}">
<img src="data:image/png;base64,{b64(out / f'dotplot_markers_r{r}.png')}" alt="top marker dotplot at resolution {r}">"""
    nav = " · ".join(f'<a href="#r{r}">{r}</a>' for r in RESOLUTIONS)
    pc = ", ".join(str(c) for c in par["clusters"])
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{day} {name.capitalize()} Subclusters</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1240px;margin:0 auto;padding:24px 16px 60px}} a{{color:inherit}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:34px 0 8px}}
.sub{{color:var(--mut)}} .tw{{overflow-x:auto;margin:8px 0}}
table.t{{border-collapse:collapse;font-size:12.5px}} .t th,.t td{{padding:3px 7px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
.t td:last-child{{text-align:left}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>{day}: {name} subclusters</h1>
<p class=sub>{a.n_obs:,} cells = 7.1b clusters {pc} at resolution {par['resolution']}. SCT B refit on these cells only
(SCTransform v2, S&minus;G2M regressed, 3,000 variable features without MT- genes); 30 PCs, 15 neighbours; UMAP 2D and
3D from the same graph; Leiden (igraph, 2 iterations, seed {SEED}). Markers: Wilcoxon, subcluster vs the rest of the
subset, unregressed log1p. % columns are of the subset's ctr or KO cells. Interactive 3D UMAP: <code>umap3d.html</code> in
this folder. Resolutions: {nav}</p>
<img src="data:image/png;base64,{b64(out / 'umap_parent.png')}" alt="subset UMAP by parent cluster, sample, genotype">
<img src="data:image/png;base64,{b64(out / 'umap_subclusters.png')}" alt="subclusters at all resolutions">
<img src="data:image/png;base64,{b64(out / 'umap_genes.png')}" alt="endocrine genes on the subset UMAP">
{secs}
</main></body></html>"""
    (out / f"{day}_{name}_subclusters.html").write_text(html)
    print("wrote", out)


if __name__ == "__main__":
    stage, day, name = sys.argv[1], sys.argv[2], sys.argv[3]
    if stage == "prep":
        prep(day, name, float(sys.argv[4]), [int(c) for c in sys.argv[5].split(",")])
    elif stage == "cluster":
        cluster(day, name)
    else:
        sys.exit(f"unknown stage {stage}")
