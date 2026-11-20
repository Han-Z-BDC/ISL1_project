"""Step 7.3c: RNA velocity of the D15 endocrine cells with scVelo; arrows on the 2D and 3D UMAPs (after 7.3b).

Decisions (user, 2026-11-19): velocity next to Slingshot on D15; arrows on both the 2D and the 3D UMAP.
  counts     spliced / unspliced from the 4 velocyto looms (7.3b), the 5,994 cells of the D15 endocrine subset
  graph      the subset's own SCT B PCs 1-30 (as 7.1d and Slingshot), 30 neighbours for the moments
  genes      scVelo filter_and_normalize: >= 20 shared counts, top 2,000 genes
  model      dynamical model (recover_dynamics), velocity, velocity graph, latent time, velocity confidence
  layouts    velocities projected onto the 7.1d 2D UMAP and 3D UMAP (same graph as the clustering)
Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
  scvelo.h5ad                 spliced / unspliced / velocity layers, graph, latent time (input of 7.3d)
  scvelo_cells.csv            barcode, label, genotype, latent time, velocity confidence, length
  scvelo_stream_2d.png        streamlines on the 2D UMAP: all cells | control | KO, coloured by subcluster name
  scvelo_grid_2d.png          averaged arrows on a grid, same panels
  scvelo_latent_time.png      latent time and velocity confidence on the 2D UMAP
  scvelo_3d.html              interactive 3D UMAP: cells by subcluster name + cone arrows averaged on a 3D grid;
                              all / control / KO (plotly.js from cdnjs)
Run: sbatch analysis/scripts/7.3c_scvelo.sb   (CELLRANK_ENV python: scvelo 0.3.4)
"""
import json
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
import scvelo as scv

OUT_DIR, REPO_DIR = Path(os.environ["OUT_DIR"]), Path(os.environ["REPO_DIR"])
DD = OUT_DIR / "day_analysis" / "D15"
SUB, TR = DD / "subcluster_endocrine", DD / "trajectory"
SAMPLES = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
SEED, NCPU = 1234, int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948', '#17becf',
       '#8c564b', '#7f7f7f', '#bcbd22']


def load():
    looms = []
    for s in SAMPLES:
        l = ad.read_loom(OUT_DIR / "velocity" / "D15" / s / f"{s}.loom")   # spliced as X; scvelo 0.3 has no scv.read
        # velocyto names cells "<sample>:<barcode without -1>x"
        l.obs_names = [f"Diff37_D15-{s}#" + n.split(":")[1][:-1] + "-1" for n in l.obs_names]
        l.var_names_make_unique()
        looms.append(l)
    lo = ad.concat(looms, join="inner")
    sc_ = pd.read_csv(SUB / "subclusters.csv", index_col=0)
    pca = pd.read_csv(SUB / "SCT_B_pca.csv.gz", index_col=0).reindex(sc_.index)
    nm = pd.read_csv(REPO_DIR / "study_config" / "annotation" / "D15_endocrine_r0.6.tsv", sep="\t", dtype=str,
                     keep_default_na=False)
    missing = sc_.index.difference(lo.obs_names)
    assert len(missing) == 0, f"{len(missing)} subset cells not in the looms, e.g. {list(missing[:3])}"
    a = lo[sc_.index].copy()
    a.obs["label"] = pd.Categorical(sc_.sub_r0p6.astype(str).map(dict(zip(nm.subcluster, nm.label))),
                                    categories=list(nm.label))
    a.obs["genotype"] = np.where(a.obs_names.str.contains("-KO_"), "KO", "ctr")
    a.obsm["X_pca"] = pca.iloc[:, :30].to_numpy()
    a.obsm["X_umap_2d"] = sc_[["UMAP1", "UMAP2"]].to_numpy()
    a.obsm["X_umap_3d"] = sc_[["UMAP3D_1", "UMAP3D_2", "UMAP3D_3"]].to_numpy()
    a.uns["label_colors"] = PAL[:len(nm)]
    return a


def grid_arrows(X, V, n=14, min_cells=8):
    """Average velocity of the cells in each grid box (any number of dimensions)."""
    lo, hi = X.min(0), X.max(0)
    idx = np.floor((X - lo) / (hi - lo + 1e-9) * n).astype(int)
    key = pd.Series([tuple(r) for r in idx])
    rows = []
    for k, ix in key.groupby(key).groups.items():
        ix = np.asarray(ix)
        if len(ix) >= min_cells:
            rows.append(np.r_[X[ix].mean(0), V[ix].mean(0)])
    return np.array(rows)


def main():
    TR.mkdir(parents=True, exist_ok=True)
    a = load()
    print(a)
    scv.pp.filter_and_normalize(a, min_shared_counts=20, n_top_genes=2000)
    sc.pp.neighbors(a, use_rep="X_pca", n_neighbors=30, random_state=SEED)
    scv.pp.moments(a, n_pcs=None, n_neighbors=None)
    scv.tl.recover_dynamics(a, n_jobs=NCPU)
    scv.tl.velocity(a, mode="dynamical")
    scv.tl.velocity_graph(a, n_jobs=NCPU)
    scv.tl.latent_time(a)
    scv.tl.velocity_confidence(a)
    for b in ["umap_2d", "umap_3d"]:
        scv.tl.velocity_embedding(a, basis=b)
    print("velocity genes:", int(a.var["velocity_genes"].sum()))
    a.write_h5ad(TR / "scvelo.h5ad")
    a.obs[["label", "genotype", "latent_time", "velocity_confidence", "velocity_length"]].to_csv(TR / "scvelo_cells.csv")

    ko = (a.obs.genotype == "KO").to_numpy()
    panels = [("all cells", np.ones(a.n_obs, bool)), ("control", ~ko), ("ISL1 KO", ko)]
    for kind in ["stream", "grid"]:
        fig, axs = plt.subplots(1, 3, figsize=(27, 8.5))
        for ax, (t, m) in zip(axs, panels):
            b = a[m].copy()
            f = scv.pl.velocity_embedding_stream if kind == "stream" else scv.pl.velocity_embedding_grid
            f(b, basis="umap_2d", color="label", ax=ax, show=False, title=f"D15 endocrine: {t} ({m.sum():,})",
              legend_loc="right margin" if t == "ISL1 KO" else "on data", legend_fontsize=8, size=25, alpha=0.6)
        fig.tight_layout(); fig.savefig(TR / f"scvelo_{kind}_2d.png", dpi=90); plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(17, 7.5))
    scv.pl.scatter(a, basis="umap_2d", color="latent_time", color_map="viridis", ax=axs[0], show=False, size=25,
                   title="latent time")
    scv.pl.scatter(a, basis="umap_2d", color="velocity_confidence", color_map="coolwarm", ax=axs[1], show=False,
                   size=25, title="velocity confidence")
    fig.tight_layout(); fig.savefig(TR / "scvelo_latent_time.png", dpi=90); plt.close(fig)

    # ---- 3D: cells + cone arrows averaged on a 3D grid, per panel
    X3, V3 = a.obsm["X_umap_3d"], a.obsm["velocity_umap_3d"]
    labs = list(a.obs.label.cat.categories)
    lab_i = a.obs.label.cat.codes.to_numpy()
    traces, groups = [], []
    for t, m in panels:
        g = []
        for i, l in enumerate(labs):
            s = m & (lab_i == i)
            if s.any():
                traces.append({"type": "scatter3d", "mode": "markers", "name": f"{l} ({s.sum()})", "legendgroup": l,
                               "x": X3[s, 0].round(2).tolist(), "y": X3[s, 1].round(2).tolist(),
                               "z": X3[s, 2].round(2).tolist(), "hoverinfo": "name",
                               "marker": {"size": 1.8, "color": PAL[i % len(PAL)], "opacity": 0.45}})
                g.append(len(traces) - 1)
        G = grid_arrows(X3[m], V3[m])
        if len(G):
            traces.append({"type": "cone", "x": G[:, 0].round(3).tolist(), "y": G[:, 1].round(3).tolist(),
                           "z": G[:, 2].round(3).tolist(), "u": G[:, 3].round(4).tolist(),
                           "v": G[:, 4].round(4).tolist(), "w": G[:, 5].round(4).tolist(), "sizemode": "scaled",
                           "sizeref": 1.6, "anchor": "tail", "showscale": False, "name": "velocity",
                           "colorscale": [[0, "#222"], [1, "#222"]], "hoverinfo": "skip"})
            g.append(len(traces) - 1)
        groups.append(g)
    for i, tr in enumerate(traces):
        tr["visible"] = i in groups[0]
    buttons = [{"label": t, "method": "update",
                "args": [{"visible": [i in groups[k] for i in range(len(traces))]}]} for k, (t, _) in enumerate(panels)]
    layout = {"height": 820, "legend": {"itemsizing": "constant"}, "margin": {"l": 0, "r": 0, "t": 50, "b": 0},
              "title": "D15 endocrine: RNA velocity on the 3D UMAP (cones = mean velocity of the cells in a 3D grid box)",
              "updatemenus": [{"buttons": buttons, "x": 0, "y": 1.06, "xanchor": "left", "type": "buttons",
                               "direction": "right"}],
              "scene": {"xaxis": {"title": "UMAP1"}, "yaxis": {"title": "UMAP2"}, "zaxis": {"title": "UMAP3"}}}
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>D15 Velocity 3D</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js"></script>
<style>body{{margin:0;font:14px system-ui,sans-serif;background:#fff}}</style></head>
<body><div id="p"></div><script>
Plotly.newPlot("p", {json.dumps(traces, separators=(",", ":"))}, {json.dumps(layout)}, {{responsive: true, displaylogo: false}});
</script></body></html>"""
    (TR / "scvelo_3d.html").write_text(html)
    print(a.obs.groupby(["label", "genotype"], observed=True)[["latent_time", "velocity_confidence"]].median().round(2).to_string())
    print("wrote", TR)


if __name__ == "__main__":
    main()
