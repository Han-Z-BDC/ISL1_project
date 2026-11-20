"""Step 7.4c: RNA velocity of the D15 + D22 endocrine integration with scVelo, on the scVI UMAPs (after 7.4 prep and
7.3b velocyto for D15 and D22).

Decisions (user, 2026-11-19): the 7.3 trajectory set on the 7.2 integration, scVI UMAP as the primary UMAP. Same
settings as the D15 run (7.3c) except where the integration changes them:
  counts     spliced / unspliced from the 8 velocyto looms (D15 and D22, 7.3b), the 21,063 cells of 7.2
  graph      scVI latent (30 dimensions), 30 neighbours for the moments (D15 used its own SCT B PCs; here the
             integrated space, so the moments of a cell can mix D15 and D22 neighbours)
  genes      scVelo filter_and_normalize: >= 20 shared counts, top 2,000 genes
  model      dynamical model, velocity, velocity graph, latent time, velocity confidence
  layouts    velocities projected onto the scVI 2D UMAP (7.2) and the scVI 3D UMAP (7.4 prep)
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/trajectory/ (reruns overwrite):
  scvelo.h5ad, scvelo_cells.csv (barcode, label, stage, genotype, latent time, velocity confidence, length)
  scvelo_stream_2d.png / scvelo_grid_2d.png   rows = all cells | D15 | D22; columns = all | control | KO
  scvelo_latent_time.png, scvelo_3d.html (cells by label + cone arrows on a 3D grid; panels as the 2D rows x cols)
Run: sbatch analysis/scripts/7.4c_scvelo.sb   (CELLRANK_ENV python: scvelo 0.3.4)
With INT_SUBSET=<name> (user, 2026-11-19: alpha_beta_EC), everything reads and writes the 7.2 subset folder
integration_D15_D22_endocrine_<name>/ instead: its cells, scVI latent and scVI UMAP.
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
SUBSET = os.environ.get("INT_SUBSET", "")   # 7.2 subset run (see 7.2a); empty = the full integration
INT = OUT_DIR / "day_analysis" / ("integration_D15_D22_endocrine" + (f"_{SUBSET}" if SUBSET else ""))
TR = INT / "trajectory"
SAMPLES = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
SEED, NCPU = 1234, int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
D15_PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948', '#17becf',
           '#8c564b', '#7f7f7f', '#bcbd22']


def palette(labels):
    """D22 labels keep their study colours (label_colours.tsv); D15 states get the 7.3c palette."""
    lc = pd.read_csv(REPO_DIR / "study_config" / "annotation" / "label_colours.tsv", sep="\t")
    known = dict(zip(lc.label, lc.colour))
    d15 = [l for l in labels if l.startswith("D15 ")]
    return [known.get(l) or D15_PAL[d15.index(l) % len(D15_PAL)] for l in labels]


def load():
    looms = []
    for day in ["D15", "D22"]:
        for s in SAMPLES:
            l = scv.read(OUT_DIR / "velocity" / day / s / f"{s}.loom", cache=False)
            # velocyto names cells "<sample>:<barcode without -1>x"
            l.obs_names = [f"Diff37_{day}-{s}#" + n.split(":")[1][:-1] + "-1" for n in l.obs_names]
            l.var_names_make_unique()
            looms.append(l)
    lo = ad.concat(looms, join="inner")
    ce = pd.read_csv(TR / "cells.csv", index_col=0)
    lat = pd.read_csv(TR / "scvi_latent.csv.gz", index_col=0).reindex(ce.index)
    missing = ce.index.difference(lo.obs_names)
    assert len(missing) == 0, f"{len(missing)} cells not in the looms, e.g. {list(missing[:3])}"
    a = lo[ce.index].copy()
    d15 = sorted(l for l in ce.label.unique() if l.startswith("D15 "))
    d22 = [l for l in ["SC_alpha", "SC_beta", "SC_delta", "GHRL+", "ENP", "SC_EC1", "SC_EC2"] if l in set(ce.label)]
    cats = d15 + d22 + sorted(set(ce.label) - set(d15) - set(d22))
    a.obs["label"] = pd.Categorical(ce.label.to_numpy(), categories=cats)
    a.obs["stage"] = ce.stage.to_numpy(); a.obs["genotype"] = ce.genotype.to_numpy()
    a.obsm["X_scvi"] = lat.to_numpy()
    a.obsm["X_umap_2d"] = ce[["UMAP1", "UMAP2"]].to_numpy()
    a.obsm["X_umap_3d"] = ce[["UMAP3D_1", "UMAP3D_2", "UMAP3D_3"]].to_numpy()
    a.uns["label_colors"] = palette(cats)
    return a


def grid_arrows(X, V, n=16, min_cells=8):
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
    sc.pp.neighbors(a, use_rep="X_scvi", n_neighbors=30, random_state=SEED)
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
    a.obs[["label", "stage", "genotype", "latent_time", "velocity_confidence", "velocity_length"]].to_csv(
        TR / "scvelo_cells.csv")

    ko = (a.obs.genotype == "KO").to_numpy()
    st = a.obs.stage.to_numpy()
    rows = [("all days", np.ones(a.n_obs, bool)), ("D15", st == "D15"), ("D22", st == "D22")]
    cols = [("all", np.ones(a.n_obs, bool)), ("control", ~ko), ("ISL1 KO", ko)]
    for kind in ["stream", "grid"]:
        fig, axs = plt.subplots(3, 3, figsize=(27, 25))
        for i, (rt, rm) in enumerate(rows):
            for j, (ct, cm) in enumerate(cols):
                m = rm & cm
                b = a[m].copy()
                f = scv.pl.velocity_embedding_stream if kind == "stream" else scv.pl.velocity_embedding_grid
                f(b, basis="umap_2d", color="label", ax=axs[i, j], show=False,
                  title=f"{rt}, {ct} ({m.sum():,})", legend_loc="right margin" if j == 2 else "none",
                  legend_fontsize=8, size=12, alpha=0.6)
        fig.tight_layout(); fig.savefig(TR / f"scvelo_{kind}_2d.png", dpi=80); plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(25, 7.5))
    scv.pl.scatter(a, basis="umap_2d", color="latent_time", color_map="viridis", ax=axs[0], show=False, size=12,
                   title="latent time")
    scv.pl.scatter(a, basis="umap_2d", color="velocity_confidence", color_map="coolwarm", ax=axs[1], show=False,
                   size=12, title="velocity confidence")
    scv.pl.scatter(a, basis="umap_2d", color="stage", ax=axs[2], show=False, size=12, title="day")
    fig.tight_layout(); fig.savefig(TR / "scvelo_latent_time.png", dpi=80); plt.close(fig)

    # ---- 3D: cells + cone arrows averaged on a 3D grid, one panel per day x genotype choice
    X3, V3 = a.obsm["X_umap_3d"], a.obsm["velocity_umap_3d"]
    labs = list(a.obs.label.cat.categories); pal = a.uns["label_colors"]
    lab_i = a.obs.label.cat.codes.to_numpy()
    traces, groups, names = [], [], []
    for rt, rm in rows:
        for ct, cm in cols:
            m = rm & cm
            g = []
            for i, l in enumerate(labs):
                s = m & (lab_i == i)
                if s.any():
                    traces.append({"type": "scatter3d", "mode": "markers", "name": f"{l} ({s.sum()})",
                                   "legendgroup": l, "x": X3[s, 0].round(2).tolist(), "y": X3[s, 1].round(2).tolist(),
                                   "z": X3[s, 2].round(2).tolist(), "hoverinfo": "name",
                                   "marker": {"size": 1.6, "color": pal[i], "opacity": 0.4}})
                    g.append(len(traces) - 1)
            G = grid_arrows(X3[m], V3[m])
            if len(G):
                traces.append({"type": "cone", "x": G[:, 0].round(3).tolist(), "y": G[:, 1].round(3).tolist(),
                               "z": G[:, 2].round(3).tolist(), "u": G[:, 3].round(4).tolist(),
                               "v": G[:, 4].round(4).tolist(), "w": G[:, 5].round(4).tolist(), "sizemode": "scaled",
                               "sizeref": 1.6, "anchor": "tail", "showscale": False, "name": "velocity",
                               "colorscale": [[0, "#222"], [1, "#222"]], "hoverinfo": "skip"})
                g.append(len(traces) - 1)
            groups.append(g); names.append(f"{rt}, {ct}")
    for i, tr in enumerate(traces):
        tr["visible"] = i in groups[0]
    buttons = [{"label": n, "method": "update", "args": [{"visible": [i in groups[k] for i in range(len(traces))]}]}
               for k, n in enumerate(names)]
    layout = {"height": 840, "legend": {"itemsizing": "constant"}, "margin": {"l": 0, "r": 0, "t": 60, "b": 0},
              "title": "D15 + D22 endocrine (scVI): RNA velocity on the 3D UMAP (cones = mean velocity per 3D grid box)",
              "updatemenus": [{"buttons": buttons, "x": 0, "y": 1.07, "xanchor": "left", "type": "dropdown"}],
              "scene": {"xaxis": {"title": "UMAP1"}, "yaxis": {"title": "UMAP2"}, "zaxis": {"title": "UMAP3"}}}
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>D15 D22 Velocity 3D</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js"></script>
<style>body{{margin:0;font:14px system-ui,sans-serif;background:#fff}}</style></head>
<body><div id="p"></div><script>
Plotly.newPlot("p", {json.dumps(traces, separators=(",", ":"))}, {json.dumps(layout)}, {{responsive: true, displaylogo: false}});
</script></body></html>"""
    (TR / "scvelo_3d.html").write_text(html)
    print(a.obs.groupby(["stage", "label", "genotype"], observed=True)[["latent_time", "velocity_confidence"]]
          .median().round(2).to_string())
    print("wrote", TR)


if __name__ == "__main__":
    main()
