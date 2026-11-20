"""Step 7.4d: CellRank 2 on the D15 + D22 endocrine integration: fates, terminal states, driver genes (after 7.4c
scVelo and 7.4a Slingshot).

Decisions (user, 2026-11-19): the 7.3 set on the 7.2 integration (scVI). Same design as the D15 run (7.3d):
  kernels    velocity   0.8 VelocityKernel (7.4c) + 0.2 ConnectivityKernel
             pseudotime 0.8 PseudotimeKernel on Slingshot pseudotime (7.4a, scVI latent; one value per cell =
                        lineage-weight-averaged pseudotime, scaled to 0-1) + 0.2 ConnectivityKernel
             (same 30-neighbour graph on the scVI latent for both)
  states     data-driven: GPCCA macrostates (Schur vectors, n 4-12) and predicted terminal states
             fixed:       terminal states = the clusters that end the 7.4a Slingshot lineages (D22 SC_beta, SC_alpha,
                          SC_delta, SC_EC1, SC_EC2 by design, plus any free ends Slingshot chose)
  fates      fate probabilities toward the fixed terminal states; summary = mean per label x genotype
  drivers    Pearson r between log1p expression (SoupX counts normalised to 10,000 per cell; genes detected in
             >= 10 cells) and the fate probability toward SC_alpha and toward SC_beta, in the early D15 cells
             (ENP1, ENP2, ENP_FEV_high), control cells only and all cells; p from Fisher z, FDR (BH)
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/trajectory/ (reruns overwrite):
  cellrank_{kernel}_fates.csv, _fate_summary.csv, _macrostates.csv, _drivers.csv, _fates.png (scVI 2D UMAP)
Run: sbatch analysis/scripts/7.4d_cellrank.sb
With INT_SUBSET=<name> (user, 2026-11-19: alpha_beta_EC), everything reads and writes the 7.2 subset folder
integration_D15_D22_endocrine_<name>/ instead: its cells, scVI latent and scVI UMAP.
"""
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import cellrank as cr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

OUT_DIR = Path(os.environ["OUT_DIR"])
SUBSET = os.environ.get("INT_SUBSET", "")   # 7.2 subset run (see 7.2a); empty = the full integration
INT = OUT_DIR / "day_analysis" / ("integration_D15_D22_endocrine" + (f"_{SUBSET}" if SUBSET else ""))
TR = INT / "trajectory"
EARLY = ["D15 ENP1", "D15 ENP2", "D15 ENP_FEV_high"]
DRIVE = ["SC_alpha", "SC_beta"]
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", ["#e8e8e8", "#fdbb84", "#e34a33", "#7f0000"])


def slingshot_ends():
    ln = pd.read_csv(TR / "slingshot_scvi_lineages.csv")
    return list(dict.fromkeys(c.split(" -> ")[-1] for c in ln.clusters))


def slingshot_time(a):
    s = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0).reindex(a.obs_names)
    pt = np.nan_to_num(s.filter(like="pt_").to_numpy()); w = np.nan_to_num(s.filter(like="w_").to_numpy())
    t = (pt * w).sum(1) / np.maximum(w.sum(1), 1e-12)
    return t / t.max()


def drivers(expr, genes, fates, mask, tag):
    rows = []
    X = expr[mask]
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    Xc = X - X.mean(0); sx = np.sqrt((Xc ** 2).sum(0))
    n = mask.sum()
    for l in DRIVE:
        if l not in fates:
            continue
        f = fates[l].to_numpy()[mask]; fc = f - f.mean()
        r = (Xc * fc[:, None]).sum(0) / (sx * np.sqrt((fc ** 2).sum()) + 1e-12)
        z = np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(n - 3)
        p = 2 * stats.norm.sf(np.abs(z))
        rows.append(pd.DataFrame({"gene": genes, "lineage": l, "cells": tag, "n": n, "r": r.round(4), "p": p,
                                  "FDR": stats.false_discovery_control(p)}))
    return pd.concat(rows)


def run(name, a, kernel, full, ends):
    g = cr.estimators.GPCCA(kernel)
    g.compute_schur(n_components=20)
    try:
        g.compute_macrostates(n_states=[4, 12], cluster_key="label")
        g.predict_terminal_states()
        ms = a.obs.assign(m=g.macrostates).dropna(subset=["m"])
        mt = ms.groupby("m", observed=True).agg(cells=("label", "size"),
                                                 main_label=("label", lambda x: x.value_counts().index[0]))
        mt["terminal"] = mt.index.isin(g.terminal_states.cat.categories)
        mt.to_csv(TR / f"cellrank_{name}_macrostates.csv")
        print(name, "data-driven macrostates:\n", mt.to_string())
    except Exception as e:   # keep going to the fixed terminal states
        print(name, "macrostates failed:", e)
    lab = a.obs.label.astype(str)
    g.set_terminal_states(states=lab.where(lab.isin(ends)).astype("category"), cluster_key="label")
    g.compute_fate_probabilities(n_jobs=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    fp = g.fate_probabilities
    fates = pd.DataFrame(np.asarray(fp), index=a.obs_names, columns=list(fp.names))
    out = pd.concat([a.obs[["label", "stage", "genotype"]], fates.round(4)], axis=1)
    out.to_csv(TR / f"cellrank_{name}_fates.csv")
    summ = out.groupby(["label", "genotype"], observed=True)[list(fates.columns)].mean().round(3)
    summ.insert(0, "cells", out.groupby(["label", "genotype"], observed=True).size())
    summ.to_csv(TR / f"cellrank_{name}_fate_summary.csv")
    print(name, "mean fate probability by label x genotype:\n", summ.to_string())
    ko = (a.obs.genotype == "KO").to_numpy()
    early = lab.isin(EARLY).to_numpy()
    dr = pd.concat([drivers(full.X, full.var_names, fates, early & ~ko, "control early"),
                    drivers(full.X, full.var_names, fates, early, "all early")])
    dr = dr.sort_values(["lineage", "cells", "r"], ascending=[True, True, False], kind="stable")
    dr.to_csv(TR / f"cellrank_{name}_drivers.csv", index=False)
    for l in DRIVE:
        top = dr[(dr.lineage == l) & (dr.cells == "control early")].head(20).gene.tolist()
        print(name, l, "top drivers (control early D15 cells):", ", ".join(top))
    u = a.obsm["X_umap_2d"]
    fig, axs = plt.subplots(len(fates.columns), 2, figsize=(11, 4.6 * len(fates.columns)), squeeze=False)
    for i, l in enumerate(fates.columns):
        for j, (t, m) in enumerate([("control", ~ko), ("ISL1 KO", ko)]):
            ax = axs[i, j]
            ax.scatter(u[~m, 0], u[~m, 1], c="#f0f0f0", s=2, linewidths=0, rasterized=True)
            v = fates[l].to_numpy(); o = np.where(m)[0]; o = o[np.argsort(v[o], kind="stable")]
            s_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=3, cmap=FEAT_CMAP, vmin=0, vmax=1, linewidths=0,
                            rasterized=True)
            ax.set_title(f"{name}: fate -> {l}, {t}", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
            if j == 1:
                plt.colorbar(s_, ax=ax, fraction=0.04)
    fig.tight_layout(); fig.savefig(TR / f"cellrank_{name}_fates.png", dpi=80); plt.close(fig)


def main():
    a = ad.read_h5ad(TR / "scvelo.h5ad")
    a.obs["label"] = a.obs["label"].astype("category")
    a.obs["slingshot_time"] = slingshot_time(a)
    ends = slingshot_ends()
    print("terminal states (Slingshot lineage ends):", ends)
    full = ad.read_h5ad(INT / "merged_counts.h5ad")[a.obs_names].copy()
    keep = np.asarray((full.X > 0).sum(0)).ravel() >= 10
    sc.pp.normalize_total(full, target_sum=1e4); sc.pp.log1p(full)     # totals over all genes, then filter
    full = full[:, keep].copy()
    print(a, "\ndriver genes tested:", full.n_vars)
    ck = cr.kernels.ConnectivityKernel(a).compute_transition_matrix()
    vk = cr.kernels.VelocityKernel(a).compute_transition_matrix()
    pk = cr.kernels.PseudotimeKernel(a, time_key="slingshot_time").compute_transition_matrix()
    run("velocity", a, 0.8 * vk + 0.2 * ck, full, ends)
    run("pseudotime", a, 0.8 * pk + 0.2 * ck, full, ends)
    print("wrote", TR)


if __name__ == "__main__":
    main()
