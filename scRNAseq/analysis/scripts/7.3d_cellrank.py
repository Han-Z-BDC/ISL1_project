"""Step 7.3d: CellRank 2 on the D15 endocrine cells: fate probabilities, terminal states, driver genes (after 7.3c).

Decisions (user, 2026-11-19): CellRank 2 in its own env (CELLRANK_ENV); both kernels, velocity and Slingshot
pseudotime, so the result does not rest on one method.
Decision (user, 2026-11-20): besides the run on all cells, a control-only and a KO-only run. In the shared run the
graph links KO cells to control neighbours, so KO fates can borrow control paths (e.g. pre_beta1 -> pre_beta2); the
per-genotype runs use a graph built within the genotype, and their data-driven terminal states show where each
genotype's flow ends. The Slingshot trajectory stays shared (one pseudotime for both genotypes).
  cells      all       the 5,994 cells (shared graph)
             ctr / KO  one genotype; 30-neighbour graph rebuilt on that genotype's cells (same SCT B PCs 1-30);
                       velocities from the shared scVelo model (7.3c); pseudotime = the shared Slingshot pseudotime
  kernels    velocity   0.8 VelocityKernel (scVelo, 7.3c) + 0.2 ConnectivityKernel
             pseudotime 0.8 PseudotimeKernel on Slingshot pseudotime (7.3a, PCs 1-30, fixed ends; one value per cell
                        = lineage-weight-averaged pseudotime) + 0.2 ConnectivityKernel
  states     data-driven: GPCCA macrostates (Schur vectors, n chosen by eigengap, 4-12) and predicted terminal states;
                          when the chosen n would split a complex-conjugate eigenvalue pair, the search is retried with
                          the upper bound lowered by one (12 -> 11 -> ...); the bound used is in the macrostates table
             fixed:       all  the Slingshot ends pre_beta2, pre_alpha, pre_delta, pre_EC, ENP_PHOX2A+, plus KO_only
                          ctr  pre_beta2, pre_alpha, pre_delta, pre_EC, ENP_PHOX2A+ (KO_only has 11 control cells)
                          KO   pre_EC, ENP_PHOX2A+, KO_only, pre_alpha_KO (no KO cells in pre_alpha, pre_beta2,
                               pre_delta)
  fates      fate probabilities toward the fixed terminal states for every cell; summary = mean per label x genotype
  drivers    Pearson r between each gene's log1p expression (all genes detected in >= 10 cells; 7.1d X) and a fate
             probability, in the early cells (ENP1, ENP2, ENP_FEV_high); p from Fisher z, FDR (BH)
             all       fates pre_alpha and pre_beta2; control early cells and all early cells
             ctr / KO  every fixed terminal state; that genotype's early cells
Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite); {run} = {kernel} (all cells) or {kernel}_ctr /
{kernel}_KO:
  cellrank_{run}_fates.csv          barcode, label, genotype, data-driven macrostate, fate probability per terminal state
  cellrank_{run}_fate_summary.csv   mean fate probability per label x genotype
  cellrank_{run}_macrostates.csv    data-driven macrostates and terminal states (cells per state, by genotype, main
                                    label, search bound)
  cellrank_{run}_drivers.csv        gene, lineage, cells, r, p, FDR
  cellrank_{run}_fates.png          fate probabilities on the 2D UMAP (all: control | KO)
Run: sbatch analysis/scripts/7.3d_cellrank.sb
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

# CellRank 2.0.7 (VelocityKernel) calls np.testing.assert_array_equal(x=..., y=...); numpy >= 2 names these
# actual / desired, so pass them on by position (the env has numpy 2.4.6; left unchanged).
_assert_array_equal = np.testing.assert_array_equal


def _assert_array_equal_xy(*args, x=None, y=None, **kw):
    if x is not None or y is not None:
        args = (x, y) + args
    return _assert_array_equal(*args, **kw)


np.testing.assert_array_equal = _assert_array_equal_xy

OUT_DIR = Path(os.environ["OUT_DIR"])
DD = OUT_DIR / "day_analysis" / "D15"
SUB, TR = DD / "subcluster_endocrine", DD / "trajectory"
SEED = 1234
NCPU = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
ENDS = {"all": ["pre_beta2", "pre_alpha", "pre_delta", "pre_EC", "ENP_PHOX2A+", "KO_only"],
        "ctr": ["pre_beta2", "pre_alpha", "pre_delta", "pre_EC", "ENP_PHOX2A+"],
        "KO": ["pre_EC", "ENP_PHOX2A+", "KO_only", "pre_alpha_KO"]}
EARLY = ["ENP1", "ENP2", "ENP_FEV_high"]
DRIVE = ["pre_alpha", "pre_beta2"]
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", ["#e8e8e8", "#fdbb84", "#e34a33", "#7f0000"])


def slingshot_time(a):
    s = pd.read_csv(TR / "slingshot_pc30_cells.csv", index_col=0).reindex(a.obs_names)
    pt = s.filter(like="pt_").to_numpy(); w = s.filter(like="w_").to_numpy()
    pt = np.nan_to_num(pt); w = np.nan_to_num(w)
    t = (pt * w).sum(1) / np.maximum(w.sum(1), 1e-12)
    return t / t.max()


def drivers(expr, genes, fates, mask, tag, lineages):
    rows = []
    X = expr[mask]
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    Xc = X - X.mean(0); sx = np.sqrt((Xc ** 2).sum(0))
    n = mask.sum()
    for l in lineages:
        f = fates[l].to_numpy()[mask]; fc = f - f.mean()
        r = (Xc * fc[:, None]).sum(0) / (sx * np.sqrt((fc ** 2).sum()) + 1e-12)
        z = np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(n - 3)
        p = 2 * stats.norm.sf(np.abs(z))
        q = stats.false_discovery_control(p)
        rows.append(pd.DataFrame({"gene": genes, "lineage": l, "cells": tag, "n": n, "r": r.round(4), "p": p, "FDR": q}))
    return pd.concat(rows)


def macrostates(g):
    """Eigengap search over 4..hi states; lower hi when the chosen n would split a complex-conjugate pair."""
    for hi in range(12, 4, -1):
        try:
            g.compute_macrostates(n_states=[4, hi], cluster_key="label")
            return hi
        except ValueError as e:
            if "complex conjugate" not in str(e):
                raise
            print("  macrostates: upper bound", hi, "splits a complex-conjugate pair; retrying with", hi - 1)
    raise ValueError("no macrostate number between 4 and 5 avoids splitting complex-conjugate eigenvalues")


def run(name, a, kernel, full, cells):
    g = cr.estimators.GPCCA(kernel)
    g.compute_schur(n_components=20)
    a.obs["macrostate"] = pd.Series(np.nan, index=a.obs_names, dtype=object)
    try:
        hi = macrostates(g)
        g.predict_terminal_states()
        a.obs["macrostate"] = g.macrostates.astype(object).to_numpy()
        ms = a.obs.dropna(subset=["macrostate"])
        mt = ms.groupby("macrostate").agg(cells=("label", "size"),
                                         ctr=("genotype", lambda x: int((x == "ctr").sum())),
                                         KO=("genotype", lambda x: int((x == "KO").sum())),
                                         main_label=("label", lambda x: x.value_counts().index[0]))
        mt["terminal"] = mt.index.isin(g.terminal_states.cat.categories)
        mt["search_max"] = hi
        mt.to_csv(TR / f"cellrank_{name}_macrostates.csv")
        print(name, f"data-driven macrostates ({len(mt)}, search 4-{hi}):\n", mt.to_string())
    except Exception as e:   # keep going to the fixed terminal states
        (TR / f"cellrank_{name}_macrostates.csv").unlink(missing_ok=True)
        print(name, "macrostates failed:", e)
    ends = ENDS[cells]
    term = a.obs.label.astype(str).where(a.obs.label.astype(str).isin(ends))
    g.set_terminal_states(states=term.astype("category"), cluster_key="label")
    g.compute_fate_probabilities(n_jobs=NCPU, use_petsc=False)   # petsc4py in the env lacks Mat.getVecs (CellRank 2.0.7); scipy solver instead
    fp = g.fate_probabilities
    fates = pd.DataFrame(np.asarray(fp), index=a.obs_names, columns=list(fp.names))
    out = pd.concat([a.obs[["label", "genotype", "macrostate"]], fates.round(4)], axis=1)
    out.to_csv(TR / f"cellrank_{name}_fates.csv")
    summ = out.groupby(["label", "genotype"], observed=True)[list(fates.columns)].mean().round(3)
    summ.insert(0, "cells", out.groupby(["label", "genotype"], observed=True).size())
    summ.to_csv(TR / f"cellrank_{name}_fate_summary.csv")
    print(name, "mean fate probability by label x genotype:\n", summ.to_string())
    lab = a.obs.label.astype(str).to_numpy(); ko = (a.obs.genotype == "KO").to_numpy()
    early = np.isin(lab, EARLY)
    if cells == "all":
        dr = pd.concat([drivers(full.X, full.var_names, fates, early & ~ko, "control early", DRIVE),
                        drivers(full.X, full.var_names, fates, early, "all early", DRIVE)])
        show = [(l, "control early") for l in DRIVE]
    else:
        tag = f"{cells} early"
        dr = drivers(full.X, full.var_names, fates, early, tag, list(fates.columns))
        show = [(l, tag) for l in fates.columns]
    dr = dr.sort_values(["lineage", "cells", "r"], ascending=[True, True, False], kind="stable")
    dr.to_csv(TR / f"cellrank_{name}_drivers.csv", index=False)
    for l, t in show:
        top = dr[(dr.lineage == l) & (dr.cells == t)].head(20).gene.tolist()
        print(name, l, f"top drivers ({t} cells):", ", ".join(top))
    u = a.obsm["X_umap_2d"]
    panels = [("control", ~ko), ("ISL1 KO", ko)] if cells == "all" else \
             [("control" if cells == "ctr" else "ISL1 KO", np.ones(a.n_obs, bool))]
    nc = len(panels)
    fig, axs = plt.subplots(len(fates.columns), nc, figsize=(5.5 * nc, 4.6 * len(fates.columns)), squeeze=False)
    for i, l in enumerate(fates.columns):
        for j, (t, m) in enumerate(panels):
            ax = axs[i, j]
            ax.scatter(u[~m, 0], u[~m, 1], c="#f0f0f0", s=2, linewidths=0, rasterized=True)
            v = fates[l].to_numpy(); o = np.where(m)[0]; o = o[np.argsort(v[o], kind="stable")]
            s_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=3, cmap=FEAT_CMAP, vmin=0, vmax=1, linewidths=0, rasterized=True)
            ax.set_title(f"{name}: fate -> {l}, {t}", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
            if j == nc - 1:
                plt.colorbar(s_, ax=ax, fraction=0.04)
    fig.tight_layout(); fig.savefig(TR / f"cellrank_{name}_fates.png", dpi=80); plt.close(fig)


def kernels(a):
    ck = cr.kernels.ConnectivityKernel(a).compute_transition_matrix()
    vk = cr.kernels.VelocityKernel(a).compute_transition_matrix()
    pk = cr.kernels.PseudotimeKernel(a, time_key="slingshot_time").compute_transition_matrix()
    return {"velocity": 0.8 * vk + 0.2 * ck, "pseudotime": 0.8 * pk + 0.2 * ck}


def main():
    a = ad.read_h5ad(TR / "scvelo.h5ad")
    a.obs["label"] = a.obs["label"].astype("category")
    a.obs["slingshot_time"] = slingshot_time(a)
    full = ad.read_h5ad(SUB / "D15_endocrine.h5ad")[a.obs_names]
    full = full[:, np.asarray((full.layers["soupx_counts"] > 0).sum(0)).ravel() >= 10].copy()
    print(a, "\ndriver genes tested:", full.n_vars)
    for k, kern in kernels(a).items():
        run(k, a, kern, full, "all")
    for gt in ["ctr", "KO"]:
        b = a[(a.obs.genotype == gt).to_numpy()].copy()
        sc.pp.neighbors(b, use_rep="X_pca", n_neighbors=30, random_state=SEED)   # graph within the genotype
        print(f"\n== {gt} only: {b.n_obs} cells; labels:", b.obs.label.value_counts().to_dict())
        fb = full[b.obs_names]
        for k, kern in kernels(b).items():
            run(f"{k}_{gt}", b, kern, fb, gt)
    print("wrote", TR)


if __name__ == "__main__":
    main()
