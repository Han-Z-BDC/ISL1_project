"""Step 7.4d: CellRank 2 on the D15 + D22 endocrine integration: fates, terminal states, driver genes (after 7.4c
scVelo and 7.4a Slingshot).

Decisions (user, 2026-11-19): the 7.3 set on the 7.2 integration (scVI). Same design as the D15 run (7.3d):
  kernels    velocity   0.8 VelocityKernel (7.4c) + 0.2 ConnectivityKernel
             pseudotime 0.8 PseudotimeKernel on Slingshot pseudotime (7.4a, scVI latent; one value per cell =
                        lineage-weight-averaged pseudotime, scaled to 0-1) + 0.2 ConnectivityKernel
             (same 30-neighbour graph on the scVI latent for both)
  states     data-driven: GPCCA macrostates (Schur vectors, n 4-12) and predicted terminal states; when the chosen
             n would split a complex-conjugate eigenvalue pair, the search is retried with the upper bound lowered
             by one (12 -> 11 -> ...); the bound used is in the macrostates table (search_max)
             fixed:       terminal states = the clusters that end the 7.4a Slingshot lineages (D22 SC_beta, SC_alpha,
                          SC_delta, SC_EC1, SC_EC2 by design, plus any free ends Slingshot chose)
             per genotype (ctr / KO runs): the Slingshot ends with >= 50 cells of that genotype; KO also
             D15 pre_alpha_KO (as in 7.3d). alpha_beta_EC: ctr SC_beta, SC_alpha, SC_EC1, SC_EC2;
             KO SC_EC1, SC_EC2, D15 pre_alpha_KO (SC_alpha 2 and SC_beta 1 KO cells)
  cells      (user, 2026-11-20, as 7.3d) all cells (shared graph), and a control-only and a KO-only run: 30-neighbour
             graph rebuilt on that genotype's cells (same scVI latent); velocities from the shared scVelo model (7.4c);
             pseudotime = the shared Slingshot pseudotime. The shared graph lets KO cells borrow control paths; the
             per-genotype data-driven terminal states show where each genotype's flow ends.
  fates      fate probabilities toward the fixed terminal states; summary = mean per label x genotype
  drivers    Pearson r between log1p expression (SoupX counts normalised to 10,000 per cell; genes detected in
             >= 10 cells) and a fate probability, in the early D15 cells (ENP1, ENP2, ENP_FEV_high); p from Fisher z,
             FDR (BH). all: SC_alpha and SC_beta, control early cells and all early cells; ctr / KO: every fixed
             terminal state, that genotype's early cells
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/trajectory/ (reruns overwrite); {run} = {kernel} (all
cells) or {kernel}_ctr / {kernel}_KO:
  cellrank_{run}_fates.csv (with the data-driven macrostate per cell), _fate_summary.csv, _macrostates.csv (cells per
  state, by genotype), _drivers.csv, _fates.png (scVI 2D UMAP)
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

# CellRank 2.0.7 (VelocityKernel) calls np.testing.assert_array_equal(x=..., y=...); numpy >= 2 names these
# actual / desired, so pass them on by position (the env has numpy 2.4.6; left unchanged).
_assert_array_equal = np.testing.assert_array_equal


def _assert_array_equal_xy(*args, x=None, y=None, **kw):
    if x is not None or y is not None:
        args = (x, y) + args
    return _assert_array_equal(*args, **kw)


np.testing.assert_array_equal = _assert_array_equal_xy

OUT_DIR = Path(os.environ["OUT_DIR"])
SUBSET = os.environ.get("INT_SUBSET", "")   # 7.2 subset run (see 7.2a); empty = the full integration
INT = OUT_DIR / "day_analysis" / ("integration_D15_D22_endocrine" + (f"_{SUBSET}" if SUBSET else ""))
TR = INT / "trajectory"
EARLY = ["D15 ENP1", "D15 ENP2", "D15 ENP_FEV_high"]
DRIVE = ["SC_alpha", "SC_beta"]
MIN_END = 50                                  # per-genotype fixed ends: >= 50 cells of that genotype
GENO_EXTRA = {"KO": ["D15 pre_alpha_KO"]}     # where KO alpha-side cells stop (as 7.3d)
SEED, NCPU = 1234, int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", ["#e8e8e8", "#fdbb84", "#e34a33", "#7f0000"])


def slingshot_ends():
    ln = pd.read_csv(TR / "slingshot_scvi_lineages.csv")
    return list(dict.fromkeys(c.split(" -> ")[-1] for c in ln.clusters))


def slingshot_time(a):
    s = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0).reindex(a.obs_names)
    pt = np.nan_to_num(s.filter(like="pt_").to_numpy()); w = np.nan_to_num(s.filter(like="w_").to_numpy())
    t = (pt * w).sum(1) / np.maximum(w.sum(1), 1e-12)
    return t / t.max()


def drivers(expr, genes, fates, mask, tag, lineages=DRIVE):
    rows = []
    X = expr[mask]
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    Xc = X - X.mean(0); sx = np.sqrt((Xc ** 2).sum(0))
    n = mask.sum()
    for l in lineages:
        if l not in fates:
            continue
        f = fates[l].to_numpy()[mask]; fc = f - f.mean()
        r = (Xc * fc[:, None]).sum(0) / (sx * np.sqrt((fc ** 2).sum()) + 1e-12)
        z = np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(n - 3)
        p = 2 * stats.norm.sf(np.abs(z))
        rows.append(pd.DataFrame({"gene": genes, "lineage": l, "cells": tag, "n": n, "r": r.round(4), "p": p,
                                  "FDR": stats.false_discovery_control(p)}))
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


def run(name, a, kernel, full, ends, cells):
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
    lab = a.obs.label.astype(str)
    g.set_terminal_states(states=lab.where(lab.isin(ends)).astype("category"), cluster_key="label")
    g.compute_fate_probabilities(n_jobs=NCPU, use_petsc=False)   # petsc4py in the env lacks Mat.getVecs (CellRank 2.0.7); scipy solver instead
    fp = g.fate_probabilities
    fates = pd.DataFrame(np.asarray(fp), index=a.obs_names, columns=list(fp.names))
    out = pd.concat([a.obs[["label", "stage", "genotype", "macrostate"]], fates.round(4)], axis=1)
    out.to_csv(TR / f"cellrank_{name}_fates.csv")
    summ = out.groupby(["label", "genotype"], observed=True)[list(fates.columns)].mean().round(3)
    summ.insert(0, "cells", out.groupby(["label", "genotype"], observed=True).size())
    summ.to_csv(TR / f"cellrank_{name}_fate_summary.csv")
    print(name, "mean fate probability by label x genotype:\n", summ.to_string())
    ko = (a.obs.genotype == "KO").to_numpy()
    early = lab.isin(EARLY).to_numpy()
    if cells == "all":
        dr = pd.concat([drivers(full.X, full.var_names, fates, early & ~ko, "control early"),
                        drivers(full.X, full.var_names, fates, early, "all early")])
        show = [(l, "control early") for l in DRIVE if l in fates]
    else:
        tag = f"{cells} early"
        dr = drivers(full.X, full.var_names, fates, early, tag, list(fates.columns))
        show = [(l, tag) for l in fates.columns]
    dr = dr.sort_values(["lineage", "cells", "r"], ascending=[True, True, False], kind="stable")
    dr.to_csv(TR / f"cellrank_{name}_drivers.csv", index=False)
    for l, t in show:
        top = dr[(dr.lineage == l) & (dr.cells == t)].head(20).gene.tolist()
        print(name, l, f"top drivers ({t} D15 cells):", ", ".join(top))
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
            s_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=3, cmap=FEAT_CMAP, vmin=0, vmax=1, linewidths=0,
                            rasterized=True)
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
    ends = slingshot_ends()
    print("terminal states (Slingshot lineage ends):", ends)
    full = ad.read_h5ad(INT / "merged_counts.h5ad")[a.obs_names].copy()
    keep = np.asarray((full.X > 0).sum(0)).ravel() >= 10
    sc.pp.normalize_total(full, target_sum=1e4); sc.pp.log1p(full)     # totals over all genes, then filter
    full = full[:, keep].copy()
    print(a, "\ndriver genes tested:", full.n_vars)
    for k, kern in kernels(a).items():
        run(k, a, kern, full, ends, "all")
    for gt in ["ctr", "KO"]:
        b = a[(a.obs.genotype == gt).to_numpy()].copy()
        # scanpy 1.11.5 + scipy 1.17.1 (CELLRANK_ENV): the default pynndescent neighbours fail above a few thousand cells
        # (csr_array has no getnnz), so exact neighbours with sklearn (user, 2026-11-21)
        sc.pp.neighbors(b, use_rep="X_scvi", n_neighbors=30, random_state=SEED, transformer="sklearn")   # graph within genotype
        n = b.obs.label.astype(str).value_counts()
        e = [x for x in ends if n.get(x, 0) >= MIN_END] + [x for x in GENO_EXTRA.get(gt, []) if n.get(x, 0) >= MIN_END]
        print(f"\n== {gt} only: {b.n_obs} cells; fixed ends: {e}; labels:", n.to_dict())
        fb = full[b.obs_names]
        for k, kern in kernels(b).items():
            run(f"{k}_{gt}", b, kern, fb, list(dict.fromkeys(e)), gt)
    print("wrote", TR)


if __name__ == "__main__":
    main()
