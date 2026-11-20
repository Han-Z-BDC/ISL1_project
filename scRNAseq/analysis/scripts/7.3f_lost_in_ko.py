"""Step 7.3f: genes of the control alpha / beta fate programs that KO cells do not switch on where they stop (after 7.3d).

Decisions (user, 2026-11-20): CellRank drivers (7.3d) plus a "lost in KO" comparison; tradeSeq (7.3e) covers the
curve-based genes. Thresholds as proposed and accepted:
  control program  genes with r >= 0.2 and FDR < 0.05 toward the fate in the control early cells (ENP1, ENP2,
                   ENP_FEV_high) of the control-only CellRank runs (7.3d *_ctr_drivers.csv); main list = both kernels
                   (velocity and pseudotime); genes passing in one kernel only are kept and flagged (program = one)
                   beta: fate pre_beta2; alpha: fate pre_alpha
  comparisons      KO vs control at matched points, per lineage:
                   beta   stop    KO pre_beta1  vs control pre_beta1
                   alpha  stop    KO pre_alpha_KO  vs control pre_alpha (the step control cells reach next)
                   both   early   KO ENP_FEV_high vs control ENP_FEV_high (the shared step before the split)
                   per gene: mean log1p expression (7.1d X), log2FC KO / control (of the mean normalised expression,
                   +0.01 pseudocount), % cells expressing, Wilcoxon rank-sum p, FDR (BH, over all tested genes of
                   that comparison); replicate check on per-replicate pseudobulk (summed SoupX counts -> log2 CPM + 1):
                   lost = both KO replicates below both control replicates, gained = both above
  lost in KO       control-program gene with log2FC <= -0.5, FDR < 0.05 and the replicate check passed
  gained in KO     any gene with log2FC >= +0.5, FDR < 0.05 and the replicate check passed (side table; control
                   program membership flagged)
  extra            r toward the same fate in KO early cells in the shared (all-cell) CellRank runs, next to r in
                   control early cells of the same runs
Genes: detected (SoupX counts > 0) in >= 10 of the 5,994 cells (the 7.3d driver gene set).
Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
  cellrank_lost_in_KO.csv     gene x lineage x comparison: program membership and r per kernel, expression stats,
                              replicate pseudobulk, call (lost / gained / -)
Run: sbatch analysis/scripts/7.3f_lost_in_ko.sb   (CELLRANK_ENV python)
"""
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

OUT_DIR = Path(os.environ["OUT_DIR"])
DD = OUT_DIR / "day_analysis" / "D15"
SUB, TR = DD / "subcluster_endocrine", DD / "trajectory"
KERNELS = ["velocity", "pseudotime"]
EARLY = ["ENP1", "ENP2", "ENP_FEV_high"]
R_MIN, FDR_MAX, LFC = 0.2, 0.05, 0.5
LINEAGES = {"beta": "pre_beta2", "alpha": "pre_alpha"}
COMPARE = [("beta", "stop", "pre_beta1", "pre_beta1"),
           ("alpha", "stop", "pre_alpha_KO", "pre_alpha"),
           ("beta", "early", "ENP_FEV_high", "ENP_FEV_high"),
           ("alpha", "early", "ENP_FEV_high", "ENP_FEV_high")]


def dense(m):
    return m.toarray() if sparse.issparse(m) else np.asarray(m)


def corr(X, f):
    Xc = X - X.mean(0); fc = f - f.mean()
    r = (Xc * fc[:, None]).sum(0) / (np.sqrt((Xc ** 2).sum(0)) * np.sqrt((fc ** 2).sum()) + 1e-12)
    return r


def program(genes):
    """Control-program membership and r per kernel (control-only runs, control early cells)."""
    out = pd.DataFrame(index=genes)
    for lin, fate in LINEAGES.items():
        hit = []
        for k in KERNELS:
            d = pd.read_csv(TR / f"cellrank_{k}_ctr_drivers.csv")
            d = d[(d.lineage == fate) & (d.cells == "ctr early")].set_index("gene").reindex(genes)
            out[f"{lin}_r_ctr_{k}"] = d.r
            out[f"{lin}_FDR_ctr_{k}"] = d.FDR
            hit.append((d.r >= R_MIN) & (d.FDR < FDR_MAX))
        n = hit[0].astype(int) + hit[1].astype(int)
        out[f"{lin}_program"] = np.select([n == 2, n == 1], ["both", "one"], "-")
    return out


def main():
    a = ad.read_h5ad(SUB / "D15_endocrine.h5ad")
    keep = np.asarray((a.layers["soupx_counts"] > 0).sum(0)).ravel() >= 10
    a = a[:, keep].copy()
    lab = pd.read_csv(TR / "cellrank_velocity_fates.csv", index_col=0, keep_default_na=False)
    a = a[lab.index].copy()
    a.obs["label"] = lab.label.to_numpy()
    genes = a.var_names
    geno = a.obs.genotype.astype(str).to_numpy(); label = a.obs.label.to_numpy()
    rep = (a.obs.genotype.astype(str) + "_" + a.obs.replicate.astype(str)).to_numpy()
    print(a.n_obs, "cells,", len(genes), "genes; replicates:", pd.Series(rep).value_counts().to_dict())
    prog = program(genes)

    # r toward the fate in KO vs control early cells of the shared runs
    early = np.isin(label, EARLY)
    Xe = {g: dense(a.X[early & (geno == g)]) for g in ["ctr", "KO"]}
    for k in KERNELS:
        f = pd.read_csv(TR / f"cellrank_{k}_fates.csv", index_col=0, keep_default_na=False).reindex(a.obs_names)
        for lin, fate in LINEAGES.items():
            for g in ["ctr", "KO"]:
                prog[f"{lin}_r_shared_{g}_{k}"] = corr(Xe[g], f[fate].to_numpy()[early & (geno == g)]).round(4)

    rows = []
    for lin, step, ko_lab, ctr_lab in COMPARE:
        mk = (geno == "KO") & (label == ko_lab); mc = (geno == "ctr") & (label == ctr_lab)
        Xk, Xc = dense(a.X[mk]), dense(a.X[mc])
        ek, ec = np.expm1(Xk).mean(0), np.expm1(Xc).mean(0)
        lfc = np.log2((ek + 0.01) / (ec + 0.01))
        p = stats.mannwhitneyu(Xk, Xc, axis=0, alternative="two-sided").pvalue
        p = np.where(np.isnan(p), 1.0, p)
        q = stats.false_discovery_control(p)
        d = pd.DataFrame({"gene": genes, "lineage": lin, "step": step,
                          "KO_cells": f"KO {ko_lab} ({mk.sum()})", "ctr_cells": f"ctr {ctr_lab} ({mc.sum()})",
                          "mean_KO": Xk.mean(0).round(3), "mean_ctr": Xc.mean(0).round(3),
                          "pct_KO": (100 * (Xk > 0).mean(0)).round(1), "pct_ctr": (100 * (Xc > 0).mean(0)).round(1),
                          "log2FC": lfc.round(3), "p": p, "FDR": q})
        # per-replicate pseudobulk, log2 CPM + 1
        pb = {}
        for r_ in sorted(set(rep[mk]) | set(rep[mc])):
            m = (mk | mc) & (rep == r_)
            c = np.asarray(a.layers["soupx_counts"][m].sum(0)).ravel()
            pb[r_] = np.log2(c / max(c.sum(), 1) * 1e6 + 1)
            d[f"pb_{r_}"] = pb[r_].round(3)
            d[f"cells_{r_}"] = int(m.sum())
        ko_r = [v for k_, v in pb.items() if k_.startswith("KO_")]
        ct_r = [v for k_, v in pb.items() if k_.startswith("ctr_")]
        lo = np.minimum.reduce(ct_r) > np.maximum.reduce(ko_r)
        hi = np.maximum.reduce(ct_r) < np.minimum.reduce(ko_r)
        d["rep_lower_in_KO"], d["rep_higher_in_KO"] = lo, hi
        pr = prog.loc[genes, f"{lin}_program"].to_numpy()
        lost = (pr == "both") & (lfc <= -LFC) & (q < FDR_MAX) & lo
        lost_one = (pr == "one") & (lfc <= -LFC) & (q < FDR_MAX) & lo
        gained = (lfc >= LFC) & (q < FDR_MAX) & hi
        d["call"] = np.select([lost, lost_one, gained], ["lost", "lost (one kernel)", "gained"], "-")
        cols = [c for c in prog.columns if c.startswith(lin + "_")]
        d = pd.concat([d.reset_index(drop=True), prog.loc[genes, cols].reset_index(drop=True)], axis=1)
        rows.append(d)
        for call in ["lost", "lost (one kernel)", "gained"]:
            s = d[d.call == call].sort_values("log2FC", ascending=(call != "gained"), kind="stable")
            print(f"\n{lin} {step} ({d.KO_cells[0]} vs {d.ctr_cells[0]}): {call}: {len(s)}")
            if len(s):
                print("  ", ", ".join(f"{g} ({v:+.2f})" for g, v in zip(s.gene.head(40), s.log2FC.head(40))))
    for lin in LINEAGES:
        print(f"\n{lin} control program: both kernels", int((prog[f'{lin}_program'] == 'both').sum()),
              "| one kernel", int((prog[f'{lin}_program'] == 'one').sum()))
    out = pd.concat(rows)
    out = out.sort_values(["lineage", "step", "call", "log2FC"], kind="stable")
    out.to_csv(TR / "cellrank_lost_in_KO.csv", index=False)
    print("\nwrote", TR / "cellrank_lost_in_KO.csv")


if __name__ == "__main__":
    main()
