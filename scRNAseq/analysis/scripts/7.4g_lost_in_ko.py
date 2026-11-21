"""Step 7.4g: D15 + D22 (alpha_beta_EC subset): KO vs control where KO stops on the beta and alpha lineages, and genes
control switches on right after the stop (the 7.3f trajectory version on the 7.4 trajectory; after 7.4a, 7.4d).

Decisions (user, 2026-11-21): repeat the D15 trajectory-based KO vs control analysis (7.3f pseudotime match, 7.3g
windows, 7.3h patterns) on the D15 + D22 integration. Changes from D15, all following from the data:
  on a lineage   curve weight >= W_MIN AND the cell's label is on that lineage's Slingshot path (7.4a, scVI latent).
                 Weight alone puts many D22 EC cells on the beta / alpha curves (e.g. beta: KO 6,218 with weight >= 0.5,
                 2,243 on the path)
  programs       CellRank drivers of the pseudotime kernel only (7.4d *_ctr_drivers.csv, control early cells; r >= 0.2,
                 FDR < 0.05): the velocity runs backwards on D15 + D22 too (DECISIONS 7.4d)
  stage          D15 and D22 are separate sequencing batches, so each comparison stays within one stage: the stop
                 windows are D15 cells; "on after the stop" compares control cells of the same stage (D15) after vs at
                 the stop
Otherwise as 7.3f (label matches not used): stop window = KO cells' 75th-95th percentile of pseudotime on the lineage;
KO vs control cells in the window (same stage): log2FC (mean normalised expression +0.01), Wilcoxon, FDR (BH), per-
replicate pseudobulk; lost = control-program gene, log2FC <= -0.5, FDR < 0.05, both KO replicates below both control
replicates; gained = any gene, the same upward. On after the stop: control after (first half of the path past the window)
vs control in the window: rises = log2FC >= 0.5, FDR < 0.05, up in each control replicate; missing in KO = rises and KO
cells in the window below control after (log2FC <= -0.5, FDR < 0.05, both KO replicates below both control replicates).
Every row carries the control beta / alpha / EC1 / EC2 program flags and r.
Expression: SoupX counts (7.2a merged_counts.h5ad of the subset) normalised to 10,000 per cell, log1p; genes detected in
>= 10 cells.
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine_<subset>/trajectory/ (reruns overwrite):
  cellrank_lost_in_KO.csv, cellrank_on_after_KO_stop.csv
Run: sbatch analysis/scripts/7.4g_lost_in_ko.sb   (CELLRANK_ENV python; INT_SUBSET, default alpha_beta_EC)
"""
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats

OUT_DIR = Path(os.environ["OUT_DIR"])
SUBSET = os.environ.get("INT_SUBSET", "alpha_beta_EC")
INT = OUT_DIR / "day_analysis" / f"integration_D15_D22_endocrine_{SUBSET}"
TR = INT / "trajectory"
R_MIN, FDR_MAX, LFC, W_MIN, Q_LO, Q_HI = 0.2, 0.05, 0.5, 0.5, 0.75, 0.95
LINEAGES = {"beta": "SC_beta", "alpha": "SC_alpha"}
PROGRAMS = {"beta": "SC_beta", "alpha": "SC_alpha", "EC1": "SC_EC1", "EC2": "SC_EC2"}


def dense(m):
    return m.toarray() if sparse.issparse(m) else np.asarray(m)


def load():
    a = ad.read_h5ad(INT / "merged_counts.h5ad")
    keep = np.asarray((a.X > 0).sum(0)).ravel() >= 10
    a.layers["soupx_counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    a = a[:, keep].copy()
    sl = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0)
    a = a[sl.index].copy()
    a.obs["label"] = sl.label.to_numpy()
    a.obs["rep"] = (a.obs.genotype.astype(str) + "_" + a.obs["sample"].astype(str)).to_numpy()
    return a, sl


def programs(genes):
    out = pd.DataFrame(index=genes)
    d = pd.read_csv(TR / "cellrank_pseudotime_ctr_drivers.csv")
    for p, fate in PROGRAMS.items():
        x = d[(d.lineage == fate) & (d.cells == "ctr early")].set_index("gene").reindex(genes)
        out[f"{p}_r_ctr"], out[f"{p}_FDR_ctr"] = x.r.round(4), x.FDR
        out[f"{p}_program"] = np.where((x.r >= R_MIN) & (x.FDR < FDR_MAX), "yes", "-")
    return out


def on_lineage(sl, L, path):
    return ((sl[f"w_{L}"].fillna(0) >= W_MIN) & sl.label.isin(path) & sl[f"pt_{L}"].notna()).to_numpy()


def pseudobulk(a, m):
    c = np.asarray(a.layers["soupx_counts"][m].sum(0)).ravel()
    return np.log2(c / max(c.sum(), 1) * 1e6 + 1)


def wilcox_lfc(X1, X0):
    lfc = np.log2((np.expm1(X1).mean(0) + 0.01) / (np.expm1(X0).mean(0) + 0.01))
    p = stats.mannwhitneyu(X1, X0, axis=0, alternative="two-sided").pvalue
    p = np.where(np.isnan(p), 1.0, p)
    return lfc, stats.false_discovery_control(p)


def main():
    a, sl = load()
    genes = a.var_names
    geno, rep = a.obs.genotype.astype(str).to_numpy(), a.obs.rep.to_numpy()
    stage, label = a.obs.stage.astype(str).to_numpy(), a.obs.label.to_numpy()
    print(a.n_obs, "cells,", len(genes), "genes; replicates:", pd.Series(rep).value_counts().to_dict())
    prog = programs(genes)
    for p in PROGRAMS:
        print(f"control {p} program (pseudotime kernel): {int((prog[f'{p}_program'] == 'yes').sum())} genes")
    pcols = list(prog.columns)
    ln = pd.read_csv(TR / "slingshot_scvi_lineages.csv")
    paths = {r.lineage: r.clusters.split(" -> ") for r in ln.itertuples()}
    rows, after = [], []
    for lin, fate in LINEAGES.items():
        L = [k for k, v in paths.items() if v[-1] == fate][0]
        pt = sl[f"pt_{L}"].reindex(a.obs_names).to_numpy()
        on = on_lineage(sl.reindex(a.obs_names), L, paths[L])
        lo_, hi_ = np.quantile(pt[on & (geno == "KO")], [Q_LO, Q_HI])
        win = on & (pt >= lo_) & (pt <= hi_)
        st = pd.Series(stage[win]).value_counts().index[0]
        win &= stage == st
        mk, mc = win & (geno == "KO"), win & (geno == "ctr")
        print(f"\n== {lin}: {L} ({' -> '.join(paths[L])}); stop window pseudotime {lo_:.2f}-{hi_:.2f} (KO "
              f"{Q_LO:.0%}-{Q_HI:.0%}), stage {st}")
        for g, m in [("KO", mk), ("ctr", mc)]:
            print(f"   {g} cells in window: {m.sum()}; labels:", pd.Series(label[m]).value_counts().to_dict())
        # KO vs control at the stop
        Xk, Xc = dense(a.X[mk]), dense(a.X[mc])
        lfc, q = wilcox_lfc(Xk, Xc)
        d = pd.DataFrame({"gene": genes, "lineage": lin, "slingshot_lineage": L, "stage": st,
                          "window": f"{lo_:.2f}-{hi_:.2f}", "KO_cells": int(mk.sum()), "ctr_cells": int(mc.sum()),
                          "mean_KO": Xk.mean(0).round(3), "mean_ctr": Xc.mean(0).round(3),
                          "pct_KO": (100 * (Xk > 0).mean(0)).round(1), "pct_ctr": (100 * (Xc > 0).mean(0)).round(1),
                          "log2FC": lfc.round(3), "FDR": q})
        pb = {r_: pseudobulk(a, (mk | mc) & (rep == r_)) for r_ in sorted(set(rep[mk | mc]))}
        for r_, v in pb.items():
            d[f"pb_{r_}"] = v.round(3)
        ko_r = [v for k_, v in pb.items() if k_.startswith("KO_")]; ct_r = [v for k_, v in pb.items() if k_.startswith("ctr_")]
        lo_rep = np.minimum.reduce(ct_r) > np.maximum.reduce(ko_r)
        hi_rep = np.maximum.reduce(ct_r) < np.minimum.reduce(ko_r)
        inprog = (prog[f"{lin}_program"] == "yes").to_numpy()
        d["call"] = np.select([inprog & (lfc <= -LFC) & (q < FDR_MAX) & lo_rep, (lfc >= LFC) & (q < FDR_MAX) & hi_rep],
                              ["lost", "gained"], "-")
        d = pd.concat([d.reset_index(drop=True), prog.reset_index(drop=True)], axis=1)
        rows.append(d)
        for call in ["lost", "gained"]:
            s = d[d.call == call].sort_values("log2FC", ascending=(call == "lost"))
            print(f"   {call}: {len(s)}:", ", ".join(s.gene.head(40)))
        # control after the stop (same stage) vs control at the stop; KO at the stop vs control after
        beyond = on & (geno == "ctr") & (pt > hi_) & (stage == st)
        mid = np.median(pt[beyond])
        ca = beyond & (pt <= mid)
        print(f"   control after the stop (stage {st}), pseudotime {hi_:.2f}-{mid:.2f}: {ca.sum()} cells;",
              pd.Series(label[ca]).value_counts().to_dict())
        Xa = dense(a.X[ca])
        lfc_up, q_up = wilcox_lfc(Xa, Xc)
        lfc_ko, q_ko = wilcox_lfc(Xk, Xa)
        reps_c = sorted(set(rep[ca | mc]))
        pa = {r_: pseudobulk(a, ca & (rep == r_)) for r_ in reps_c}
        ps = {r_: pseudobulk(a, mc & (rep == r_)) for r_ in reps_c}
        pk = {r_: pseudobulk(a, mk & (rep == r_)) for r_ in sorted(set(rep[mk]))}
        rise_rep = np.logical_and.reduce([pa[r_] > ps[r_] for r_ in reps_c])
        ko_low = np.maximum.reduce(list(pk.values())) < np.minimum.reduce(list(pa.values()))
        rises = (lfc_up >= LFC) & (q_up < FDR_MAX) & rise_rep
        missing = rises & (lfc_ko <= -LFC) & (q_ko < FDR_MAX) & ko_low
        e = pd.DataFrame({"gene": genes, "lineage": lin, "slingshot_lineage": L, "stage": st,
                          "window_stop": f"{lo_:.2f}-{hi_:.2f}", "window_after": f"{hi_:.2f}-{mid:.2f}",
                          "mean_ctr_stop": Xc.mean(0).round(3), "mean_ctr_after": Xa.mean(0).round(3),
                          "mean_KO_stop": Xk.mean(0).round(3), "pct_ctr_stop": (100 * (Xc > 0).mean(0)).round(1),
                          "pct_ctr_after": (100 * (Xa > 0).mean(0)).round(1), "pct_KO_stop": (100 * (Xk > 0).mean(0)).round(1),
                          "log2FC_ctr_after_vs_stop": lfc_up.round(3), "FDR_ctr_after_vs_stop": q_up,
                          "log2FC_KO_stop_vs_ctr_after": lfc_ko.round(3), "FDR_KO_stop_vs_ctr_after": q_ko,
                          "rep_rise_in_ctr": rise_rep, "rep_KO_below_ctr_after": ko_low,
                          "call": np.select([missing, rises], ["on after stop, missing in KO", "on after stop"], "-")})
        e = pd.concat([e.reset_index(drop=True), prog.reset_index(drop=True)], axis=1)
        after.append(e)
        for call in ["on after stop, missing in KO", "on after stop"]:
            s = e[e.call == call].sort_values("log2FC_ctr_after_vs_stop", ascending=False)
            print(f"   {call}: {len(s)} (in the control {lin} program: {int((s[f'{lin}_program'] == 'yes').sum())}):",
                  ", ".join(s.gene.head(40)))
    pd.concat(rows).to_csv(TR / "cellrank_lost_in_KO.csv", index=False)
    pd.concat(after).to_csv(TR / "cellrank_on_after_KO_stop.csv", index=False)
    print("\nwrote", TR / "cellrank_lost_in_KO.csv", "and cellrank_on_after_KO_stop.csv")


if __name__ == "__main__":
    main()
