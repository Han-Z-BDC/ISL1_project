"""Step 7.3f: genes of the control alpha / beta fate programs that KO cells do not switch on where they stop (after 7.3d).

Decisions (user, 2026-11-20): CellRank drivers (7.3d) plus a "lost in KO" comparison; tradeSeq (7.3e) covers the
curve-based genes. Thresholds as proposed and accepted:
  control program  genes with r >= 0.2 and FDR < 0.05 toward the fate in the control early cells (ENP1, ENP2,
                   ENP_FEV_high) of the control-only CellRank runs (7.3d *_ctr_drivers.csv); main list = both kernels
                   (velocity and pseudotime); genes passing in one kernel only are kept and flagged (program = one)
                   beta: fate pre_beta2; alpha: fate pre_alpha
  comparisons      KO vs control, per lineage, in two ways (user, 2026-11-20: both):
                   by label (match = "label")
                     beta   stop    KO pre_beta1  vs control pre_beta1
                     alpha  stop    KO pre_alpha_KO  vs control pre_alpha (the step control cells reach next)
                     both   early   KO ENP_FEV_high vs control ENP_FEV_high (the shared step before the split)
                   by Slingshot pseudotime (match = "pseudotime"; 7.3a PCs 1-30, the shared trajectory)
                     cells on the lineage = curve weight >= 0.5 (beta: the lineage ending in pre_beta2, alpha: the
                     lineage ending in pre_alpha); stop window = the KO cells' 75th to 95th percentile of pseudotime on
                     that lineage; KO vs control cells on the lineage inside the window (any label)
                   per gene: mean log1p expression (7.1d X), log2FC KO / control (of the mean normalised expression,
                   +0.01 pseudocount), % cells expressing, Wilcoxon rank-sum p, FDR (BH, over all tested genes of
                   that comparison); replicate check on per-replicate pseudobulk (summed SoupX counts -> log2 CPM + 1):
                   lost = both KO replicates below both control replicates, gained = both above
  lost in KO       control-program gene with log2FC <= -0.5, FDR < 0.05 and the replicate check passed
  gained in KO     any gene with log2FC >= +0.5, FDR < 0.05 and the replicate check passed (side table; control
                   program membership flagged)
  extra            r toward the same fate in KO early cells in the shared (all-cell) CellRank runs, next to r in
                   control early cells of the same runs
  profile          mean expression of the program genes (both + one kernel) along each lineage, control and KO, in
                   12 equal-width pseudotime bins (cells on the lineage), for later plots
  on after stop    (user, 2026-11-20) genes control cells switch on right after the point where KO stops, per
                   lineage: control cells on the lineage in the first half of the control path past the stop window
                   (pseudotime from the window's end to the median of the control cells beyond it) vs control cells in
                   the stop window; rises = log2FC >= 0.5, FDR < 0.05, and higher after the stop in each control
                   replicate (pseudobulk); missing in KO = rises and KO cells in the stop window are below the control
                   cells after it (log2FC <= -0.5, FDR < 0.05, both KO replicates below both control replicates);
                   all genes tested, control-program membership flagged
Genes: detected (SoupX counts > 0) in >= 10 of the 5,994 cells (the 7.3d driver gene set).
Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
  cellrank_lost_in_KO.csv          gene x lineage x match x comparison: program membership and r per kernel,
                                   expression stats, replicate pseudobulk, call (lost / lost (one kernel) / gained / -)
  cellrank_lost_in_KO_profile.csv  gene, lineage, genotype, bin, pseudotime (bin centre), cells, mean log1p
  cellrank_on_after_KO_stop.csv    gene x lineage: control stop / after and KO stop expression, rise in control,
                                   KO vs control after, replicate pseudobulk, call
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
W_MIN, Q_LO, Q_HI, NBIN = 0.5, 0.75, 0.95, 12
LINEAGES = {"beta": "pre_beta2", "alpha": "pre_alpha"}
COMPARE = [("beta", "stop", "pre_beta1", "pre_beta1"),
           ("alpha", "stop", "pre_alpha_KO", "pre_alpha"),
           ("beta", "early", "ENP_FEV_high", "ENP_FEV_high"),
           ("alpha", "early", "ENP_FEV_high", "ENP_FEV_high")]


def dense(m):
    return m.toarray() if sparse.issparse(m) else np.asarray(m)


def corr(X, f):
    Xc = X - X.mean(0); fc = f - f.mean()
    return (Xc * fc[:, None]).sum(0) / (np.sqrt((Xc ** 2).sum(0)) * np.sqrt((fc ** 2).sum()) + 1e-12)


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


def compare(a, mk, mc, rep, prog, lin, match, step, ko_desc, ctr_desc):
    """KO (mk) vs control (mc) cells: expression stats, replicate pseudobulk, call."""
    genes = a.var_names
    Xk, Xc = dense(a.X[mk]), dense(a.X[mc])
    ek, ec = np.expm1(Xk).mean(0), np.expm1(Xc).mean(0)
    lfc = np.log2((ek + 0.01) / (ec + 0.01))
    p = stats.mannwhitneyu(Xk, Xc, axis=0, alternative="two-sided").pvalue
    p = np.where(np.isnan(p), 1.0, p)
    q = stats.false_discovery_control(p)
    d = pd.DataFrame({"gene": genes, "lineage": lin, "match": match, "step": step,
                      "KO_cells": f"{ko_desc} ({mk.sum()})", "ctr_cells": f"{ctr_desc} ({mc.sum()})",
                      "mean_KO": Xk.mean(0).round(3), "mean_ctr": Xc.mean(0).round(3),
                      "pct_KO": (100 * (Xk > 0).mean(0)).round(1), "pct_ctr": (100 * (Xc > 0).mean(0)).round(1),
                      "log2FC": lfc.round(3), "p": p, "FDR": q})
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
    for call in ["lost", "lost (one kernel)", "gained"]:
        s = d[d.call == call].sort_values("log2FC", ascending=(call != "gained"), kind="stable")
        print(f"\n{lin} {match} {step} ({d.KO_cells[0]} vs {d.ctr_cells[0]}): {call}: {len(s)}")
        if len(s):
            print("  ", ", ".join(f"{g} ({v:+.2f})" for g, v in zip(s.gene.head(40), s.log2FC.head(40))))
    return d


def wilcox_lfc(X1, X0):
    """log2FC (mean normalised expression, +0.01) and Wilcoxon FDR of X1 vs X0 (cells x genes, log1p)."""
    lfc = np.log2((np.expm1(X1).mean(0) + 0.01) / (np.expm1(X0).mean(0) + 0.01))
    p = stats.mannwhitneyu(X1, X0, axis=0, alternative="two-sided").pvalue
    p = np.where(np.isnan(p), 1.0, p)
    return lfc, stats.false_discovery_control(p)


def pseudobulk(a, m):
    c = np.asarray(a.layers["soupx_counts"][m].sum(0)).ravel()
    return np.log2(c / max(c.sum(), 1) * 1e6 + 1)


def on_after(a, pt, on, geno, label, rep, prog, lin, L, lo_, hi_):
    """Genes control cells switch on right after the KO stop window, and whether KO cells express them."""
    genes = a.var_names
    beyond = on & (geno == "ctr") & (pt > hi_)
    mid = np.median(pt[beyond])
    ca = beyond & (pt <= mid)                                   # control, first half of the path past the KO stop
    cs = on & (geno == "ctr") & (pt >= lo_) & (pt <= hi_)        # control, stop window
    ks = on & (geno == "KO") & (pt >= lo_) & (pt <= hi_)         # KO, stop window
    kb = on & (geno == "KO") & (pt > hi_)                        # the few KO cells past the stop
    print(f"\n== {lin}: control after the KO stop, pseudotime {hi_:.1f}-{mid:.1f}: {ca.sum()} cells; labels:",
          pd.Series(label[ca]).value_counts().to_dict(), f"| KO cells past the stop: {kb.sum()}")
    Xa, Xs, Xk = dense(a.X[ca]), dense(a.X[cs]), dense(a.X[ks])
    lfc_up, q_up = wilcox_lfc(Xa, Xs)                            # rise in control after vs at the stop
    lfc_ko, q_ko = wilcox_lfc(Xk, Xa)                            # KO at the stop vs control after
    reps = sorted(set(rep[ca | cs]))
    pa = {r_: pseudobulk(a, ca & (rep == r_)) for r_ in reps}
    ps = {r_: pseudobulk(a, cs & (rep == r_)) for r_ in reps}
    rise_rep = np.logical_and.reduce([pa[r_] > ps[r_] for r_ in reps])
    pk = {r_: pseudobulk(a, ks & (rep == r_)) for r_ in sorted(set(rep[ks]))}
    ko_low = np.maximum.reduce(list(pk.values())) < np.minimum.reduce(list(pa.values()))
    rises = (lfc_up >= LFC) & (q_up < FDR_MAX) & rise_rep
    missing = rises & (lfc_ko <= -LFC) & (q_ko < FDR_MAX) & ko_low
    d = pd.DataFrame({"gene": genes, "lineage": lin, "slingshot_lineage": L,
                      "window_stop": f"{lo_:.1f}-{hi_:.1f}", "window_after": f"{hi_:.1f}-{mid:.1f}",
                      "mean_ctr_stop": Xs.mean(0).round(3), "mean_ctr_after": Xa.mean(0).round(3),
                      "mean_KO_stop": Xk.mean(0).round(3),
                      "pct_ctr_stop": (100 * (Xs > 0).mean(0)).round(1), "pct_ctr_after": (100 * (Xa > 0).mean(0)).round(1),
                      "pct_KO_stop": (100 * (Xk > 0).mean(0)).round(1),
                      "log2FC_ctr_after_vs_stop": lfc_up.round(3), "FDR_ctr_after_vs_stop": q_up,
                      "log2FC_KO_stop_vs_ctr_after": lfc_ko.round(3), "FDR_KO_stop_vs_ctr_after": q_ko,
                      "rep_rise_in_ctr": rise_rep, "rep_KO_below_ctr_after": ko_low,
                      "call": np.select([missing, rises], ["on after stop, missing in KO", "on after stop"], "-")})
    for r_ in reps:
        d[f"pb_after_{r_}"] = pa[r_].round(3); d[f"pb_stop_{r_}"] = ps[r_].round(3)
    for r_, v in pk.items():
        d[f"pb_stop_{r_}"] = v.round(3)
    cols = [c for c in prog.columns if c.startswith(lin + "_")]
    d = pd.concat([d.reset_index(drop=True), prog.loc[genes, cols].reset_index(drop=True)], axis=1)
    d = d.sort_values(["call", "log2FC_ctr_after_vs_stop"], ascending=[True, False], kind="stable")
    for call in ["on after stop, missing in KO", "on after stop"]:
        s = d[d.call == call]
        print(f"{lin}: {call}: {len(s)} (control program both kernels: "
              f"{int((s[f'{lin}_program'] == 'both').sum())})")
        print("  ", ", ".join(f"{g} ({v:+.2f})" for g, v in zip(s.gene.head(50), s.log2FC_ctr_after_vs_stop.head(50))))
    return d


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
    # 1. matched by label
    for lin, step, ko_lab, ctr_lab in COMPARE:
        mk = (geno == "KO") & (label == ko_lab); mc = (geno == "ctr") & (label == ctr_lab)
        rows.append(compare(a, mk, mc, rep, prog, lin, "label", step, f"KO {ko_lab}", f"ctr {ctr_lab}"))

    # 2. matched by Slingshot pseudotime on the lineage
    sl = pd.read_csv(TR / "slingshot_pc30_cells.csv", index_col=0).reindex(a.obs_names)
    ln = pd.read_csv(TR / "slingshot_pc30_lineages.csv")
    end = {r.lineage: r.clusters.split(" -> ")[-1] for r in ln.itertuples()}
    prof_genes = np.isin(genes, prog.index[(prog.beta_program != "-") | (prog.alpha_program != "-")])
    profile, after = [], []
    for lin, fate in LINEAGES.items():
        L = [k for k, v in end.items() if v == fate]
        assert len(L) == 1, f"{fate}: lineages {L}"
        L = L[0]
        pt, w = sl[f"pt_{L}"].to_numpy(), sl[f"w_{L}"].fillna(0).to_numpy()
        on = (w >= W_MIN) & ~np.isnan(pt)
        lo_, hi_ = np.quantile(pt[on & (geno == "KO")], [Q_LO, Q_HI])
        win = on & (pt >= lo_) & (pt <= hi_)
        mk, mc = win & (geno == "KO"), win & (geno == "ctr")
        print(f"\n== {lin}: {L} ({' -> '.join(ln.set_index('lineage').clusters[L].split(' -> ')[-2:])}); "
              f"stop window pseudotime {lo_:.1f}-{hi_:.1f} (KO {Q_LO:.0%}-{Q_HI:.0%})")
        for g, m in [("KO", mk), ("ctr", mc)]:
            print(f"   {g} cells in window: {m.sum()}; labels:", pd.Series(label[m]).value_counts().to_dict())
        rows.append(compare(a, mk, mc, rep, prog, lin, "pseudotime", "stop",
                            f"KO {L} pt {lo_:.0f}-{hi_:.0f}", f"ctr {L} pt {lo_:.0f}-{hi_:.0f}"))
        after.append(on_after(a, pt, on, geno, label, rep, prog, lin, L, lo_, hi_))
        # profile along the lineage
        edges = np.linspace(np.nanmin(pt[on]), np.nanmax(pt[on]), NBIN + 1)
        b = np.clip(np.digitize(pt, edges) - 1, 0, NBIN - 1)
        Xp = a.X[:, prof_genes]
        for g in ["ctr", "KO"]:
            for i in range(NBIN):
                m = on & (geno == g) & (b == i)
                if m.sum() == 0:
                    continue
                profile.append(pd.DataFrame({"gene": genes[prof_genes], "lineage": lin, "genotype": g, "bin": i + 1,
                                             "pseudotime": round((edges[i] + edges[i + 1]) / 2, 1),
                                             "cells": int(m.sum()),
                                             "mean": np.asarray(Xp[m].mean(0)).ravel().round(4)}))
    for lin in LINEAGES:
        print(f"\n{lin} control program: both kernels", int((prog[f'{lin}_program'] == 'both').sum()),
              "| one kernel", int((prog[f'{lin}_program'] == 'one').sum()))
    out = pd.concat(rows)
    out = out.sort_values(["lineage", "match", "step", "call", "log2FC"], kind="stable")
    out.to_csv(TR / "cellrank_lost_in_KO.csv", index=False)
    pd.concat(profile).to_csv(TR / "cellrank_lost_in_KO_profile.csv", index=False)
    pd.concat(after).to_csv(TR / "cellrank_on_after_KO_stop.csv", index=False)
    # agreement of the two matchings at the stop
    for lin in LINEAGES:
        s = {m: set(out[(out.lineage == lin) & (out.match == m) & (out.step == "stop") & (out.call == "lost")].gene)
             for m in ["label", "pseudotime"]}
        print(f"\n{lin} lost at the stop (both kernels): label {len(s['label'])}, pseudotime {len(s['pseudotime'])}, "
              f"both {len(s['label'] & s['pseudotime'])}: {', '.join(sorted(s['label'] & s['pseudotime']))}")
        print(f"  pseudotime only: {', '.join(sorted(s['pseudotime'] - s['label']))}")
    print("\nwrote", TR / "cellrank_lost_in_KO.csv", "and _profile.csv")


if __name__ == "__main__":
    main()
