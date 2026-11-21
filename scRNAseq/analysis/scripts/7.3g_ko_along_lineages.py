"""Step 7.3g: KO vs control along the beta, alpha and EC lineages in matched pseudotime windows (after 7.3a, 7.3d, 7.3f).

Decisions (user, 2026-11-21): the control vs KO comparison uses the trajectory (Slingshot pseudotime) version, and the
pre_EC lineage is added. On the EC lineage KO cells do not stop (they reach further than control: median pseudotime
143 vs 137, 95th percentile 184 vs 177; 712 KO vs 233 control pre_EC cells), so the KO stop window of 7.3f does not
apply; instead each lineage is compared along its whole length. 7.3f (stop window, on after the stop) is unchanged
apart from an EC-program flag.
  lineages   beta (Slingshot lineage ending in pre_beta2), alpha (pre_alpha), EC (pre_EC; ENP1 -> ENP2 ->
             ENP_FEV_high -> pre_beta1 -> pre_EC); 7.3a PCs 1-30, one shared trajectory for both genotypes; a cell is
             on a lineage at curve weight >= W_MIN
  windows    NWIN windows of equal cell numbers along each lineage (pseudotime quantiles of the cells on it); a window
             is tested when it has >= MIN_CELLS cells of each genotype; the replicate check needs >= MIN_REP cells in
             every replicate. Part = "shared path" when >= half of the window's cells are also on another of the three
             lineages (curve weight >= W_MIN), else "branch"; only branch windows are specific to one lineage
  test       KO vs control cells in the window, all genes detected in >= 10 cells (as 7.3d / 7.3f): mean log1p,
             % > 0, log2FC of the mean normalised expression (+0.01), Wilcoxon rank-sum, FDR (BH, per window); per-
             replicate pseudobulk (summed SoupX counts -> log2 CPM + 1). higher in KO = log2FC >= LFC, FDR < 0.05 and
             both KO replicates above both control replicates; lower in KO = the same downward (as 7.3f)
  programs   control fate programs from the control-only CellRank runs (7.3d *_ctr_drivers.csv, control early cells):
             r >= 0.2 and FDR < 0.05 toward pre_beta2 / pre_alpha / pre_EC; both = in both kernels, one = in one.
             KO EC program from the KO-only runs (*_KO_drivers.csv, KO early cells, toward pre_EC), same rule.
  readouts   per window and genotype: cells, replicates, label mix, % KO, mean CellRank fate toward pre_EC / pre_beta2 /
             pre_alpha (shared runs, both kernels), and program scores (mean log1p of each control program's
             both-kernel genes): is the EC program on earlier / stronger in KO, do KO cells on the EC branch keep
             beta-program genes, where KO cells pile up
  profile    mean log1p of all program genes (beta, alpha, EC control programs; KO EC program) along each lineage,
             12 equal-width pseudotime bins, per genotype
Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
  ko_windows.csv              lineage x window x genotype: pseudotime range, part, cells (per replicate), labels, % KO,
                              mean fates, program scores, tested
  ko_vs_ctr_windows.csv.gz    gene x lineage x window: stats, replicate pseudobulk, call, program flags (beta / alpha /
                              EC control, EC KO)
  ec_programs.csv             gene: r / FDR toward pre_EC in control and KO early cells, per kernel; control and KO
                              program membership; overlap (both / control only / KO only)
  ko_lineage_profile.csv      gene, lineage, genotype, bin, pseudotime, cells, mean log1p
Run: sbatch analysis/scripts/7.3g_ko_along_lineages.sb   (CELLRANK_ENV python)
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
R_MIN, FDR_MAX, LFC = 0.2, 0.05, 0.5
W_MIN, NWIN, MIN_CELLS, MIN_REP, NBIN = 0.5, 8, 30, 5, 12
LINEAGES = {"beta": "pre_beta2", "alpha": "pre_alpha", "EC": "pre_EC"}
FATES = ["pre_EC", "pre_beta2", "pre_alpha"]


def dense(m):
    return m.toarray() if sparse.issparse(m) else np.asarray(m)


def drivers(geno, fate, genes):
    """r / FDR toward a fate in one genotype's early cells (7.3d per-genotype runs), both kernels; membership."""
    out = pd.DataFrame(index=genes)
    hit = []
    for k in KERNELS:
        d = pd.read_csv(TR / f"cellrank_{k}_{geno}_drivers.csv")
        d = d[(d.lineage == fate) & (d.cells == f"{geno} early")].set_index("gene").reindex(genes)
        assert d.r.notna().any(), f"no {geno} drivers toward {fate} ({k})"
        out[f"r_{geno}_{k}"], out[f"FDR_{geno}_{k}"] = d.r, d.FDR
        hit.append(((d.r >= R_MIN) & (d.FDR < FDR_MAX)).to_numpy())
    n = hit[0].astype(int) + hit[1].astype(int)
    out["program"] = np.select([n == 2, n == 1], ["both", "one"], "-")
    return out


def pseudobulk(a, m):
    c = np.asarray(a.layers["soupx_counts"][m].sum(0)).ravel()
    return np.log2(c / max(c.sum(), 1) * 1e6 + 1)


def test_window(a, mk, mc, rep):
    Xk, Xc = dense(a.X[mk]), dense(a.X[mc])
    lfc = np.log2((np.expm1(Xk).mean(0) + 0.01) / (np.expm1(Xc).mean(0) + 0.01))
    p = stats.mannwhitneyu(Xk, Xc, axis=0, alternative="two-sided").pvalue
    p = np.where(np.isnan(p), 1.0, p)
    q = stats.false_discovery_control(p)
    d = pd.DataFrame({"mean_KO": Xk.mean(0).round(3), "mean_ctr": Xc.mean(0).round(3),
                      "pct_KO": (100 * (Xk > 0).mean(0)).round(1), "pct_ctr": (100 * (Xc > 0).mean(0)).round(1),
                      "log2FC": lfc.round(3), "p": p, "FDR": q}, index=a.var_names)
    reps = sorted(set(rep[mk | mc]))
    n_rep = {r_: int(((mk | mc) & (rep == r_)).sum()) for r_ in reps}
    rep_ok = len([r_ for r_ in reps if r_.startswith("KO_")]) >= 2 and \
        len([r_ for r_ in reps if r_.startswith("ctr_")]) >= 2 and min(n_rep.values()) >= MIN_REP
    pb = {r_: pseudobulk(a, (mk | mc) & (rep == r_)) for r_ in reps}
    for r_ in reps:
        d[f"pb_{r_}"] = pb[r_].round(3)
    if rep_ok:
        ko_r = [v for k_, v in pb.items() if k_.startswith("KO_")]
        ct_r = [v for k_, v in pb.items() if k_.startswith("ctr_")]
        lo = np.minimum.reduce(ct_r) > np.maximum.reduce(ko_r)
        hi = np.maximum.reduce(ct_r) < np.minimum.reduce(ko_r)
    else:
        lo = hi = np.zeros(len(d), bool)
    d["rep_lower_in_KO"], d["rep_higher_in_KO"] = lo, hi
    d["call"] = np.select([(lfc >= LFC) & (q < FDR_MAX) & hi, (lfc <= -LFC) & (q < FDR_MAX) & lo],
                          ["higher in KO", "lower in KO"], "-")
    return d, rep_ok, n_rep


def main():
    a = ad.read_h5ad(SUB / "D15_endocrine.h5ad")
    keep = np.asarray((a.layers["soupx_counts"] > 0).sum(0)).ravel() >= 10
    a = a[:, keep].copy()
    lab = pd.read_csv(TR / "cellrank_velocity_fates.csv", index_col=0, keep_default_na=False)
    a = a[lab.index].copy()
    genes = a.var_names
    label = lab.label.to_numpy()
    geno = a.obs.genotype.astype(str).to_numpy()
    rep = (a.obs.genotype.astype(str) + "_" + a.obs.replicate.astype(str)).to_numpy()
    print(a.n_obs, "cells,", len(genes), "genes; replicates:", pd.Series(rep).value_counts().to_dict())

    # ---- programs
    prog = pd.DataFrame(index=genes)
    for lin, fate in LINEAGES.items():
        d = drivers("ctr", fate, genes)
        prog[f"{lin}_program"] = d.program
        for k in KERNELS:
            prog[f"{lin}_r_ctr_{k}"] = d[f"r_ctr_{k}"].round(4)
    ec_c, ec_k = drivers("ctr", "pre_EC", genes), drivers("KO", "pre_EC", genes)
    prog["EC_KO_program"] = ec_k.program
    ec = ec_c.drop(columns="program").join(ec_k.drop(columns="program")).round(4)
    ec["ctr_program"], ec["KO_program"] = ec_c.program, ec_k.program
    inc, ink = ec.ctr_program != "-", ec.KO_program != "-"
    ec["overlap"] = np.select([inc & ink, inc, ink], ["both", "control only", "KO only"], "-")
    ec = ec[ec.overlap != "-"].rename_axis("gene").reset_index()
    ec = ec.sort_values(["overlap", "r_ctr_pseudotime"], ascending=[True, False], kind="stable")
    ec.to_csv(TR / "ec_programs.csv", index=False)
    for lin in LINEAGES:
        print(f"{lin} control program: both kernels {int((prog[f'{lin}_program'] == 'both').sum())}, "
              f"one kernel {int((prog[f'{lin}_program'] == 'one').sum())}")
    print("EC program, control vs KO (any kernel):", ec.overlap.value_counts().to_dict())
    bb = ec[(ec.ctr_program == "both") & (ec.KO_program == "both")].gene.tolist()
    print(f"  in both genotypes' programs, both kernels ({len(bb)}):", ", ".join(bb[:40]))
    for o in ["control only", "KO only"]:
        s = ec[(ec.overlap == o) & ((ec.ctr_program == "both") | (ec.KO_program == "both"))]
        print(f"  {o}, both kernels ({len(s)}):", ", ".join(s.gene.head(40)))

    # ---- scores and fates per cell
    X = a.X
    score = {}
    for lin in LINEAGES:
        g = np.isin(genes, prog.index[prog[f"{lin}_program"] == "both"])
        score[f"{lin}_program_score"] = np.asarray(X[:, g].mean(1)).ravel()
    fates = {}
    for k in KERNELS:
        f = pd.read_csv(TR / f"cellrank_{k}_fates.csv", index_col=0, keep_default_na=False).reindex(a.obs_names)
        for fa in FATES:
            fates[f"fate_{fa}_{k}"] = f[fa].to_numpy()

    # ---- windows along each lineage
    sl = pd.read_csv(TR / "slingshot_pc30_cells.csv", index_col=0).reindex(a.obs_names)
    ln = pd.read_csv(TR / "slingshot_pc30_lineages.csv")
    end = {r.lineage: r.clusters.split(" -> ")[-1] for r in ln.itertuples()}
    Ls = {}
    for lin, fate in LINEAGES.items():
        L = [k for k, v in end.items() if v == fate]
        assert len(L) == 1, f"{fate}: lineages {L}"
        Ls[lin] = L[0]
    onl = {lin: (sl[f"w_{L}"].fillna(0).to_numpy() >= W_MIN) & sl[f"pt_{L}"].notna().to_numpy() for lin, L in Ls.items()}
    win_rows, tests, profile = [], [], []
    prof_g = np.isin(genes, prog.index[(prog[[f"{l}_program" for l in LINEAGES] + ["EC_KO_program"]] != "-").any(axis=1)])
    for lin, L in Ls.items():
        pt = sl[f"pt_{L}"].to_numpy(); on = onl[lin]
        other = np.logical_or.reduce([onl[o] for o in Ls if o != lin])
        edges = np.quantile(pt[on], np.linspace(0, 1, NWIN + 1))
        print(f"\n== {lin}: {L} ({ln.set_index('lineage').clusters[L]}); cells on it: KO {(on & (geno == 'KO')).sum()}, "
              f"ctr {(on & (geno == 'ctr')).sum()}")
        for i in range(NWIN):
            m = on & (pt >= edges[i]) & ((pt <= edges[i + 1]) if i == NWIN - 1 else (pt < edges[i + 1]))
            part = "shared path" if other[m].mean() >= 0.5 else "branch"
            mk, mc = m & (geno == "KO"), m & (geno == "ctr")
            tested = mk.sum() >= MIN_CELLS and mc.sum() >= MIN_CELLS
            labs = pd.Series(label[m]).value_counts()
            base = {"lineage": lin, "slingshot_lineage": L, "window": i + 1,
                    "pt_from": round(edges[i], 1), "pt_to": round(edges[i + 1], 1), "part": part,
                    "pct_KO": round(100 * mk.sum() / max(m.sum(), 1), 1),
                    "labels": "; ".join(f"{k} {v}" for k, v in labs.items())}
            rep_ok, n_rep = False, {}
            if tested:
                d, rep_ok, n_rep = test_window(a, mk, mc, rep)
                d = d.rename_axis("gene").reset_index()
                d.insert(1, "lineage", lin); d.insert(2, "window", i + 1); d.insert(3, "part", part)
                d.insert(4, "pt_range", f"{edges[i]:.1f}-{edges[i + 1]:.1f}")
                d.insert(5, "cells", f"KO {mk.sum()} / ctr {mc.sum()}")
                d = d.merge(prog[[c for c in prog.columns if c.endswith("program")]], left_on="gene", right_index=True)
                tests.append(d)
                c = d.call.value_counts()
                hk = d[d.call == "higher in KO"].sort_values("log2FC", ascending=False)
                lk = d[d.call == "lower in KO"].sort_values("log2FC")
                print(f"  w{i + 1} pt {edges[i]:.0f}-{edges[i + 1]:.0f} {part}: KO {mk.sum()} / ctr {mc.sum()}"
                      f"{'' if rep_ok else ' (replicate check not possible)'}; higher in KO {c.get('higher in KO', 0)}, "
                      f"lower {c.get('lower in KO', 0)}; labels: {base['labels']}")
                print("     higher:", ", ".join(hk.gene.head(15)), "| lower:", ", ".join(lk.gene.head(15)))
            else:
                print(f"  w{i + 1} pt {edges[i]:.0f}-{edges[i + 1]:.0f} {part}: KO {mk.sum()} / ctr {mc.sum()}: not tested")
            for g, mg in [("KO", mk), ("ctr", mc)]:
                r = dict(base, genotype=g, cells=int(mg.sum()), tested=tested, replicate_check=rep_ok)
                for r_ in sorted(set(rep)):
                    if r_.startswith(g + "_"):
                        r[f"cells_{r_}"] = int((mg & (rep == r_)).sum())
                for k_, v in {**score, **fates}.items():
                    r[k_] = round(float(np.mean(v[mg])), 4) if mg.any() else np.nan
                win_rows.append(r)
        # profile, equal-width bins
        e2 = np.linspace(np.nanmin(pt[on]), np.nanmax(pt[on]), NBIN + 1)
        b = np.clip(np.digitize(pt, e2) - 1, 0, NBIN - 1)
        for g in ["ctr", "KO"]:
            for i in range(NBIN):
                m = on & (geno == g) & (b == i)
                if m.any():
                    profile.append(pd.DataFrame({"gene": genes[prof_g], "lineage": lin, "genotype": g, "bin": i + 1,
                                                 "pseudotime": round((e2[i] + e2[i + 1]) / 2, 1), "cells": int(m.sum()),
                                                 "mean": np.asarray(X[m][:, prof_g].mean(0)).ravel().round(4)}))
    w = pd.DataFrame(win_rows)
    w.to_csv(TR / "ko_windows.csv", index=False)
    t = pd.concat(tests)
    t.to_csv(TR / "ko_vs_ctr_windows.csv.gz", index=False)
    pd.concat(profile).to_csv(TR / "ko_lineage_profile.csv", index=False)

    # ---- summaries
    cols = ["lineage", "window", "part", "pt_from", "pt_to", "genotype", "cells", "pct_KO",
            "EC_program_score", "beta_program_score", "alpha_program_score", "fate_pre_EC_pseudotime",
            "fate_pre_beta2_pseudotime", "fate_pre_alpha_pseudotime"]
    print("\nwindows (program scores = mean log1p of the both-kernel control program genes):\n",
          w[cols].round(3).to_string(index=False))
    s = t[t.call != "-"]
    for flag in ["EC_program", "beta_program", "alpha_program", "EC_KO_program"]:
        x = s[s[flag] == "both"].groupby(["lineage", "window", "call"]).size().unstack(fill_value=0)
        print(f"\ncalls among {flag} genes (both kernels), per lineage x window:\n", x.to_string())
    print("\nwrote", TR / "ko_windows.csv", ", ko_vs_ctr_windows.csv.gz, ec_programs.csv, ko_lineage_profile.csv")


if __name__ == "__main__":
    main()
