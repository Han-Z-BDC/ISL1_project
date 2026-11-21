"""Step 7.4h: D15 + D22 (alpha_beta_EC subset): KO vs control in matched pseudotime windows along the beta, alpha, EC2
and EC1 lineages (the 7.3g design on the 7.4 trajectory; after 7.4a, 7.4d, 7.4g).

Decisions (user, 2026-11-21): repeat the D15 trajectory analysis on D15 + D22. As 7.3g, with the 7.4g changes:
  lineages   the four 7.4a Slingshot lineages (scVI latent): beta (-> SC_beta, through D15 pre_EC, pre_beta1, pre_beta2),
             alpha (-> SC_alpha, through D15 pre_alpha_KO, pre_alpha), EC2 (-> SC_EC2, through D15 pre_EC), EC1
             (-> SC_EC1); a cell is on a lineage at curve weight >= W_MIN and with a label on the lineage's path
  windows    NWIN windows of equal cell numbers along each lineage; within a window, KO vs control is tested
             separately in each stage (D15, D22: separate sequencing batches) when that stage has >= MIN_CELLS cells
             of each genotype; replicate check needs >= MIN_REP cells in each of the 2 replicates per genotype.
             Part = "shared path" when >= half of the window's cells are also on another lineage, else "branch"
  test       as 7.3g (log2FC 0.5, Wilcoxon FDR < 0.05, both KO replicates beyond both control replicates)
  programs   control beta / alpha / EC1 / EC2 programs and KO EC1 / EC2 programs from the pseudotime-kernel CellRank
             drivers (7.4d, r >= 0.2, FDR < 0.05; velocity runs backwards on this data set)
  readouts   per window x stage x genotype: cells per replicate, labels, % KO, mean pseudotime-kernel fate toward
             SC_beta / SC_alpha / SC_EC1 / SC_EC2, program scores (mean log1p of each control program's genes)
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine_<subset>/trajectory/ (reruns overwrite):
  ko_windows.csv, ko_vs_ctr_windows.csv.gz, ec_programs.csv (control vs KO EC1 / EC2 programs)
Run: sbatch analysis/scripts/7.4h_ko_along_lineages.sb   (CELLRANK_ENV python; INT_SUBSET, default alpha_beta_EC)
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
R_MIN, FDR_MAX, LFC = 0.2, 0.05, 0.5
W_MIN, NWIN, MIN_CELLS, MIN_REP = 0.5, 8, 30, 5
LINEAGES = {"beta": "SC_beta", "alpha": "SC_alpha", "EC2": "SC_EC2", "EC1": "SC_EC1"}
FATES = ["SC_beta", "SC_alpha", "SC_EC1", "SC_EC2"]


def dense(m):
    return m.toarray() if sparse.issparse(m) else np.asarray(m)


def drivers(geno, fate, genes):
    d = pd.read_csv(TR / f"cellrank_pseudotime_{geno}_drivers.csv")
    x = d[(d.lineage == fate) & (d.cells == f"{geno} early")].set_index("gene").reindex(genes)
    assert x.r.notna().any(), f"no {geno} drivers toward {fate}"
    return x.r, x.FDR, np.where((x.r >= R_MIN) & (x.FDR < FDR_MAX), "yes", "-")


def pseudobulk(a, m):
    c = np.asarray(a.layers["soupx_counts"][m].sum(0)).ravel()
    return np.log2(c / max(c.sum(), 1) * 1e6 + 1)


def test(a, mk, mc, rep):
    Xk, Xc = dense(a.X[mk]), dense(a.X[mc])
    lfc = np.log2((np.expm1(Xk).mean(0) + 0.01) / (np.expm1(Xc).mean(0) + 0.01))
    p = stats.mannwhitneyu(Xk, Xc, axis=0, alternative="two-sided").pvalue
    q = stats.false_discovery_control(np.where(np.isnan(p), 1.0, p))
    d = pd.DataFrame({"mean_KO": Xk.mean(0).round(3), "mean_ctr": Xc.mean(0).round(3),
                      "pct_KO": (100 * (Xk > 0).mean(0)).round(1), "pct_ctr": (100 * (Xc > 0).mean(0)).round(1),
                      "log2FC": lfc.round(3), "FDR": q}, index=a.var_names)
    reps = sorted(set(rep[mk | mc]))
    n_rep = {r_: int(((mk | mc) & (rep == r_)).sum()) for r_ in reps}
    ok = sum(r_.startswith("KO_") for r_ in reps) >= 2 and sum(r_.startswith("ctr_") for r_ in reps) >= 2 \
        and min(n_rep.values()) >= MIN_REP
    pb = {r_: pseudobulk(a, (mk | mc) & (rep == r_)) for r_ in reps}
    for r_ in reps:
        d[f"pb_{r_}"] = pb[r_].round(3)
    if ok:
        ko_r = [v for k_, v in pb.items() if k_.startswith("KO_")]; ct_r = [v for k_, v in pb.items() if k_.startswith("ctr_")]
        lo, hi = np.minimum.reduce(ct_r) > np.maximum.reduce(ko_r), np.maximum.reduce(ct_r) < np.minimum.reduce(ko_r)
    else:
        lo = hi = np.zeros(len(d), bool)
    d["call"] = np.select([(lfc >= LFC) & (q < FDR_MAX) & hi, (lfc <= -LFC) & (q < FDR_MAX) & lo],
                          ["higher in KO", "lower in KO"], "-")
    return d, ok


def main():
    a = ad.read_h5ad(INT / "merged_counts.h5ad")
    keep = np.asarray((a.X > 0).sum(0)).ravel() >= 10
    a.layers["soupx_counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    a = a[:, keep].copy()
    sl = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0)
    a = a[sl.index].copy()
    genes = a.var_names
    label, stage = sl.label.to_numpy(), a.obs.stage.astype(str).to_numpy()
    geno = a.obs.genotype.astype(str).to_numpy()
    rep = (a.obs.genotype.astype(str) + "_" + a.obs["sample"].astype(str)).to_numpy()
    print(a.n_obs, "cells,", len(genes), "genes; replicates:", pd.Series(rep).value_counts().to_dict())

    prog = pd.DataFrame(index=genes)
    for lin, fate in LINEAGES.items():
        r, _, m = drivers("ctr", fate, genes)
        prog[f"{lin}_program"], prog[f"{lin}_r_ctr"] = m, r.round(4)
    ec_rows = []
    for lin in ["EC1", "EC2"]:
        rc, qc, mc_ = drivers("ctr", LINEAGES[lin], genes)
        rk, qk, mk_ = drivers("KO", LINEAGES[lin], genes)
        prog[f"{lin}_KO_program"] = mk_
        e = pd.DataFrame({"gene": genes, "fate": LINEAGES[lin], "r_ctr": rc.round(4).to_numpy(), "FDR_ctr": qc.to_numpy(),
                          "r_KO": rk.round(4).to_numpy(), "FDR_KO": qk.to_numpy(), "ctr_program": mc_, "KO_program": mk_})
        e["overlap"] = np.select([(e.ctr_program == "yes") & (e.KO_program == "yes"), e.ctr_program == "yes",
                                  e.KO_program == "yes"], ["both", "control only", "KO only"], "-")
        e = e[e.overlap != "-"]
        ec_rows.append(e)
        print(f"{lin} program, control vs KO:", e.overlap.value_counts().to_dict())
        for o in ["both", "control only", "KO only"]:
            s = e[e.overlap == o].sort_values("r_ctr" if o != "KO only" else "r_KO", ascending=False)
            print(f"   {o}:", ", ".join(s.gene.head(30)))
    pd.concat(ec_rows).to_csv(TR / "ec_programs.csv", index=False)
    for lin in LINEAGES:
        print(f"control {lin} program: {int((prog[f'{lin}_program'] == 'yes').sum())} genes")

    X = a.X
    score = {f"{lin}_program_score": np.asarray(X[:, (prog[f"{lin}_program"] == "yes").to_numpy()].mean(1)).ravel()
             for lin in LINEAGES}
    f = pd.read_csv(TR / "cellrank_pseudotime_fates.csv", index_col=0).reindex(a.obs_names)
    fates = {f"fate_{fa}": f[fa].to_numpy() for fa in FATES}

    ln = pd.read_csv(TR / "slingshot_scvi_lineages.csv")
    paths = {r.lineage: r.clusters.split(" -> ") for r in ln.itertuples()}
    Ls = {lin: [k for k, v in paths.items() if v[-1] == fate][0] for lin, fate in LINEAGES.items()}
    onl = {lin: ((sl[f"w_{L}"].fillna(0) >= W_MIN) & sl.label.isin(paths[L]) & sl[f"pt_{L}"].notna()).to_numpy()
           for lin, L in Ls.items()}
    win_rows, tests = [], []
    for lin, L in Ls.items():
        pt = sl[f"pt_{L}"].to_numpy(); on = onl[lin]
        other = np.logical_or.reduce([onl[o] for o in Ls if o != lin])
        edges = np.quantile(pt[on], np.linspace(0, 1, NWIN + 1))
        print(f"\n== {lin}: {L} ({' -> '.join(paths[L])}); on it: KO {(on & (geno == 'KO')).sum()}, "
              f"ctr {(on & (geno == 'ctr')).sum()}")
        for i in range(NWIN):
            m = on & (pt >= edges[i]) & ((pt <= edges[i + 1]) if i == NWIN - 1 else (pt < edges[i + 1]))
            part = "shared path" if other[m].mean() >= 0.5 else "branch"
            for st in ["D15", "D22"]:
                ms = m & (stage == st)
                if not ms.any():
                    continue
                mk, mc = ms & (geno == "KO"), ms & (geno == "ctr")
                tested = mk.sum() >= MIN_CELLS and mc.sum() >= MIN_CELLS
                labs = pd.Series(label[ms]).value_counts()
                base = {"lineage": lin, "slingshot_lineage": L, "window": i + 1, "stage": st,
                        "pt_from": round(edges[i], 2), "pt_to": round(edges[i + 1], 2), "part": part,
                        "pct_KO": round(100 * mk.sum() / max(ms.sum(), 1), 1),
                        "labels": "; ".join(f"{k} {v}" for k, v in labs.items())}
                rep_ok = False
                if tested:
                    d, rep_ok = test(a, mk, mc, rep)
                    d = d.rename_axis("gene").reset_index()
                    for j, (c, v) in enumerate([("lineage", lin), ("window", i + 1), ("stage", st), ("part", part),
                                                ("pt_range", f"{edges[i]:.2f}-{edges[i + 1]:.2f}"),
                                                ("cells", f"KO {mk.sum()} / ctr {mc.sum()}")]):
                        d.insert(1 + j, c, v)
                    d = d.merge(prog[[c for c in prog.columns if c.endswith("program")]], left_on="gene",
                                right_index=True)
                    tests.append(d)
                    c = d.call.value_counts()
                    hk = d[d.call == "higher in KO"].sort_values("log2FC", ascending=False).gene.head(12)
                    lk = d[d.call == "lower in KO"].sort_values("log2FC").gene.head(12)
                    print(f"  w{i + 1} {st} pt {edges[i]:.1f}-{edges[i + 1]:.1f} {part}: KO {mk.sum()} / ctr {mc.sum()}"
                          f"{'' if rep_ok else ' (no replicate check)'}; higher {c.get('higher in KO', 0)}, lower "
                          f"{c.get('lower in KO', 0)}; {base['labels'][:120]}\n     higher: {', '.join(hk)} | "
                          f"lower: {', '.join(lk)}")
                else:
                    print(f"  w{i + 1} {st} pt {edges[i]:.1f}-{edges[i + 1]:.1f} {part}: KO {mk.sum()} / ctr "
                          f"{mc.sum()}: not tested; {base['labels'][:120]}")
                for g, mg in [("KO", mk), ("ctr", mc)]:
                    r = dict(base, genotype=g, cells=int(mg.sum()), tested=tested, replicate_check=rep_ok)
                    for r_ in sorted(set(rep)):
                        if r_.startswith(g + "_"):
                            r[f"cells_{r_}"] = int((mg & (rep == r_)).sum())
                    for k_, v in {**score, **fates}.items():
                        r[k_] = round(float(np.nanmean(v[mg])), 4) if mg.any() else np.nan
                    win_rows.append(r)
    w = pd.DataFrame(win_rows); w.to_csv(TR / "ko_windows.csv", index=False)
    t = pd.concat(tests); t.to_csv(TR / "ko_vs_ctr_windows.csv.gz", index=False)
    cols = ["lineage", "window", "stage", "part", "genotype", "cells", "pct_KO"] + list(score) + list(fates)
    print("\nwindows:\n", w[cols].round(3).to_string(index=False))
    s = t[t.call != "-"]
    for flag in [c for c in t.columns if c.endswith("program")]:
        x = s[s[flag] == "yes"].groupby(["lineage", "window", "stage", "call"]).size().unstack(fill_value=0)
        print(f"\ncalls among {flag} genes:\n", x.to_string())
    n = s.groupby(["gene", "call"]).size().unstack(fill_value=0)
    ntest = t.groupby(["lineage", "window", "stage"]).ngroups
    lw = n[n.max(axis=1) >= 2 / 3 * ntest]
    print(f"\ngenes called in >= 2/3 of the {ntest} tests (KO-line-wide candidates): {len(lw)}:", ", ".join(lw.index[:60]))
    print("\nwrote", TR / "ko_windows.csv", ", ko_vs_ctr_windows.csv.gz, ec_programs.csv")


if __name__ == "__main__":
    main()
