"""Step 7.4f: alpha- vs beta-competent progenitors in the D15 + D22 trajectory (after 7.4c scVelo and 7.4d CellRank).

*** DRAFT (2026-11-20): plan under review by the user -- do NOT run. Thresholds (BIAS_HI, MIN_MASS), contrasts and
*** programs are proposals. Not written yet: 7.4f_pseudobulk_edger.R (edgeR QL with an optional sample block) and
*** 7.4f_competency.sb. See DECISIONS.md, "Competency question".

Question (user, 2026-11-20): experiments show the progenitors competent for alpha and for beta differ and their
competency changes over time: alpha competency is kept until ~D18 and then lost, beta competency until ~D26. In the
ISL1 KO, alpha-competent cells turn on NKX6-1 (absent from control alpha-competent cells) but keep alpha competency
until D18. Can the trajectory define the two progenitor kinds and say why they differ?
Design (user, 2026-11-20): run on the alpha_beta_EC subset trajectory (INT_SUBSET, default alpha_beta_EC); the D22
progenitors ("ENP", dropped from that subset) come in by expression only (no rerun of the integration / trajectory).

  groups      early D15 progenitors (ENP1, ENP2, ENP_FEV_high), each genotype. From the CellRank fate probabilities
              of both kernels (velocity, Slingshot pseudotime): alpha bias = P(alpha) / (P(alpha) + P(beta)), endocrine
              mass = P(alpha) + P(beta). alpha-biased: bias >= BIAS_HI and mass >= MIN_MASS in both kernels;
              beta-biased: bias <= 1 - BIAS_HI and mass >= MIN_MASS in both; EC-directed: mass < MIN_MASS in both;
              everything else: mixed. Groups come from fates only, never from markers (NKX6-1 does not mean beta
              competence in KO).
  checks      marker genes per group x genotype (% of cells > 0, mean log1p): NKX6-1 expected on in KO alpha-biased
              cells and off in control alpha-biased cells.
  where/when  per group: D15 state, cell-cycle phase, scVelo latent time, Slingshot time (weight-averaged).
  why         edgeR QL on per-sample pseudobulks (7.4f_pseudobulk_edger.R; >= MIN_CELLS cells per pseudobulk):
                alpha_vs_beta_ctr   alpha- vs beta-biased, control, paired by sample
                alpha_vs_beta_KO    the same in KO
                KO_vs_ctr_alpha     KO vs control within alpha-biased cells (NKX6-1 and ISL1-loss genes)
                KO_vs_ctr_beta      KO vs control within beta-biased cells
              programs: alpha program = up in alpha_vs_beta_ctr (FDR < 0.05, log2FC > 1, top N_PROG by log2FC);
              beta program = down in it. Core alpha genes = alpha program genes also up (log2FC > 0) in
              alpha_vs_beta_KO.
  time        D22 progenitors (label "ENP" of the full 7.2 integration, mostly control) by expression only:
              program scores (scanpy score_genes) and a logistic-regression classifier (alpha- vs beta-biased D15
              control cells, standardised log1p of the program genes; 5-fold CV accuracy reported) applied to D22 ENP.
              Per program gene: mean / % > 0 in D15 alpha-biased, D15 beta-biased and D22 ENP (control only). Stage and
              sequencing batch are the same thing here, so D15 vs D22 differences are descriptive.
Expression: SoupX counts normalised to 10,000 per cell, log1p (as 7.4d).
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine_<subset>/trajectory/competency/ (reruns overwrite):
  cells.csv, group_summary.csv, markers.csv, edger_<contrast>.csv, programs.csv, d22_enp.csv, program_genes_time.csv,
  groups_umap.png, fates_scatter.png, markers_dot.png, program_scores.png
Run: sbatch analysis/scripts/7.4f_competency.sb
"""
import os
import subprocess
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
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

OUT_DIR, REPO_DIR = Path(os.environ["OUT_DIR"]), Path(os.environ["REPO_DIR"])
SUBSET = os.environ.get("INT_SUBSET", "alpha_beta_EC")
DA = OUT_DIR / "day_analysis"
FULL = DA / "integration_D15_D22_endocrine"
INT = DA / f"integration_D15_D22_endocrine_{SUBSET}"
TR = INT / "trajectory"
OUT = TR / "competency"
EARLY = ["D15 ENP1", "D15 ENP2", "D15 ENP_FEV_high"]
KERNELS = ["velocity", "pseudotime"]
BIAS_HI, MIN_MASS = 2 / 3, 0.2
MIN_CELLS, N_PROG, SEED = 20, 50, 1234
GROUPS = ["alpha-biased", "beta-biased", "mixed", "EC-directed"]
GCOL = {"alpha-biased": "#d6604d", "beta-biased": "#2166ac", "mixed": "#b8b8b8", "EC-directed": "#1b9e77"}
MARKERS = ["NKX6-1", "PDX1", "ISL1", "ARX", "IRX1", "IRX2", "MAFB", "PAX4", "NKX2-2", "NEUROG3", "NEUROD1", "FEV",
           "LMX1A", "TPH1", "HHEX", "SST", "GCG", "INS", "CHGA", "MKI67"]


def lognorm(a):
    a = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a


def classify(f):
    out = pd.Series("mixed", index=f.index)
    hi, lo, ec = [], [], []
    for k in KERNELS:
        pa, pb = f[f"{k}_SC_alpha"], f[f"{k}_SC_beta"]
        m = pa + pb; b = pa / m.where(m > 0)
        f[f"{k}_mass"], f[f"{k}_alpha_bias"] = m, b
        hi.append((b >= BIAS_HI) & (m >= MIN_MASS)); lo.append((b <= 1 - BIAS_HI) & (m >= MIN_MASS))
        ec.append(m < MIN_MASS)
    out[np.logical_and(*hi)] = "alpha-biased"; out[np.logical_and(*lo)] = "beta-biased"
    out[np.logical_and(*ec)] = "EC-directed"
    return out


def pseudobulk(a, obs, name, ref, test, block=None):
    """edgeR on per-sample sums; obs needs 'cmp' (ref / test) and 'sample'."""
    d = OUT / f"pb_{name}"; d.mkdir(exist_ok=True)
    cols, meta = [], []
    for (s, g), idx in obs.groupby(["sample", "cmp"], observed=True).groups.items():
        if len(idx) < MIN_CELLS:
            continue
        cols.append(pd.Series(np.asarray(a[idx].X.sum(0)).ravel(), index=a.var_names, name=f"{s}|{g}"))
        meta.append({"pb": f"{s}|{g}", "group": g, "block": s if block else "", "cells": len(idx)})
    meta = pd.DataFrame(meta).set_index("pb")
    print(name, "pseudobulks:\n", meta.to_string())
    if (meta.group == ref).sum() < 2 or (meta.group == test).sum() < 2:
        print(name, "skipped: fewer than 2 pseudobulks per group"); return None
    pd.concat(cols, axis=1).to_csv(d / "pb_counts.csv"); meta.to_csv(d / "pb_meta.csv")
    r = subprocess.run([os.environ["TOOLS"] + "/Rscript", str(REPO_DIR / "analysis/scripts/7.4f_pseudobulk_edger.R"),
                        str(d), ref, test], capture_output=True, text=True)
    print(r.stdout, r.stderr[-2000:])
    if r.returncode:
        raise RuntimeError(f"edgeR failed for {name}")
    res = pd.read_csv(d / "edger.csv")
    res.to_csv(OUT / f"edger_{name}.csv", index=False)
    return res


def dotplot(expr, obs, genes, rows, fn, title):
    pct = pd.DataFrame({r: (expr.loc[m, genes] > 0).mean() * 100 for r, m in rows.items()}).T
    mean = pd.DataFrame({r: expr.loc[m, genes].mean() for r, m in rows.items()}).T
    fig, ax = plt.subplots(figsize=(0.55 * len(genes) + 3, 0.5 * len(rows) + 1.6))
    yy, xx = np.meshgrid(range(len(rows)), range(len(genes)), indexing="ij")
    s = ax.scatter(xx.ravel(), yy.ravel(), s=pct.to_numpy().ravel() * 3, c=mean.to_numpy().ravel(), cmap="Reds",
                   edgecolors="#999", linewidths=0.3)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{r} ({m.sum()})" for r, m in rows.items()])
    ax.invert_yaxis(); plt.colorbar(s, ax=ax, fraction=0.03, label="mean log1p"); ax.set_title(title, fontsize=10)
    fig.tight_layout(); fig.savefig(fn, dpi=90); plt.close(fig)
    return pct.round(1), mean.round(3)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in OUT.glob("*.*"):
        p.unlink()
    # ---- fates of both kernels, groups
    f = None
    for k in KERNELS:
        x = pd.read_csv(TR / f"cellrank_{k}_fates.csv", index_col=0)
        meta, fp = x[["label", "stage", "genotype"]], x.drop(columns=["label", "stage", "genotype"]).add_prefix(f"{k}_")
        f = pd.concat([meta, fp], axis=1) if f is None else f.join(fp)
    obs = pd.read_csv(INT / "cells.tsv", sep="\t", index_col=0)
    f = f.join(obs[["sample", "phase"]])
    sv = pd.read_csv(TR / "scvelo_cells.csv", index_col=0)
    f["latent_time"] = sv.latent_time.reindex(f.index)
    sl = pd.read_csv(TR / "slingshot_scvi_cells.csv", index_col=0).reindex(f.index)
    pt, w = np.nan_to_num(sl.filter(like="pt_").to_numpy()), np.nan_to_num(sl.filter(like="w_").to_numpy())
    f["slingshot_time"] = (pt * w).sum(1) / np.maximum(w.sum(1), 1e-12)
    early = f[f.label.isin(EARLY)].copy()
    for k in KERNELS:   # distributions, to judge the thresholds
        m = early[f"{k}_SC_alpha"] + early[f"{k}_SC_beta"]
        print(k, "endocrine mass quantiles by genotype:\n",
              m.groupby(early.genotype).quantile([0.1, 0.25, 0.5, 0.75, 0.9]).unstack().round(3).to_string())
        b = early[f"{k}_SC_alpha"] / m.where(m > 0)
        print(k, "alpha-bias quantiles by genotype:\n",
              b.groupby(early.genotype).quantile([0.1, 0.25, 0.5, 0.75, 0.9]).unstack().round(3).to_string())
    early["group"] = classify(early)
    early["group"] = pd.Categorical(early.group, categories=GROUPS)
    early.to_csv(OUT / "cells.csv")
    summ = []
    for (g, geno), d in early.groupby(["group", "genotype"], observed=True):
        r = {"group": g, "genotype": geno, "cells": len(d)}
        r.update({f"{l}_pct": round(100 * (d.label == l).mean(), 1) for l in EARLY})
        r.update({f"{p}_pct": round(100 * (d.phase == p).mean(), 1) for p in ["G1", "S", "G2M"]})
        r.update({"latent_time_median": d.latent_time.median(), "slingshot_time_median": d.slingshot_time.median()})
        r.update({f"{k}_{c}_mean": d[f"{k}_{c}"].mean() for k in KERNELS for c in ["SC_alpha", "SC_beta"]})
        summ.append(r)
    summ = pd.DataFrame(summ).round(3)
    summ.to_csv(OUT / "group_summary.csv", index=False)
    print("groups (early D15 progenitors):\n", summ.to_string())
    print(pd.crosstab([early.label, early.genotype], early.group).to_string())

    # ---- expression: subset counts + D22 ENP from the full integration
    a = ad.read_h5ad(INT / "merged_counts.h5ad")
    full = ad.read_h5ad(FULL / "merged_counts.h5ad", backed="r")
    enp = full.obs.index[(full.obs.label.astype(str) == "ENP").to_numpy()]
    e22 = full[enp].to_memory(); full.file.close()
    assert list(e22.var_names) == list(a.var_names)
    print("D22 ENP cells:", e22.n_obs, e22.obs.genotype.value_counts().to_dict())
    ae = a[early.index].copy()
    ln = lognorm(ae); le = lognorm(e22)
    expr = pd.DataFrame(ln[:, [g for g in MARKERS if g in ln.var_names]].X.toarray(), index=ln.obs_names,
                        columns=[g for g in MARKERS if g in ln.var_names])
    rows = {f"{g} {geno}": ((early.group == g) & (early.genotype == geno)).to_numpy()
            for g in GROUPS for geno in ["ctr", "KO"]}
    rows = {k: v for k, v in rows.items() if v.sum() >= 5}
    pct, mean = dotplot(expr, early, list(expr.columns), rows, OUT / "markers_dot.png",
                        "Early D15 progenitors: markers by fate group x genotype (size = % > 0, colour = mean log1p)")
    pd.concat({"pct_expressing": pct, "mean_log1p": mean}, axis=1).to_csv(OUT / "markers.csv")
    print("NKX6-1 % > 0 by group x genotype:\n", pct["NKX6-1"].to_string() if "NKX6-1" in pct else "NKX6-1 absent")

    # ---- figures: groups on the UMAP, fates
    ce = pd.read_csv(TR / "cells.csv", index_col=0)
    u = ce[["UMAP1", "UMAP2"]].to_numpy(); pos = pd.Series(range(len(ce)), index=ce.index)
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.5))
    for ax, geno in zip(axs, ["ctr", "KO"]):
        ax.scatter(u[:, 0], u[:, 1], c="#ececec", s=2, linewidths=0, rasterized=True)
        for g in GROUPS:
            i = pos[early.index[(early.group == g) & (early.genotype == geno)]].to_numpy()
            ax.scatter(u[i, 0], u[i, 1], c=GCOL[g], s=5, linewidths=0, rasterized=True, label=f"{g} ({len(i)})")
        ax.set_title(f"early D15 progenitors, {geno}"); ax.legend(markerscale=3, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(OUT / "groups_umap.png", dpi=85); plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(11, 10))
    for i, k in enumerate(KERNELS):
        for j, geno in enumerate(["ctr", "KO"]):
            d = early[early.genotype == geno]; ax = axs[i, j]
            ax.scatter(d[f"{k}_SC_beta"], d[f"{k}_SC_alpha"], c=d.group.map(GCOL), s=4, linewidths=0)
            ax.plot([0, 1], [1, 0], c="#ccc", lw=0.8)
            ax.set_xlabel("P(SC_beta)"); ax.set_ylabel("P(SC_alpha)"); ax.set_title(f"{k} kernel, {geno}")
    fig.tight_layout(); fig.savefig(OUT / "fates_scatter.png", dpi=85); plt.close(fig)

    # ---- pseudobulk edgeR
    res = {}
    for name, geno in [("alpha_vs_beta_ctr", "ctr"), ("alpha_vs_beta_KO", "KO")]:
        o = early[(early.genotype == geno) & early.group.isin(["alpha-biased", "beta-biased"])].copy()
        o["cmp"] = o.group.map({"alpha-biased": "alpha", "beta-biased": "beta"})
        res[name] = pseudobulk(ae, o, name, ref="beta", test="alpha", block=True)
    for name, g in [("KO_vs_ctr_alpha", "alpha-biased"), ("KO_vs_ctr_beta", "beta-biased")]:
        o = early[early.group == g].copy(); o["cmp"] = o.genotype
        res[name] = pseudobulk(ae, o, name, ref="ctr", test="KO")

    # ---- programs
    r = res["alpha_vs_beta_ctr"]
    if r is None:
        print("no control alpha vs beta test: programs, scores and the D22 step skipped"); return
    sig = r[r.FDR < 0.05]
    progs = {"alpha": sig[sig.logFC > 1].nlargest(N_PROG, "logFC").gene.tolist(),
             "beta": sig[sig.logFC < -1].nsmallest(N_PROG, "logFC").gene.tolist()}
    pr = pd.DataFrame([(p, g) for p, gs in progs.items() for g in gs], columns=["program", "gene"])
    pr = pr.merge(r[["gene", "logFC", "FDR"]], on="gene")
    for nm in ["alpha_vs_beta_KO", "KO_vs_ctr_alpha", "KO_vs_ctr_beta"]:
        if res[nm] is not None:
            pr = pr.merge(res[nm][["gene", "logFC", "FDR"]].rename(columns={"logFC": f"logFC_{nm}", "FDR": f"FDR_{nm}"}),
                          on="gene", how="left")
    if "logFC_alpha_vs_beta_KO" in pr:
        pr["core"] = np.where(pr.program == "alpha", pr.logFC_alpha_vs_beta_KO > 0, pr.logFC_alpha_vs_beta_KO < 0)
    pr.to_csv(OUT / "programs.csv", index=False)
    print("programs:", {p: len(g) for p, g in progs.items()})
    print(pr.head(30).to_string())

    # ---- D22 ENP by expression
    ctr = (early.genotype == "ctr").to_numpy()
    tr_m = ctr & early.group.isin(["alpha-biased", "beta-biased"]).to_numpy()
    comb = ad.concat({"D15": ln, "D22": le}, label="src", join="inner")
    comb.obs["grp"] = np.r_[early.group.astype(str) + " " + early.genotype, ["D22 ENP " + g for g in le.obs.genotype]]
    for p, gs in progs.items():
        if len(gs) >= 5:
            sc.tl.score_genes(comb, gs, score_name=f"{p}_score", random_state=SEED)
    genes = [g for g in progs["alpha"] + progs["beta"] if g in comb.var_names]
    out = {}
    if len(genes) >= 5 and tr_m.sum() >= 20:
        X = comb[:, genes].X; X = X.toarray() if sparse.issparse(X) else np.asarray(X)
        y = (early.group.to_numpy()[tr_m] == "alpha-biased").astype(int)
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=0.1))
        cv = cross_val_score(clf, X[:len(early)][tr_m], y, cv=5)
        clf.fit(X[:len(early)][tr_m], y)
        comb.obs["p_alpha_classifier"] = clf.predict_proba(X)[:, 1]
        out["cv_accuracy"] = cv.mean()
        print(f"classifier: {len(genes)} genes, 5-fold CV accuracy {cv.mean():.3f} (+/- {cv.std():.3f})")
    sc_cols = [c for c in comb.obs if c.endswith("_score") or c == "p_alpha_classifier"]
    tab = comb.obs.groupby("grp")[sc_cols].agg(["median", "mean", "size"]).round(3)
    tab.to_csv(OUT / "d22_enp.csv")
    print("program scores / classifier by group:\n", tab.to_string())
    fig, axs = plt.subplots(1, len(sc_cols), figsize=(5.5 * len(sc_cols), 5), squeeze=False)
    order = [g for g in [f"{x} ctr" for x in GROUPS] + [f"{x} KO" for x in GROUPS] + ["D22 ENP ctr", "D22 ENP KO"]
             if (comb.obs.grp == g).sum() >= 5]
    for ax, c in zip(axs[0], sc_cols):
        ax.boxplot([comb.obs.loc[comb.obs.grp == g, c] for g in order], showfliers=False)
        ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=60, ha="right"); ax.set_title(c)
    fig.tight_layout(); fig.savefig(OUT / "program_scores.png", dpi=85); plt.close(fig)
    # program genes over time (control only)
    m15a = comb.obs.grp.eq("alpha-biased ctr").to_numpy(); m15b = comb.obs.grp.eq("beta-biased ctr").to_numpy()
    m22 = comb.obs.grp.eq("D22 ENP ctr").to_numpy()
    Xg = comb[:, [g for g in pr.gene if g in comb.var_names]].X
    Xg = pd.DataFrame(Xg.toarray() if sparse.issparse(Xg) else Xg, columns=[g for g in pr.gene if g in comb.var_names])
    t = pd.DataFrame({"D15_alpha_mean": Xg[m15a].mean(), "D15_beta_mean": Xg[m15b].mean(), "D22_ENP_mean": Xg[m22].mean(),
                      "D15_alpha_pct": (Xg[m15a] > 0).mean() * 100, "D15_beta_pct": (Xg[m15b] > 0).mean() * 100,
                      "D22_ENP_pct": (Xg[m22] > 0).mean() * 100}).round(3)
    t = pr[["program", "gene", "logFC", "FDR"]].merge(t, left_on="gene", right_index=True)
    t.to_csv(OUT / "program_genes_time.csv", index=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
