"""Step 7.3i: alpha vs beta bias of the endocrine progenitors, from the existing CellRank fates (after 7.3d and 7.4d).

Request (user, 2026-11-21): before clustering the progenitors by alpha / beta potential, look at the per-cell bias:
bias on the 2D UMAP, its distribution (clear separation?), bias vs the EC fate and bias vs a single-value pseudotime.
Nothing is rerun: the 7.3d / 7.4d fate probabilities (runs on all cells) are read as they are.
  bias       P(alpha) / (P(alpha) + P(beta)) per cell and kernel (pseudotime = main, velocity = check: velocity runs
             backwards on these data); S = P(alpha) + P(beta); cells with S < 0.05 marked (ratio on small numbers)
  D15        7.3d fates; alpha = pre_alpha, beta = pre_beta2, EC = pre_EC (an end state there, so pre_EC cells have no
             bias and are left out); cells ENP1, ENP2, ENP_FEV_high, pre_alpha_KO (user: include pre_alpha_KO);
             2D UMAP of the D15 endocrine subset (7.1d); Slingshot PCs 1-30 (7.3a)
  D15_D22    7.4d fates (alpha_beta_EC subset); alpha = SC_alpha, beta = SC_beta, EC = SC_EC1 + SC_EC2 (D22 ends, so the
             D15 pre_EC cells get a bias); cells D15 ENP1, ENP2, ENP_FEV_high, pre_alpha_KO, pre_EC (user); scVI 2D UMAP
             of the subset (7.2); Slingshot on the scVI latent (7.4a)
  pseudotime single value = Slingshot pseudotimes averaged with the curve weights, scaled 0-1 (as the CellRank
             pseudotime kernel used)
Writes ${OUT_DIR}/day_analysis/D15/trajectory/alpha_bias/{D15,D15_D22}/ (reruns overwrite):
  alpha_bias_cells_{kernel}.csv   per cell: label, genotype, fates, P_alpha, P_beta, P_EC, S, bias, pseudotime, UMAP
  check1.csv                      per subcluster x genotype x kernel: cells, S = 0, S < 0.01 / 0.05, S quantiles, medians
  correlations.csv                Spearman of the bias with P(EC) and with pseudotime (all cells, per subcluster)
  umap_{kernel}.png, hist_{kernel}.png, vs_ec_{kernel}.png, vs_pt_{kernel}.png
Report: analysis/scripts/build_alpha_bias_report.py on the Mac.
Run: sbatch analysis/scripts/7.3i_alpha_bias.sb   (CELLRANK_ENV python)
"""
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DA = Path(os.environ["OUT_DIR"]) / "day_analysis"
OUT = DA / "D15" / "trajectory" / "alpha_bias"
GT = [("ctr", "control"), ("KO", "ISL1 KO")]
S_MIN = 0.05
BIAS_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("bias", ["#08519c", "#7b3fa0", "#cb181d"])
LCOL = {"ENP1": "#6a3d9a", "ENP2": "#1f78b4", "ENP_FEV_high": "#e6550d", "pre_alpha_KO": "#c51b7d", "pre_EC": "#1b7837"}
INT = DA / "integration_D15_D22_endocrine_alpha_beta_EC"
SETS = {
    "D15": dict(tr=DA / "D15" / "trajectory", alpha="pre_alpha", beta="pre_beta2", ec=["pre_EC"],
                labels=["ENP1", "ENP2", "ENP_FEV_high", "pre_alpha_KO"], prefix="",
                umap=(DA / "D15" / "subcluster_endocrine" / "subclusters.csv", ["UMAP1", "UMAP2"]),
                sling="slingshot_pc30_cells.csv", title="D15 only"),
    "D15_D22": dict(tr=INT / "trajectory", alpha="SC_alpha", beta="SC_beta", ec=["SC_EC1", "SC_EC2"],
                    labels=["ENP1", "ENP2", "ENP_FEV_high", "pre_alpha_KO", "pre_EC"], prefix="D15 ",
                    umap=(INT / "umap_scvi.csv", ["UMAP1", "UMAP2"]),
                    sling="slingshot_scvi_cells.csv", title="D15 + D22 (alpha_beta_EC subset)"),
}


def sp(a, b):
    return spearmanr(a, b, nan_policy="omit").correlation if len(a) > 2 else np.nan


for name, c in SETS.items():
    od = OUT / name; od.mkdir(parents=True, exist_ok=True)
    LAB = c["labels"]
    um = pd.read_csv(c["umap"][0], index_col=0)[c["umap"][1]]
    um.columns = ["UMAP1", "UMAP2"]
    sl = pd.read_csv(c["tr"] / c["sling"], index_col=0)
    pt = np.nan_to_num(sl.filter(like="pt_").to_numpy()); w = np.nan_to_num(sl.filter(like="w_").to_numpy())
    t = (pt * w).sum(1) / np.maximum(w.sum(1), 1e-12)
    sl["pseudotime"] = t / t.max()
    check, corr = [], []
    for k in ["pseudotime", "velocity"]:
        f = pd.read_csv(c["tr"] / f"cellrank_{k}_fates.csv", index_col=0, keep_default_na=False)
        f["lab"] = f.label.astype(str).str.replace(c["prefix"], "", regex=False)
        f["P_alpha"], f["P_beta"] = f[c["alpha"]].astype(float), f[c["beta"]].astype(float)
        f["P_EC"] = f[c["ec"]].astype(float).sum(1)
        f["S"] = f.P_alpha + f.P_beta
        f["bias"] = f.P_alpha / f.S.where(f.S > 0)
        f["pseudotime"] = sl.pseudotime.reindex(f.index)
        f = f.join(um)
        e = f[f.lab.isin(LAB) & f.label.str.startswith(c["prefix"])]
        e.to_csv(od / f"alpha_bias_cells_{k}.csv")
        for (lab, g), s in e.groupby(["lab", "genotype"]):
            check.append({"kernel": k, "subcluster": lab, "genotype": g, "cells": len(s),
                          "S = 0": int((s.S == 0).sum()), "S < 0.01": int((s.S < 0.01).sum()),
                          "S < 0.05": int((s.S < 0.05).sum()), "S p10": s.S.quantile(.1), "S median": s.S.median(),
                          "S p90": s.S.quantile(.9), "P_alpha median": s.P_alpha.median(),
                          "P_beta median": s.P_beta.median(), "P_EC median": s.P_EC.median(),
                          "bias median": s.bias.median()})
        lowS = (e.S < S_MIN).to_numpy()

        # UMAP
        fig, axs = plt.subplots(1, 2, figsize=(15, 6.6))
        for ax, (g, gn) in zip(axs, GT):
            ax.scatter(f.UMAP1, f.UMAP2, s=2, c="#e3e3e3", linewidths=0, rasterized=True)
            m = (e.genotype == g).to_numpy()
            s = e[m & ~lowS].sort_values("bias", key=lambda v: (v - 0.5).abs(), kind="stable")
            sc = ax.scatter(s.UMAP1, s.UMAP2, s=6, c=s.bias, cmap=BIAS_CMAP, vmin=0, vmax=1, linewidths=0, rasterized=True)
            z = e[m & lowS]
            ax.scatter(z.UMAP1, z.UMAP2, s=7, c="#333333", marker="x", linewidths=0.6, rasterized=True,
                       label=f"P(alpha)+P(beta) < {S_MIN} ({len(z)})")
            for lab in LAB:
                q = e[e.lab == lab][["UMAP1", "UMAP2"]].median()
                ax.text(q.UMAP1, q.UMAP2, lab, fontsize=9, weight="bold", ha="center")
            ax.set_title(f"{gn}: {', '.join(LAB)} ({m.sum():,} cells)", fontsize=10)
            ax.set_xticks([]); ax.set_yticks([]); ax.legend(loc="lower left", fontsize=8, frameon=False)
        cb = fig.colorbar(sc, ax=axs, fraction=0.025)
        cb.set_label("alpha bias = P(alpha) / (P(alpha) + P(beta))   0 = beta, 1 = alpha")
        fig.suptitle(f"{c['title']}, {k} kernel: alpha vs beta bias on the 2D UMAP (other cells grey)", fontsize=12)
        fig.savefig(od / f"umap_{k}.png", dpi=90, bbox_inches="tight"); plt.close(fig)

        # histograms
        fig, axs = plt.subplots(1, len(LAB), figsize=(4.3 * len(LAB), 4.4), sharex=True)
        bins = np.linspace(0, 1, 41)
        for ax, lab in zip(axs, LAB):
            for g, gn in GT:
                v = e[(e.lab == lab) & (e.genotype == g)].bias.dropna()
                if len(v) == 0:
                    continue
                ax.hist(v, bins=bins, density=True, histtype="stepfilled" if g == "ctr" else "step",
                        alpha=0.45 if g == "ctr" else 1, lw=2, color="#2166ac" if g == "ctr" else "#d6604d",
                        label=f"{gn} ({len(v)}; median {v.median():.2f})")
            ax.set_title(lab); ax.set_xlabel("alpha bias (0 = beta, 1 = alpha)"); ax.legend(fontsize=7.5, frameon=False)
        axs[0].set_ylabel("density")
        fig.suptitle(f"{c['title']}, {k} kernel: distribution of the alpha bias per subcluster", fontsize=12)
        fig.tight_layout(); fig.savefig(od / f"hist_{k}.png", dpi=90); plt.close(fig)

        # scatter plots
        for ycol, ylab, tag in [("P_EC", "P(EC) fate", "ec"), ("pseudotime", "pseudotime (curve-weighted, 0-1)", "pt")]:
            fig, axs = plt.subplots(1, 2, figsize=(14, 5.6), sharex=True, sharey=True)
            for ax, (g, gn) in zip(axs, GT):
                v = e[e.genotype == g]
                for lab in LAB:
                    s = v[v.lab == lab]
                    ax.scatter(s.bias, s[ycol], s=5, alpha=0.5, c=LCOL[lab], linewidths=0, label=lab, rasterized=True)
                r = {lab: sp(v[v.lab == lab].bias, v[v.lab == lab][ycol]) for lab in LAB}
                ra = sp(v.bias, v[ycol])
                ax.set_title(f"{gn}: Spearman all {ra:+.2f}\n(" + ", ".join(f"{l} {x:+.2f}" for l, x in r.items()) + ")",
                             fontsize=9)
                ax.set_xlabel("alpha bias (0 = beta, 1 = alpha)"); ax.legend(fontsize=8, frameon=False, markerscale=3)
                corr.append({"kernel": k, "y": ycol, "genotype": g, "all": round(ra, 3),
                             **{lab: round(x, 3) for lab, x in r.items()}})
            axs[0].set_ylabel(ylab)
            fig.suptitle(f"{c['title']}, {k} kernel: alpha bias vs {ylab}", fontsize=12)
            fig.tight_layout(); fig.savefig(od / f"vs_{tag}_{k}.png", dpi=90); plt.close(fig)
    ch = pd.DataFrame(check)
    ch["subcluster"] = pd.Categorical(ch.subcluster, LAB)
    ch = ch.sort_values(["kernel", "subcluster", "genotype"], ascending=[False, True, False], kind="stable").round(4)
    ch.to_csv(od / "check1.csv", index=False)
    pd.DataFrame(corr).to_csv(od / "correlations.csv", index=False)
    print(f"\n===== {name}\n", ch.to_string(index=False), "\n", pd.DataFrame(corr).to_string(index=False))
print("wrote", OUT)
