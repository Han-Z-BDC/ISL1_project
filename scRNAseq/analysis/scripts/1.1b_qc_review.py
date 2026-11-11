"""Step 1.1b: review step 1.1's QC metrics before committing to cutoffs.

Reads pipeline step 1.1's {sample}.qc_annotated.h5ad (all Cell Ranger-called
cells + QC metrics) and writes, to ${OUT_DIR}/rna_qc_review/:
  - qc_distributions.png      per-sample histograms (counts, genes, %MT, %ribo)
                              with the 3/4/5-MAD cutoffs drawn, plus counts-vs-genes
                              scatter coloured by %MT
  - qc_cutoff_scenarios.csv   cells removed per sample under several cutoff choices
  - qc_quantiles.csv          per-sample quantiles of each metric
  - qc_removed_markers.csv    marker-gene detection in cells kept vs removed by the
                              default (5-MAD) filter, to check whether the filter
                              hits a specific cell type rather than junk

No model fitting or dimensionality reduction (login-node safe: reads 4 small h5ad).
Run (from repo root): source study_config/config.sh && OUT_DIR=$OUT_DIR SAMPLE_SHEET=$SAMPLE_SHEET $TOOLS/python analysis/scripts/1.1b_qc_review.py
"""
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR = Path(os.environ["OUT_DIR"])
SAMPLE_SHEET = os.environ["SAMPLE_SHEET"]
IN = OUT_DIR / "rna_qc" / "h5ad"
RES = OUT_DIR / "rna_qc_review"
RES.mkdir(parents=True, exist_ok=True)

samples = pd.read_csv(SAMPLE_SHEET, sep="\t")["sample_id"].tolist()
MARKERS = ["INS", "GCG", "SST", "TPH1", "FEV", "CHGA", "NKX6-1", "SOX9", "KRT19",
           "COL1A1", "VIM", "MKI67", "TOP2A", "ALB", "CDX2"]


def mad_low(x, n):
    med = np.median(x)
    return med - n * np.median(np.abs(x - med))


def mt_thresh(x, n, floor=3.0, cap=15.0):
    med = np.median(x)
    return float(np.clip(med + n * np.median(np.abs(x - med)), floor, cap))


scen_rows, q_rows, mk_rows = [], [], []
fig, axes = plt.subplots(len(samples), 5, figsize=(24, 4.2 * len(samples)))

for i, s in enumerate(samples):
    a = ad.read_h5ad(IN / f"{s}.qc_annotated.h5ad")
    o = a.obs
    lc, lg = o["log1p_total_counts"].to_numpy(), o["log1p_n_genes_by_counts"].to_numpy()
    mt, ribo = o["pct_counts_mt"].to_numpy(), o["pct_counts_ribo"].to_numpy()
    n = a.n_obs

    for col in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]:
        qs = np.percentile(o[col], [1, 5, 25, 50, 75, 95, 99])
        q_rows.append({"sample": s, "metric": col, **{f"p{p}": round(v, 2) for p, v in zip([1, 5, 25, 50, 75, 95, 99], qs)}})

    # cutoff scenarios
    for nm in (3, 4, 5):
        c_lo, g_lo = np.expm1(mad_low(lc, nm)), np.expm1(mad_low(lg, nm))
        mt_t = mt_thresh(mt, nm)
        low = (lc < mad_low(lc, nm)) | (lg < mad_low(lg, nm))
        hmt = mt > mt_t
        scen_rows.append({"sample": s, "scenario": f"{nm}-MAD (counts/genes + %MT)", "min_counts": round(c_lo),
                          "min_genes": round(g_lo), "max_pct_mt": round(mt_t, 2), "n_cells": n,
                          "removed_low": int(low.sum()), "removed_mt": int(hmt.sum()),
                          "removed_total": int((low | hmt).sum()), "pct_removed": round(100 * (low | hmt).mean(), 1)})
    for mc, mg, mm in [(1000, 500, 10), (2000, 1000, 10), (3000, 1500, 5)]:
        low = (o["total_counts"] < mc) | (o["n_genes_by_counts"] < mg)
        hmt = o["pct_counts_mt"] > mm
        scen_rows.append({"sample": s, "scenario": f"fixed: counts>={mc}, genes>={mg}, %MT<={mm}", "min_counts": mc,
                          "min_genes": mg, "max_pct_mt": mm, "n_cells": n, "removed_low": int(low.sum()),
                          "removed_mt": int(hmt.sum()), "removed_total": int((low | hmt).sum()),
                          "pct_removed": round(100 * (low | hmt).mean(), 1)})

    # what does the default (5-MAD) filter remove?
    keep = o["passed_qc_step1_2"].to_numpy().astype(bool) if "passed_qc_step1_2" in o else ~((lc < mad_low(lc, 5)) | (lg < mad_low(lg, 5)) | (mt > mt_thresh(mt, 5)))
    X = a[:, MARKERS].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    for j, g in enumerate(MARKERS):
        mk_rows.append({"sample": s, "gene": g,
                        "pct_detected_kept": round(100 * (X[keep, j] > 0).mean(), 1),
                        "pct_detected_removed": round(100 * (X[~keep, j] > 0).mean(), 1) if (~keep).any() else np.nan})

    # plots
    ax = axes[i]
    for k, (vals, lab, logx) in enumerate([(o["total_counts"], "UMI counts / cell", True),
                                           (o["n_genes_by_counts"], "genes / cell", True),
                                           (mt, "% mitochondrial", False), (ribo, "% ribosomal", False)]):
        v = np.log10(vals) if logx else vals
        ax[k].hist(v, bins=120, color="#6b7a8f")
        ax[k].set_xlabel(("log10 " if logx else "") + lab)
        ax[k].set_title(f"{s}: {lab}", fontsize=10)
        for nm, colr in ((3, "#d9534f"), (4, "#f0ad4e"), (5, "#2f7d4f")):
            if k == 0:
                ax[k].axvline(np.log10(np.expm1(mad_low(lc, nm))), color=colr, ls="--", label=f"{nm} MAD")
            elif k == 1:
                ax[k].axvline(np.log10(np.expm1(mad_low(lg, nm))), color=colr, ls="--", label=f"{nm} MAD")
            elif k == 2:
                ax[k].axvline(mt_thresh(mt, nm), color=colr, ls="--", label=f"{nm} MAD")
        if k < 3:
            ax[k].legend(fontsize=8)
    order = np.argsort(mt)  # high %MT drawn on top
    sc = ax[4].scatter(np.log10(o["total_counts"].to_numpy()[order]), np.log10(o["n_genes_by_counts"].to_numpy()[order]),
                       c=np.clip(mt[order], 0, 15), s=3, cmap="viridis")
    ax[4].set_xlabel("log10 UMI counts")
    ax[4].set_ylabel("log10 genes")
    ax[4].set_title(f"{s}: counts vs genes (colour = %MT, capped 15)", fontsize=10)
    plt.colorbar(sc, ax=ax[4])

fig.tight_layout()
fig.savefig(RES / "qc_distributions.png", dpi=110)
pd.DataFrame(scen_rows).to_csv(RES / "qc_cutoff_scenarios.csv", index=False)
pd.DataFrame(q_rows).to_csv(RES / "qc_quantiles.csv", index=False)
pd.DataFrame(mk_rows).to_csv(RES / "qc_removed_markers.csv", index=False)
print("wrote", RES)
