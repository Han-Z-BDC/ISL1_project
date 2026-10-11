#!/usr/bin/env python
"""
Diagnostic (not a pipeline step): pull real pct_counts_mt distributions
(pre-filter) for the 15 renamed/used samples, to see whether the fixed
mt_pct_max=3.0 cutoff in 1.1_rna_qc.py is cutting into the bulk of the
distribution (too strict) or cleanly separating a small high-MT tail
(appropriately protective) -- and compare against what a MAD-based
adaptive threshold (matching the SAME n_mads mechanism already used for
low_counts/low_genes in this same script) would recover per sample.

Reads only .obs["pct_counts_mt"] from each sample's qc_annotated.h5ad
(pre-filter checkpoint, written by 1.1_rna_qc.py before its own cutoff is
applied) via backed='r' mode -- doesn't touch the expression matrix.

Usage:
  python check_mt_distributions.py <outdir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path(sys.argv[1])
outdir.mkdir(parents=True, exist_ok=True)

RNA_QC_H5AD_DIR = "/scratch/alpine/hzhu10@xsede.org/sc_pipeline_full/out/rna_qc/h5ad"

# old_name -> (new_name, diff_batch, stage, tech_batch) -- stage order matches
# the developmental time course, for consistent plotting/reporting order
SAMPLES = [
    ("HZ180_MM1016_MM1017_H1_hPSC_rep1", "HZ015_D0_rep1", "HZ015", "D0"),
    ("HZ184_MM1008_MM1009_H1_S1_rep1", "HZ015_D4_rep1", "HZ015", "D4"),
    ("HZ185_MM1010_MM1011_H1_S1_rep2", "HZ015_D4_rep2", "HZ015", "D4"),
    ("HZ189_MM1046_MM1053_H1_S2_rep2", "HZ015_D7_rep2", "HZ015", "D7"),
    ("HZ193_MM1050_MM1057_H1_S3_rep2", "HZ015_D9_rep2", "HZ015", "D9"),
    ("HZ_13_XW002_D7_rep1", "XW002_D7_rep1", "XW002", "D7"),
    ("HZ_14_XW002_D9_rep1", "XW002_D9_rep1", "XW002", "D9"),
    ("HZ_15_XW002_D12_rep1", "XW002_D12_rep1", "XW002", "D12"),
    ("HZ_16_XW002_D12_rep2", "XW002_D12_rep2", "XW002", "D12"),
    ("HZ_01_HZ021_D15_rep1", "HZ021_D15_rep1", "HZ021", "D15"),
    ("HZ_02_HZ021_D15_rep2", "HZ021_D15_rep2", "HZ021", "D15"),
    ("HZ_03_HZ021_D22_rep1", "HZ021_D22_rep1", "HZ021", "D22"),
    ("HZ_04_HZ021_D22_rep2", "HZ021_D22_rep2", "HZ021", "D22"),
    ("HZ_05_HZ019_D45_rep1", "HZ019_D45_rep1", "HZ019", "D45"),
    ("HZ_06_HZ019_D45_rep2", "HZ019_D45_rep2", "HZ019", "D45"),
]

MT_PCT_MAX_CURRENT = 3.0
N_MADS_OPTIONS = [3, 5]

records = []
summary_rows = []

for old_name, new_name, diff_batch, stage in SAMPLES:
    path = f"{RNA_QC_H5AD_DIR}/{old_name}.qc_annotated.h5ad"
    a = sc.read_h5ad(path, backed="r")
    mt = a.obs["pct_counts_mt"].to_numpy().astype(float)
    a.file.close()

    n = len(mt)
    med, mad = np.median(mt), np.median(np.abs(mt - np.median(mt)))
    pct_over_current = 100 * np.mean(mt > MT_PCT_MAX_CURRENT)

    row = {
        "old_name": old_name, "new_name": new_name, "diff_batch": diff_batch, "stage": stage,
        "n_cells": n,
        "median_mt_pct": round(med, 3), "mad_mt_pct": round(mad, 3),
        "p50": round(np.percentile(mt, 50), 2), "p75": round(np.percentile(mt, 75), 2),
        "p90": round(np.percentile(mt, 90), 2), "p95": round(np.percentile(mt, 95), 2),
        "p99": round(np.percentile(mt, 99), 2), "max": round(mt.max(), 2),
        "n_removed_at_3.0pct_fixed": int(np.sum(mt > MT_PCT_MAX_CURRENT)),
        "pct_removed_at_3.0pct_fixed": round(pct_over_current, 1),
    }
    for nmad in N_MADS_OPTIONS:
        thresh = med + nmad * mad
        row[f"adaptive_thresh_{nmad}mad"] = round(thresh, 2)
        row[f"n_removed_at_{nmad}mad"] = int(np.sum(mt > thresh))
        row[f"pct_removed_at_{nmad}mad"] = round(100 * np.mean(mt > thresh), 1)
    summary_rows.append(row)

    for v in mt:
        records.append({"new_name": new_name, "diff_batch": diff_batch, "stage": stage, "pct_counts_mt": v})
    print(f"[{new_name}] n={n} median={med:.2f} MAD={mad:.2f} "
          f"3.0%-fixed removes {row['n_removed_at_3.0pct_fixed']} ({row['pct_removed_at_3.0pct_fixed']}%) | "
          f"5-MAD adaptive ({row['adaptive_thresh_5mad']}%) removes {row['n_removed_at_5mad']} ({row['pct_removed_at_5mad']}%)",
          flush=True)

summary = pd.DataFrame(summary_rows)
summary.to_csv(outdir / "mt_pct_distribution_summary.csv", index=False)
print(f"\nWritten: {outdir / 'mt_pct_distribution_summary.csv'}")

long_df = pd.DataFrame(records)
long_df.to_parquet(outdir / "mt_pct_long.parquet")
print(f"Written: {outdir / 'mt_pct_long.parquet'}")

# violin plot, ordered by diff_batch then stage, with the fixed 3.0% line
order = [s[1] for s in SAMPLES]
fig, ax = plt.subplots(figsize=(16, 6))
data = [long_df.loc[long_df.new_name == name, "pct_counts_mt"].clip(upper=20) for name in order]
parts = ax.violinplot(data, showmedians=True, widths=0.8)
ax.axhline(MT_PCT_MAX_CURRENT, color="red", linestyle="--", linewidth=1.5, label=f"current fixed cutoff ({MT_PCT_MAX_CURRENT}%)")
ax.set_xticks(range(1, len(order) + 1))
ax.set_xticklabels(order, rotation=45, ha="right")
ax.set_ylabel("pct_counts_mt (clipped at 20 for display)")
ax.set_title("Pre-filter %MT distribution per sample (violin), current fixed 3.0% cutoff in red")
ax.legend()
plt.tight_layout()
plt.savefig(outdir / "mt_pct_violin.png", dpi=150)
print(f"Written: {outdir / 'mt_pct_violin.png'}")
