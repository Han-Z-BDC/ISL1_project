#!/usr/bin/env python
"""
Coverage table for cross-sample RNA integration: diff_batch x stage x
tech_batch, parsed from each sample_id via split("_") -- matches the lab's
own notebook exactly (3.single.timepoint.joint.analysis.ipynb):

    diff_batch, stage, tech_batch = sample_id.split("_")

e.g. "HZ019_D45_rep1" -> diff_batch=HZ019, stage=D45, tech_batch=rep1.

This is deliberately run BEFORE any merge/integration step, not as a
formality -- which diff_batch x stage combinations exist (and which are
missing) directly determines how much a batch-correction method can be
trusted at each stage: a stage sampled from multiple diff_batches gives
Harmony/scVI real cross-batch signal to correct against; a stage sampled
from only one diff_batch gives them nothing to correct against, and
whatever they do there can't be verified from this data alone. Look at
the output before trusting correction on any sparse stage.

Usage:
  python 5.1_coverage_table.py --sample-ids-file samples.integration.txt --outdir out/
"""
import argparse
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_sample_id(sample_id: str):
    parts = sample_id.split("_")
    if len(parts) != 3:
        raise ValueError(
            f"Expected sample_id '{sample_id}' to split into exactly 3 "
            f"underscore-separated parts (diff_batch, stage, tech_batch), "
            f"matching the lab's own notebook convention -- got {len(parts)}: {parts}"
        )
    diff_batch, stage, tech_batch = parts
    return diff_batch, stage, tech_batch


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-ids-file", required=True, help="One sample_id per line, {diff_batch}_{stage}_{tech_batch} format")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sample_ids = [l.strip() for l in Path(args.sample_ids_file).read_text().splitlines() if l.strip()]
    rows = []
    for sid in sample_ids:
        diff_batch, stage, tech_batch = parse_sample_id(sid)
        rows.append(dict(sample_id=sid, diff_batch=diff_batch, stage=stage, tech_batch=tech_batch))
    meta = pd.DataFrame(rows)
    meta.to_csv(outdir / "sample_metadata.csv", index=False)

    # long-form coverage: one row per (diff_batch, stage), listing which tech_batches exist
    coverage = (
        meta.groupby(["diff_batch", "stage"])["tech_batch"]
        .apply(lambda s: ",".join(sorted(s)))
        .reset_index()
        .rename(columns={"tech_batch": "tech_batches_present"})
    )
    coverage["n_tech_batches"] = coverage["tech_batches_present"].str.count(",") + 1
    coverage.to_csv(outdir / "coverage_long.csv", index=False)

    # wide pivot: diff_batch (rows) x stage (cols), cell = n_tech_batches (0 = missing entirely)
    stage_order = sorted(meta["stage"].unique(), key=lambda s: int(s[1:]) if s[1:].isdigit() else s)
    pivot = coverage.pivot(index="diff_batch", columns="stage", values="n_tech_batches").reindex(columns=stage_order).fillna(0).astype(int)
    pivot.to_csv(outdir / "coverage_wide.csv")

    print("=== diff_batch x stage coverage (n_tech_batches per cell, 0 = missing) ===")
    print(pivot.to_string())
    print()
    n_possible = pivot.shape[0] * pivot.shape[1]
    n_present = (pivot > 0).sum().sum()
    print(f"{n_present}/{n_possible} diff_batch x stage combinations present "
          f"({100*n_present/n_possible:.0f}%) -- {n_possible - n_present} missing entirely.")
    single_tech = ((pivot == 1)).sum().sum()
    print(f"{single_tech} combinations have only 1 tech_batch (no internal replicate to cross-check against).")

    # single-diff_batch stages: the ones where no other diff_batch exists to correct against
    n_diffbatch_per_stage = (pivot > 0).sum(axis=0)
    lonely_stages = n_diffbatch_per_stage[n_diffbatch_per_stage <= 1].index.tolist()
    if lonely_stages:
        print(f"\nWARNING: stage(s) with only ONE diff_batch present -- Harmony/scVI have "
              f"nothing to correct diff_batch effect against here, trust these least: {lonely_stages}")

    # heatmap
    fig, ax = plt.subplots(figsize=(1.2 * len(stage_order) + 2, 0.6 * pivot.shape[0] + 2))
    im = ax.imshow(pivot.values, cmap="Greens", vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            ax.text(j, i, str(v) if v else "-", ha="center", va="center",
                     color="white" if v >= 2 else "black", fontsize=11)
    ax.set_xlabel("stage")
    ax.set_ylabel("diff_batch")
    ax.set_title("diff_batch x stage coverage (n_tech_batches; 0/'-' = missing)")
    fig.colorbar(im, ax=ax, label="n_tech_batches", ticks=[0, 1, 2])
    fig.tight_layout()
    fig.savefig(outdir / "coverage_heatmap.png", dpi=150)
    plt.close(fig)

    print(f"\nWritten: {outdir}/sample_metadata.csv, coverage_long.csv, coverage_wide.csv, coverage_heatmap.png")


if __name__ == "__main__":
    main()
