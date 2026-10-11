#!/usr/bin/env python
"""
RNA QC step 1.1: count/gene pre-threshold (MAD or manual) + MT%/Ribo% cutoff.

Design (agreed lab default -- see configs/qc_RNA.yaml for the full rationale):
  - Step 1 default is 5-MAD outlier detection on log1p(total_counts) and
    log1p(n_genes_by_counts), matching sc-best-practices / IGVF-UCSD /
    CellCommander. A per-sample "manual" override with fixed log1p bounds is
    supported but must be explicitly set, never silently assumed.
  - Step 2's MT% cutoff is, by default, ALSO a one-sided (upper-only)
    MAD-based adaptive threshold now -- median(pct_counts_mt) + n_mads *
    MAD(pct_counts_mt), clamped to [mt_pct_floor, mt_pct_cap] -- not a
    single fixed percentage for every sample. Replaced a fixed 3.0% after
    checking real per-cell distributions (see
    scripts/_check_mt_distributions.py): several samples' MEDIAN MT% was
    already at or above 3% (e.g. HZ015_D0_rep1's median was 8.15%) --
    meaning the fixed cutoff was slicing through the bulk of the
    distribution, not a dying-cell tail -- while other samples (e.g.
    XW002's, median ~0.1%) were untouched by it. One fixed number can't be
    right for both. A "fixed" per-sample override is still supported (see
    mt_pct_strategy_by_sample), same pattern as step 1's mad/manual choice.
    pct_counts_ribo stays a plain fixed cutoff (default OFF -- diagnostic
    only) -- no equivalent adaptive concern was found for it.
  - Steps 1 and 2 are order-independent (both are just boolean masks
    combined with AND), so they're computed together here.
  - Every run writes a before/after distribution plot and a QC summary CSV
    -- never filters silently.

Takes CellRanger's FILTERED matrix (already-called cells), not the raw
one -- matching the lab's own notebook exactly. The MAD threshold assumes
a population that's mostly real cells with rare outliers; computed
against the raw matrix instead, the median/MAD would be dominated by
empty droplets, the same failure mode ATAC QC's own design deliberately
avoids for the same reason. The raw matrix is still used, just not here
-- step 1.2 (SoupX) reads it separately for its ambient-RNA background.

This script does ONE sample per invocation, meant to be called by
1.1_rna_qc.sb as a Slurm array task (one array index = one sample), the same
pattern as cutnrun-pipeline's 01_align.sb / submit_align.sh.

Usage:
  python 1.1_rna_qc.py --sample-id D22_rep1 \
      --h5-path /path/to/D22_rep1/outs/filtered_feature_bc_matrix.h5 \
      --config configs/qc_RNA.yaml \
      --outdir /path/to/sc_pipeline_out/rna_qc
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_config(config_path: str, sample_id: str) -> dict:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    strategy = cfg.get("filtering_strategy_by_sample", {}).get(
        sample_id, cfg["filtering_strategy"]
    )
    mt_pct_strategy = cfg.get("mt_pct_strategy_by_sample", {}).get(
        sample_id, cfg.get("mt_pct_strategy", "mad")
    )
    resolved = {
        "filtering_strategy": strategy,
        "n_mads": cfg.get("n_mads", 5),
        "mt_pct_strategy": mt_pct_strategy,
        "mt_pct_n_mads": cfg.get("mt_pct_n_mads", 5),
        "mt_pct_floor": cfg.get("mt_pct_floor"),
        "mt_pct_cap": cfg.get("mt_pct_cap"),
        "mt_pct_max": cfg.get("mt_pct_max_by_sample", {}).get(
            sample_id, cfg.get("mt_pct_max")
        ),
        "ribo_pct_max": cfg.get("ribo_pct_max"),
    }

    if strategy == "manual":
        bounds = cfg.get("manual_bounds_by_sample", {}).get(sample_id)
        if bounds is None:
            sys.exit(
                f"ERROR: sample '{sample_id}' is set to filtering_strategy "
                f"'manual' but has no entry under manual_bounds_by_sample in "
                f"{config_path}. Add explicit log1p-space bounds for it -- "
                f"manual overrides are never inferred."
            )
        resolved["manual_bounds"] = bounds

    if mt_pct_strategy == "fixed" and resolved["mt_pct_max"] is None:
        sys.exit(
            f"ERROR: sample '{sample_id}' is set to mt_pct_strategy 'fixed' "
            f"but has no mt_pct_max (or an entry under mt_pct_max_by_sample) "
            f"in {config_path}. Fixed overrides are never inferred."
        )

    return resolved


def is_mad_outlier(values: np.ndarray, n_mads: float) -> np.ndarray:
    """True where `values` is more than n_mads MADs below the median.
    One-sided on purpose: step 1 is a *low-quality* filter (empty droplets,
    degraded cells), not a doublet filter -- unusually high counts/genes are
    a doublet-detection problem (step 1.3), not this step's job."""
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    if mad == 0:
        return np.zeros(values.shape, dtype=bool)
    return values < (median - n_mads * mad)


def compute_mt_threshold(mt_values: np.ndarray, cfg: dict):
    """Returns the pct_counts_mt threshold to apply (cells ABOVE it fail),
    or None if MT filtering is off entirely.

    "mad" (default): one-sided (upper only) adaptive threshold -- low %MT
    is never a quality problem, so there is no lower bound, matching
    is_mad_outlier()'s own one-sidedness for step 1. median + n_mads*MAD,
    clamped to [mt_pct_floor, mt_pct_cap]:
      - mt_pct_floor stops the threshold from going BELOW a sane minimum
        for a sample whose own MT% spread is very tight (e.g. XW002,
        MAD ~0.03) -- without it, even a "generous" n_mads multiple is a
        tiny absolute number, stricter than any fixed cutoff would be,
        which would newly flag normal cells in an already-clean sample as
        outliers. Same purpose as scater::isOutlier's own min_diff arg.
      - mt_pct_cap stops a genuinely elevated-or-skewed baseline from
        producing an implausibly permissive threshold.
    "fixed": use mt_pct_max directly, the old single-number-for-everyone
    behavior -- only for a deliberate, documented per-sample override."""
    if cfg["mt_pct_strategy"] == "fixed":
        return cfg["mt_pct_max"]  # may itself be None -> no filter
    median = np.median(mt_values)
    mad = np.median(np.abs(mt_values - median))
    thresh = median + cfg["mt_pct_n_mads"] * mad
    if cfg["mt_pct_floor"] is not None:
        thresh = max(thresh, cfg["mt_pct_floor"])
    if cfg["mt_pct_cap"] is not None:
        thresh = min(thresh, cfg["mt_pct_cap"])
    return thresh


def compute_pass_mask(adata, cfg: dict) -> tuple[np.ndarray, dict]:
    """Returns (boolean mask of cells to KEEP, dict of per-criterion removed
    counts for logging)."""
    n = adata.n_obs
    removed = {}

    # --- step 1 ---
    if cfg["filtering_strategy"] == "mad":
        low_counts = is_mad_outlier(
            adata.obs["log1p_total_counts"].to_numpy(), cfg["n_mads"]
        )
        low_genes = is_mad_outlier(
            adata.obs["log1p_n_genes_by_counts"].to_numpy(), cfg["n_mads"]
        )
        step1_fail = low_counts | low_genes
        removed["step1_low_counts_mad"] = int(low_counts.sum())
        removed["step1_low_genes_mad"] = int(low_genes.sum())
    else:  # manual
        b = cfg["manual_bounds"]
        lc = adata.obs["log1p_total_counts"].to_numpy()
        lg = adata.obs["log1p_n_genes_by_counts"].to_numpy()
        low_counts = (lc < b["log1p_total_counts_min"]) | (
            lc > b["log1p_total_counts_max"]
        )
        low_genes = (lg < b["log1p_n_genes_by_counts_min"]) | (
            lg > b["log1p_n_genes_by_counts_max"]
        )
        step1_fail = low_counts | low_genes
        removed["step1_out_of_manual_bounds_counts"] = int(low_counts.sum())
        removed["step1_out_of_manual_bounds_genes"] = int(low_genes.sum())

    # --- step 2 ---
    mt_fail = np.zeros(n, dtype=bool)
    mt_values = adata.obs["pct_counts_mt"].to_numpy()
    mt_threshold = compute_mt_threshold(mt_values, cfg)
    if mt_threshold is not None:
        mt_fail = mt_values > mt_threshold
        removed["step2_high_mt"] = int(mt_fail.sum())
    else:
        removed["step2_high_mt"] = "not filtered (mt filtering off)"
    removed["mt_pct_threshold_used"] = round(float(mt_threshold), 3) if mt_threshold is not None else None

    ribo_fail = np.zeros(n, dtype=bool)
    if cfg["ribo_pct_max"] is not None:
        ribo_fail = adata.obs["pct_counts_ribo"].to_numpy() > cfg["ribo_pct_max"]
        removed["step2_high_ribo"] = int(ribo_fail.sum())
    else:
        removed["step2_high_ribo"] = "not filtered (ribo_pct_max: null, diagnostic only)"

    keep = ~(step1_fail | mt_fail | ribo_fail)
    return keep, removed


def plot_before_after(adata, keep: np.ndarray, sample_id: str, out_png: Path):
    metrics = [
        ("log1p_total_counts", "log1p(total counts)"),
        ("log1p_n_genes_by_counts", "log1p(n genes)"),
        ("pct_counts_mt", "% MT counts"),
        ("pct_counts_ribo", "% Ribo counts"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(18, 4))
    for ax, (col, label) in zip(axes, metrics):
        vals = adata.obs[col].to_numpy()
        # Shared bin edges, computed once from the full (wider-range) "before"
        # data and reused for both histograms. Without this, hist(bins=60)
        # picks DIFFERENT bin edges per call, sized to each array's own
        # range -- "after" (outliers removed) gets a narrower range and so
        # narrower bins, which collect fewer cells per bin even where the
        # underlying density barely changed. That made the bulk of the
        # distribution look like it lost far more cells than the ~2% actually
        # removed, when the real removal was concentrated at the tails.
        bin_edges = np.linspace(vals.min(), vals.max(), 61)
        ax.hist(vals, bins=bin_edges, color="#999999", alpha=0.6, label=f"before (n={len(vals)})")
        ax.hist(
            vals[keep], bins=bin_edges, color="#2b6cb0", alpha=0.7,
            label=f"after (n={int(keep.sum())})",
        )
        ax.set_xlabel(label)
        ax.set_ylabel("n cells")
        ax.legend(fontsize=7)
    fig.suptitle(f"{sample_id}: RNA QC step 1.1, before vs. after")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--h5-path", required=True, help="Path to CellRanger's filtered_feature_bc_matrix.h5 (called cells) -- NOT the raw matrix, which would dominate the MAD threshold with empty droplets")
    ap.add_argument("--config", required=True, help="Path to configs/qc_RNA.yaml")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "summary").mkdir(parents=True, exist_ok=True)

    cfg = load_config(args.config, args.sample_id)

    adata = sc.read_10x_h5(args.h5_path)
    adata.var_names_make_unique()
    # If this h5 is a multiome (RNA+ATAC) feature matrix, keep only Gene
    # Expression features -- this script is RNA QC only.
    if "feature_types" in adata.var.columns:
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True
    )

    n_before = adata.n_obs
    keep, removed = compute_pass_mask(adata, cfg)
    n_after = int(keep.sum())

    adata.obs["passed_qc_step1_2"] = keep

    plot_before_after(
        adata, keep, args.sample_id, outdir / "plots" / f"{args.sample_id}.qc_before_after.png"
    )

    # Full annotated object (all cells, QC metrics + pass/fail flag) -- kept
    # so a later look-back doesn't need to recompute metrics from scratch.
    adata.write_h5ad(outdir / "h5ad" / f"{args.sample_id}.qc_annotated.h5ad")
    # Filtered object -- what step 1.2 (SoupX) actually consumes.
    adata[keep].copy().write_h5ad(outdir / "h5ad" / f"{args.sample_id}.qc_filtered.h5ad")

    mt_pct_threshold_used = removed.pop("mt_pct_threshold_used", None)
    summary = {
        "sample_id": args.sample_id,
        "filtering_strategy": cfg["filtering_strategy"],
        "n_mads": cfg["n_mads"] if cfg["filtering_strategy"] == "mad" else "n/a",
        "mt_pct_strategy": cfg["mt_pct_strategy"],
        "mt_pct_n_mads": cfg["mt_pct_n_mads"] if cfg["mt_pct_strategy"] == "mad" else "n/a",
        "mt_pct_floor": cfg["mt_pct_floor"] if cfg["mt_pct_strategy"] == "mad" else "n/a",
        "mt_pct_cap": cfg["mt_pct_cap"] if cfg["mt_pct_strategy"] == "mad" else "n/a",
        "mt_pct_threshold_used": mt_pct_threshold_used,
        "ribo_pct_max": cfg["ribo_pct_max"],
        "n_cells_before": n_before,
        "n_cells_after": n_after,
        "n_cells_removed_total": n_before - n_after,
        **{f"n_removed_{k}": v for k, v in removed.items()},
    }
    pd.DataFrame([summary]).to_csv(
        outdir / "summary" / f"{args.sample_id}.qc_summary.csv", index=False
    )

    mt_note = f"mt_pct threshold used: {mt_pct_threshold_used}" if mt_pct_threshold_used is not None else "mt filtering off"
    print(f"[{args.sample_id}] {n_before} -> {n_after} cells "
          f"({n_before - n_after} removed, {mt_note}). Summary + plot written to {outdir}")


if __name__ == "__main__":
    main()
