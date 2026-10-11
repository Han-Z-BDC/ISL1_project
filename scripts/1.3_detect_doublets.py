#!/usr/bin/env python
"""
RNA QC step 1.3: doublet detection with scDblFinder + Scrublet.

Design (agreed lab default -- see configs/detect_doublets.yaml for the full
rationale):
  - scDblFinder runs first; its own predicted-doublet RATE on this sample
    becomes Scrublet's expected_doublet_rate (kept deliberately dependent,
    per-sample, rather than a fixed loading-chip prior).
  - A cell is a doublet if EITHER tool flags it (union) -- not only if both
    agree. Matches the lab's own notebook and IGVF-UCSD/CellCommander.
  - CellRanger's own doublet calls are NOT used as a third vote -- the
    lab's own validated two-tool union is used, not a 3-way consensus.
  - Doublets are removed from the working object (not just flagged).

Input is step 1.2's *_soupx_corrected.h5ad. One sample per invocation --
called by 1.3_detect_doublets.sb as a Slurm array task, same pattern as the
earlier steps.

Usage:
  python 1.3_detect_doublets.py --sample-id D22_rep1 \
      --soupx-h5ad /path/to/remove_background/h5ad/D22_rep1.soupx_corrected.h5ad \
      --config configs/detect_doublets.yaml \
      --outdir /path/to/sc_pipeline_out/detect_doublets
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

anndata2ri.activate()
anndata2ri.scipy2ri.activate()


def run_scdblfinder(adata, random_state: int) -> pd.DataFrame:
    scDblFinder = importr("scDblFinder")
    SingleCellExperiment = importr("SingleCellExperiment")

    data_mat = adata.X.T.astype("float32")
    ro.globalenv["data_mat"] = data_mat
    ro.r(f"set.seed({random_state})")

    sce_obj = scDblFinder.scDblFinder(
        SingleCellExperiment.SingleCellExperiment(
            assays=ro.ListVector({"counts": ro.globalenv["data_mat"]})
        )
    )
    col_data = anndata2ri.rpy2py(sce_obj.do_slot("colData"))
    return pd.DataFrame(
        {
            "scDblFinder_doublet_score": col_data["scDblFinder.score"].values,
            "scDblFinder_predicted_doublet": (col_data["scDblFinder.class"].values == "doublet"),
        },
        index=adata.obs_names,
    )


def run_scrublet(adata, expected_doublet_rate: float, cfg: dict, random_state: int) -> pd.DataFrame:
    adata_sc = adata.copy()
    sc.pp.scrublet(
        adata_sc,
        expected_doublet_rate=expected_doublet_rate,
        sim_doublet_ratio=cfg["sim_doublet_ratio"],
        n_neighbors=cfg["n_neighbors"],
        mean_center=cfg["mean_center"],
        normalize_variance=cfg["normalize_variance"],
        n_prin_comps=cfg["n_prin_comps"],
        random_state=random_state,
    )
    return pd.DataFrame(
        {
            "scrublet_doublet_score": adata_sc.obs["doublet_score"].values,
            "scrublet_predicted_doublet": adata_sc.obs["predicted_doublet"].values,
        },
        index=adata.obs_names,
    )


def plot_diagnostics(obs: pd.DataFrame, sample_id: str, out_png: Path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    axes[0].hist(obs["scDblFinder_doublet_score"], bins=100, color="#805ad5", alpha=0.8)
    axes[0].set_xlabel("scDblFinder score")
    axes[0].set_ylabel("n cells")
    axes[0].set_title("scDblFinder score distribution")

    counts = pd.Series(
        {
            "singlet\n(both)": (~obs["scDblFinder_predicted_doublet"] & ~obs["scrublet_predicted_doublet"]).sum(),
            "scDblFinder\nonly": (obs["scDblFinder_predicted_doublet"] & ~obs["scrublet_predicted_doublet"]).sum(),
            "Scrublet\nonly": (~obs["scDblFinder_predicted_doublet"] & obs["scrublet_predicted_doublet"]).sum(),
            "both\n(agree)": (obs["scDblFinder_predicted_doublet"] & obs["scrublet_predicted_doublet"]).sum(),
        }
    )
    axes[1].bar(counts.index, counts.values, color=["#999999", "#4299e1", "#ed8936", "#e53e3e"])
    axes[1].set_ylabel("n cells")
    axes[1].set_title("Doublet calls by tool (union = all but leftmost bar)")

    fig.suptitle(f"{sample_id}: doublet detection")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--soupx-h5ad", required=True, help="Output of 1.2_remove_background.py: {sample}.soupx_corrected.h5ad")
    ap.add_argument("--config", required=True, help="Path to configs/detect_doublets.yaml")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "summary").mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    adata = sc.read_h5ad(args.soupx_h5ad)
    n_before = adata.n_obs

    scdbl = run_scdblfinder(adata, cfg["random_state"])
    doublet_rate = float(scdbl["scDblFinder_predicted_doublet"].mean())

    scrub = run_scrublet(adata, doublet_rate, cfg["scrublet"], cfg["random_state"])

    obs = scdbl.join(scrub)
    obs["doublet_filter_union"] = obs[["scDblFinder_predicted_doublet", "scrublet_predicted_doublet"]].any(axis=1)
    obs["doublet_filter_consensus"] = obs[["scDblFinder_predicted_doublet", "scrublet_predicted_doublet"]].all(axis=1)

    for col in obs.columns:
        adata.obs[col] = obs[col].to_numpy()

    if cfg["save_diagnostic_plot"]:
        plot_diagnostics(obs, args.sample_id, outdir / "plots" / f"{args.sample_id}.doublet_diagnostics.png")

    # Full annotated object (all cells, both tools' scores + calls) --
    # kept so the union/consensus decision can be revisited without rerunning
    # scDblFinder/Scrublet from scratch.
    adata.write_h5ad(outdir / "h5ad" / f"{args.sample_id}.doublets_annotated.h5ad")

    adata_filtered = adata[~adata.obs["doublet_filter_union"]].copy()
    adata_filtered.write_h5ad(outdir / "h5ad" / f"{args.sample_id}.doublets_removed.h5ad")

    n_after = adata_filtered.n_obs
    summary = {
        "sample_id": args.sample_id,
        "n_cells_before": n_before,
        "n_cells_after": n_after,
        "n_doublets_removed_union": int(n_before - n_after),
        "n_scdblfinder_doublets": int(obs["scDblFinder_predicted_doublet"].sum()),
        "n_scrublet_doublets": int(obs["scrublet_predicted_doublet"].sum()),
        "n_consensus_doublets_both_agree": int(obs["doublet_filter_consensus"].sum()),
        "scdblfinder_doublet_rate": round(doublet_rate, 4),
        "scrublet_expected_doublet_rate_used": round(doublet_rate, 4),
    }
    pd.DataFrame([summary]).to_csv(outdir / "summary" / f"{args.sample_id}.doublets_summary.csv", index=False)

    print(f"[{args.sample_id}] {n_before} -> {n_after} cells "
          f"({n_before - n_after} doublets removed, union rule). "
          f"scDblFinder rate {doublet_rate:.2%} used as Scrublet's expected_doublet_rate.")


if __name__ == "__main__":
    main()
