#!/usr/bin/env python
"""
ATAC QC, AMULET step -- run separately from everything else, on purpose.

AMULET is slow and its runtime doesn't depend on QC thresholds -- it scans
the ENTIRE raw fragment file genome-wide, for every barcode, regardless of
how tightly configs/atac_qc.yaml's general-QC thresholds are set (this is
inherent to what it computes, not something --qc-only in 2.2_atac_qc.py can
speed up). The lab's own production job for this (jobs/templates/
single_cell/batch_amulet.sb, calling software/python_scripts/
run_Amulet.py) requests 3 days -- confirming this is a genuinely long,
independent cost, not something to bury inside a single combined script.

This script runs the exact same AMULET call as that validated script and
as the lab's own notebook: amulet() from the R scDblFinder package, on the
RAW (unfiltered) fragment file -- never the QC-filtered one, since
AMULET's genome-wide coverage computation needs the full droplet
population, not just called cells -- excluding the lab's own blacklist BED
plus chrM/chrX/chrY, q-value < 0.05.

Unlike the lab's own script (a single sequential Python loop over every
sample, hence the 3-day multi-sample budget), this runs ONE sample per
invocation as a Slurm array task -- so per-task wall time only needs to
cover one sample's AMULET run (observed elsewhere to be on the order of a
few hours), not all of them end to end.

Writes {sample_id}.amulet_result.csv (barcode, p.value, q.value,
amulet_predicted_doublet) to <outdir>/amulet/ -- 2.2_atac_qc.py reads this
file rather than calling AMULET itself, so tuning QC thresholds with
qc_only_atac.sh never re-triggers it, and a failed/slow AMULET run never
takes the rest of ATAC QC down with it.

Usage:
  python 2.1_atac_amulet.py --sample-id D22_rep1 \
      --frag-path /path/to/D22_rep1/outs/atac_fragments.tsv.gz \
      --blacklist-bed /path/to/blacklist.bed \
      --config configs/atac_qc.yaml \
      --outdir /path/to/sc_pipeline_out/atac_qc
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

anndata2ri.activate()
anndata2ri.scipy2ri.activate()


def run_amulet(frag_path: str, blacklist_bed: str, random_state: int) -> pd.DataFrame:
    """Runs AMULET (bundled in the scDblFinder R package) on raw
    fragments. Returns a DataFrame indexed by RAW barcode (no sample
    prefix) with p.value / q.value columns."""
    importr("scDblFinder")
    GenomicRanges = importr("GenomicRanges")
    importr("rtracklayer")
    ro.r(f"set.seed({random_state})")
    ro.globalenv["frag_path"] = frag_path
    ro.globalenv["blacklist_path"] = blacklist_bed
    ro.r("""
        repeats <- rtracklayer::import(blacklist_path)
        otherChroms <- GenomicRanges::GRanges(c("chrM", "chrX", "chrY", "MT"), IRanges::IRanges(1L, width=10^8))
        toExclude <- suppressWarnings(c(repeats, otherChroms))
        amulet_result <- amulet(frag_path, regionsToExclude=toExclude)
    """)
    return anndata2ri.rpy2py(ro.globalenv["amulet_result"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--frag-path", required=True, help="Path to the RAW atac_fragments.tsv.gz -- never a QC-filtered file")
    ap.add_argument("--blacklist-bed", required=True)
    ap.add_argument("--config", required=True, help="Path to configs/atac_qc.yaml")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "amulet").mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    print(f"[{args.sample_id}] Running AMULET on the raw fragment file -- this can take a while.")
    amulet_result = run_amulet(args.frag_path, args.blacklist_bed, cfg["random_state"])

    q_threshold = cfg["amulet"]["q_value_threshold"]
    amulet_result["amulet_predicted_doublet"] = amulet_result["q.value"] < q_threshold

    out_csv = outdir / "amulet" / f"{args.sample_id}.amulet_result.csv"
    amulet_result.to_csv(out_csv, index_label="barcode")

    with plt.rc_context():
        neg_log10_q = (-1 * np.log10(amulet_result["q.value"])).clip(upper=10)
        _, ax = plt.subplots(figsize=(4, 4))
        ax.hist(neg_log10_q, bins=100, color="#ed8936")
        ax.set_xlabel("AMULET score (-log10 q-value, clipped at 10)")
        ax.set_ylabel("n barcodes")
        ax.set_title(f"{args.sample_id}: AMULET")
        plt.tight_layout()
        plt.savefig(outdir / "amulet" / f"{args.sample_id}.amulet_score.png", dpi=150)
        plt.close()

    n_doublets = int(amulet_result["amulet_predicted_doublet"].sum())
    print(f"[{args.sample_id}] AMULET done: {n_doublets} / {len(amulet_result)} barcodes "
          f"flagged (q < {q_threshold}). Written to {out_csv}")


if __name__ == "__main__":
    main()
