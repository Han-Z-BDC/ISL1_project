#!/bin/bash
# Submits 2.2_atac_qc_only.sb (--qc-only mode) as a Slurm array sized to
# however many samples are currently in samples.tsv. Run this BEFORE
# scripts/submit_2.2_atac_qc.sh -- generates every sample's
# threshold_review.png for review, without paying for the full pipeline.
# Independent of RNA QC and AMULET -- can run in parallel with both.
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "${P}/config.sh"

n=$(tail -n +2 "${SAMPLE_SHEET}" | grep -c .)
(( n > 0 )) || { echo "ERROR: no sample rows found in ${SAMPLE_SHEET}" >&2; exit 1; }

echo "Submitting 2.2_atac_qc_only.sb as an array of ${n} task(s) (--array=0-$((n-1)))"
sbatch --array="0-$((n-1))" "${P}/scripts/2.2_atac_qc_only.sb"
