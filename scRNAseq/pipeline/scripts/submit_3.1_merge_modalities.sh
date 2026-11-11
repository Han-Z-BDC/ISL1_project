#!/bin/bash
# Submits 3.1_merge_modalities.sb as a Slurm array sized to however many
# samples are currently in samples.tsv -- run this instead of
# `sbatch scripts/3.1_merge_modalities.sb` directly. Needs both RNA QC
# (through 1.3_detect_doublets) and ATAC QC (through 2.2_atac_qc) already
# finished for every sample.
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "${P}/config.sh"

n=$(tail -n +2 "${SAMPLE_SHEET}" | grep -c .)
(( n > 0 )) || { echo "ERROR: no sample rows found in ${SAMPLE_SHEET}" >&2; exit 1; }

echo "Submitting 3.1_merge_modalities.sb as an array of ${n} task(s) (--array=0-$((n-1)))"
sbatch --array="0-$((n-1))" "${P}/scripts/3.1_merge_modalities.sb"
