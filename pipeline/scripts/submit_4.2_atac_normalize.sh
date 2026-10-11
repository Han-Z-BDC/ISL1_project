#!/bin/bash
# Submits 4.2_atac_normalize.sb as a Slurm array sized to however many
# samples are currently in samples.tsv -- run this instead of
# `sbatch scripts/4.2_atac_normalize.sb` directly. Needs 3.1_merge_modalities
# already finished for every sample.
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "${P}/config.sh"

n=$(tail -n +2 "${SAMPLE_SHEET}" | grep -c .)
(( n > 0 )) || { echo "ERROR: no sample rows found in ${SAMPLE_SHEET}" >&2; exit 1; }

echo "Submitting 4.2_atac_normalize.sb as an array of ${n} task(s) (--array=0-$((n-1)))"
sbatch --array="0-$((n-1))" "${P}/scripts/4.2_atac_normalize.sb"
