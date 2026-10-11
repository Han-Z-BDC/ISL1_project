#!/bin/bash
# Submits 2.1_atac_amulet.sb as a Slurm array sized to however many samples
# are currently in samples.tsv. Run this FIRST, before submit_2.2_atac_qc.sh
# -- AMULET is slow (see 2.1_atac_amulet.py), so start it while you're
# still tuning general-QC thresholds with qc_only_atac.sh, not after.
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "${P}/config.sh"

n=$(tail -n +2 "${SAMPLE_SHEET}" | grep -c .)
(( n > 0 )) || { echo "ERROR: no sample rows found in ${SAMPLE_SHEET}" >&2; exit 1; }

echo "Submitting 2.1_atac_amulet.sb as an array of ${n} task(s) (--array=0-$((n-1)))"
sbatch --array="0-$((n-1))" "${P}/scripts/2.1_atac_amulet.sb"
