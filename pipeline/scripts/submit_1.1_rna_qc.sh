#!/bin/bash
# Submits 1.1_rna_qc.sb as a Slurm array sized to however many samples are
# currently in samples.tsv -- run this instead of `sbatch scripts/1.1_rna_qc.sb`
# directly, so the --array range always matches your sample sheet.
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "${P}/config.sh"

n=$(tail -n +2 "${SAMPLE_SHEET}" | grep -c .)
(( n > 0 )) || { echo "ERROR: no sample rows found in ${SAMPLE_SHEET}" >&2; exit 1; }

echo "Submitting 1.1_rna_qc.sb as an array of ${n} task(s) (--array=0-$((n-1)))"
sbatch --array="0-$((n-1))" "${P}/scripts/1.1_rna_qc.sb"
