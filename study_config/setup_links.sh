#!/bin/bash
# One-time setup after cloning (safe to re-run): lets the vendored pipeline/ run this study's config.
#   pipeline/config.sh -> study_config/config.sh
#   pipeline/logs      -> ${LOGS}/pipeline   (pipeline .sb files log to logs/ relative to pipeline/)
#   ${RAW_LINKS}/<library>_R{1,2}_001.fastq.gz -> the fastq in the Admera delivery folder
# (the pipeline reads ${RAW_DIR}, the 10x-filtered copies written by scripts/00f_filter_10x.sb)
# config.sh and logs/ are gitignored inside pipeline/, so it stays identical to cutnrun-pipeline.
set -euo pipefail
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IGG_METHOD=matched
source "$R/study_config/config.sh"

ln -sfn ../study_config/config.sh "$R/pipeline/config.sh"
mkdir -p "${LOGS}/pipeline" "${OUT_DIR}/bam" "${RAW_LINKS}"
ln -sfn "${LOGS}/pipeline" "$R/pipeline/logs"

n=0
while IFS=$'\t' read -r library _; do
  for mate in R1 R2; do
    f=${library}_${mate}_001.fastq.gz
    hits=()
    for b in "${RAW_BATCHES[@]}"; do
      [[ -s "${RAW_BASE}/${b}/${f}" ]] && hits+=("${RAW_BASE}/${b}/${f}")
    done
    (( ${#hits[@]} == 1 )) || { echo "ERROR: ${f} found in ${#hits[@]} batches (need exactly 1)" >&2; exit 1; }
    ln -sfn "${hits[0]}" "${RAW_LINKS}/${f}"
    n=$((n + 1))
  done
done < <(tail -n +2 "${SAMPLE_SHEET}")

echo "${n} fastq links in ${RAW_LINKS}; pipeline outputs go to ${OUT_DIR}"
echo "next: sbatch --array=0-61%20 scripts/00f_filter_10x.sb (writes ${RAW_DIR}), then pipeline/scripts/validate_samples.sh"
