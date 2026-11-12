#!/bin/bash
# One-time setup on Alpine (safe to re-run): lets the vendored pipeline/ run this study's config.
#   pipeline/config.sh -> study_config/config.sh;  pipeline/logs -> ${LOGS}/pipeline
#   batch 3 fastqs (TENX_DELIVERY): ${RAW_LINKS}/<library>_R{1,2}_001.fastq.gz -> raw file
#       (step 01 writes the 10x-filtered copies to ${RAW_DIR})
#   batch 1/2 fastqs: ${RAW_DIR}/<library>_R{1,2}_001.fastq.gz -> raw file (aligned as delivered)
# config.sh and logs/ are gitignored, so pipeline/ stays identical to cutnrun-pipeline.
set -euo pipefail
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$R/study_config/config.sh"

ln -sfn ../study_config/config.sh "$R/pipeline/config.sh"
mkdir -p "${LOGS}/pipeline" "${RAW_LINKS}" "${RAW_DIR}" "${OUT_DIR}"
ln -sfn "${LOGS}/pipeline" "$R/pipeline/logs"

n=0
while IFS=$'\t' read -r _s library _d _r _a _b delivery _o; do
  for mate in R1 R2; do
    f=${library}_${mate}_001.fastq.gz
    src=${RAW_BASE}/${delivery}/${f}
    [[ -s "${src}" ]] || { echo "ERROR: missing ${src}" >&2; exit 1; }
    if [[ "${delivery}" == "${TENX_DELIVERY}" ]]; then ln -sfn "${src}" "${RAW_LINKS}/${f}"; else ln -sfn "${src}" "${RAW_DIR}/${f}"; fi
    n=$((n + 1))
  done
done < <(tail -n +2 "${METADATA}")
echo "${n} fastq links: $(ls "${RAW_LINKS}" | wc -l) for the 10x filter (${RAW_LINKS}), $(ls "${RAW_DIR}" | wc -l) already in ${RAW_DIR}"
echo "next: bash scripts/run_all.sh"
