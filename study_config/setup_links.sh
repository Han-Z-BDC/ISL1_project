#!/bin/bash
# One-time setup after cloning (safe to re-run): lets the vendored pipeline/ run this study's config.
#   pipeline/config.sh -> study_config/config.sh
#   pipeline/logs      -> ${LOGS}/pipeline   (pipeline .sb files log to logs/ relative to pipeline/)
#   ${IN_DIR}/<sample>/outs/{raw,filtered}_feature_bc_matrix.h5 -> Cell Ranger multi per-sample h5
# config.sh and logs/ are gitignored inside pipeline/, so the subtree stays identical to single-cell-pipeline.
set -euo pipefail
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$R/study_config/config.sh"

ln -sfn ../study_config/config.sh "$R/pipeline/config.sh"
mkdir -p "${LOGS}/pipeline" "${OUT_DIR}"
ln -sfn "${LOGS}/pipeline" "$R/pipeline/logs"

tail -n +2 "${SAMPLE_SHEET}" | while IFS=$'\t' read -r s _; do
  src="${CELLRANGER_SAMPLES}/${s}/count"
  [[ -s "${src}/sample_raw_feature_bc_matrix.h5" ]] || { echo "ERROR: ${src}/sample_raw_feature_bc_matrix.h5 missing" >&2; exit 1; }
  mkdir -p "${IN_DIR}/${s}/outs"
  ln -sfn "${src}/sample_raw_feature_bc_matrix.h5"      "${IN_DIR}/${s}/outs/raw_feature_bc_matrix.h5"
  ln -sfn "${src}/sample_filtered_feature_bc_matrix.h5" "${IN_DIR}/${s}/outs/filtered_feature_bc_matrix.h5"
done
echo "links in place; pipeline outputs go to ${OUT_DIR}"
