#!/bin/bash
# Fast, interactive threshold-tuning helper for ATAC QC -- run this from an
# acompile session (NOT sbatch) to iterate on configs/atac_qc.yaml's
# thresholds for one sample without paying for the full pipeline (tile
# matrix, feature selection, all three doublet detectors -- especially
# AMULET, which is slow).
#
# Usage: bash scripts/qc_only_atac.sh <sample_id>
#
# Look at {sample_id}.threshold_review.png (and the printed pass/fail
# count) in $OUT_DIR/atac_qc/plots -- if the cutoff lines don't look
# right, edit configs/atac_qc.yaml and rerun this. The expensive part
# (parsing the fragment file) is cached to
# {sample_id}.import_tsse_checkpoint.h5ad after the first run, so repeat
# runs here stay fast -- only delete that file if the fragment file
# itself changed. Once the thresholds look right, run the full pipeline
# with scripts/submit_2.2_atac_qc.sh (which does NOT pass --qc-only).
set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "${P}/config.sh"

sample_id="${1:?Usage: bash scripts/qc_only_atac.sh <sample_id>}"

row=$(awk -F'\t' -v s="${sample_id}" '$1==s {print; f=1} END{if(!f) exit 1}' "${SAMPLE_SHEET}") \
  || { echo "ERROR: sample_id '${sample_id}' not found in ${SAMPLE_SHEET}" >&2; exit 1; }
h5_rel=$(echo "${row}" | cut -f2)

outs_dir=$(dirname "${h5_rel}")
frag_path="${IN_DIR}/${outs_dir}/atac_fragments.tsv.gz"
[[ -s "${frag_path}" ]] || { echo "ERROR: ${frag_path} not found -- check samples.tsv / IN_DIR layout" >&2; exit 1; }

mkdir -p "${OUT_DIR}/atac_qc"

"${TOOLS}/python" "${P}/scripts/2.2_atac_qc.py" \
  --sample-id "${sample_id}" \
  --frag-path "${frag_path}" \
  --blacklist-bed "${BLACKLIST_BED}" \
  --config "${CONFIG_DIR}/atac_qc.yaml" \
  --outdir "${OUT_DIR}/atac_qc" \
  --qc-only
