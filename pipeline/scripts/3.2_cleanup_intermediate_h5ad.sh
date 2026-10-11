#!/bin/bash
# Deletes intermediate h5ad files, per sample, ONLY once that sample's
# final merged output already exists (merged/h5ad/{sample}.rna_final.h5ad
# and {sample}.atac_final.h5ad) -- i.e. only after scripts/
# submit_3.1_merge_modalities.sh has succeeded for that sample and you've had
# a chance to review the diagnostic plots/summaries along the way.
#
# This is DESTRUCTIVE and NOT easily reversible -- rerunning from scratch
# means redoing everything each deleted file represents, including
# AMULET (2.1_atac_amulet.py), which is slow. Defaults to a dry run that
# only lists what would be deleted and the total size; pass --yes to
# actually delete.
#
# Deletes, once a sample's merged output exists:
#   rna_qc/h5ad/{sample}.qc_annotated.h5ad
#   rna_qc/h5ad/{sample}.qc_filtered.h5ad
#   remove_background/h5ad/{sample}.soupx_corrected.h5ad
#   detect_doublets/h5ad/{sample}.doublets_annotated.h5ad
#   detect_doublets/h5ad/{sample}.doublets_removed.h5ad   -- superseded by merged/h5ad/{sample}.rna_final.h5ad
#   atac_qc/h5ad/{sample}.import_tsse_checkpoint.h5ad     -- a performance cache, not a decision record
#   atac_qc/h5ad/{sample}.atac_qc_annotated.h5ad
#   atac_qc/h5ad/{sample}.atac_qc_final.h5ad              -- superseded by merged/h5ad/{sample}.atac_final.h5ad
#
# NEVER deletes: any plots/, any summary/, AMULET's amulet_result.csv, or
# anything under merged/ -- the final output this whole thing protects.
#
# Usage:
#   bash scripts/3.2_cleanup_intermediate_h5ad.sh                          # dry run, all samples
#   bash scripts/3.2_cleanup_intermediate_h5ad.sh --yes                    # delete, all samples
#   bash scripts/3.2_cleanup_intermediate_h5ad.sh --sample D22_rep1        # dry run, one sample
#   bash scripts/3.2_cleanup_intermediate_h5ad.sh --sample D22_rep1 --yes  # delete, one sample

set -euo pipefail
P=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "${P}/config.sh"

DO_DELETE=false
ONLY_SAMPLE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --yes) DO_DELETE=true; shift ;;
    --sample) ONLY_SAMPLE="${2:?--sample needs a sample_id}"; shift 2 ;;
    *) echo "ERROR: unrecognized argument '$1'" >&2; exit 1 ;;
  esac
done

human_size() {
  numfmt --to=iec-i --suffix=B "$1" 2>/dev/null || echo "$1 bytes"
}

mapfile -t SAMPLES < <(tail -n +2 "${SAMPLE_SHEET}" | cut -f1)

total_bytes=0
n_ready=0
n_skipped=0

for sample_id in "${SAMPLES[@]}"; do
  if [[ -n "${ONLY_SAMPLE}" && "${sample_id}" != "${ONLY_SAMPLE}" ]]; then
    continue
  fi

  rna_final="${OUT_DIR}/merged/h5ad/${sample_id}.rna_final.h5ad"
  atac_final="${OUT_DIR}/merged/h5ad/${sample_id}.atac_final.h5ad"

  if [[ ! -s "${rna_final}" || ! -s "${atac_final}" ]]; then
    echo "SKIP ${sample_id}: merged output not found -- run scripts/submit_3.1_merge_modalities.sh for this sample first"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  n_ready=$((n_ready + 1))
  candidates=(
    "${OUT_DIR}/rna_qc/h5ad/${sample_id}.qc_annotated.h5ad"
    "${OUT_DIR}/rna_qc/h5ad/${sample_id}.qc_filtered.h5ad"
    "${OUT_DIR}/remove_background/h5ad/${sample_id}.soupx_corrected.h5ad"
    "${OUT_DIR}/detect_doublets/h5ad/${sample_id}.doublets_annotated.h5ad"
    "${OUT_DIR}/detect_doublets/h5ad/${sample_id}.doublets_removed.h5ad"
    "${OUT_DIR}/atac_qc/h5ad/${sample_id}.import_tsse_checkpoint.h5ad"
    "${OUT_DIR}/atac_qc/h5ad/${sample_id}.atac_qc_annotated.h5ad"
    "${OUT_DIR}/atac_qc/h5ad/${sample_id}.atac_qc_final.h5ad"
  )

  for f in "${candidates[@]}"; do
    [[ -f "${f}" ]] || continue
    size=$(stat -c%s "${f}" 2>/dev/null || stat -f%z "${f}" 2>/dev/null || echo 0)
    total_bytes=$((total_bytes + size))
    if [[ "${DO_DELETE}" == true ]]; then
      rm -f "${f}"
      echo "DELETED ${f}"
    else
      echo "would delete: ${f} ($(human_size "${size}"))"
    fi
  done
done

echo
echo "Samples with a confirmed merged final output: ${n_ready}"
echo "Samples skipped (no merged output yet): ${n_skipped}"
if [[ "${DO_DELETE}" == true ]]; then
  echo "Freed: $(human_size "${total_bytes}")"
else
  echo "Would free: $(human_size "${total_bytes}") -- rerun with --yes to actually delete"
fi
