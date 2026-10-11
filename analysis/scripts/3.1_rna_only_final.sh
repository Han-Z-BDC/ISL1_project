#!/bin/bash
# RNA-only stand-in for pipeline step 3.1 (merge_modalities), which reconciles RNA + ATAC
# barcodes for multiome data. This project is RNA-only, so the "final" RNA object per sample
# is simply step 1.3's doublets_removed.h5ad. Pipeline steps 4.1/5.2 read
# ${OUT_DIR}/merged/h5ad/{sample}.rna_final.h5ad, so link it there (no copies).
# Run from the repo root: bash analysis/scripts/3.1_rna_only_final.sh
set -euo pipefail
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
source "$R/study_config/config.sh"

mkdir -p "${OUT_DIR}/merged/h5ad"
tail -n +2 "${SAMPLE_SHEET}" | while IFS=$'\t' read -r s _; do
  src="${OUT_DIR}/detect_doublets/h5ad/${s}.doublets_removed.h5ad"
  [[ -s "${src}" ]] || { echo "ERROR: ${src} missing -- run pipeline step 1.3 first" >&2; exit 1; }
  ln -sfn "${src}" "${OUT_DIR}/merged/h5ad/${s}.rna_final.h5ad"
  echo "${s}: rna_final -> ${src}"
done
