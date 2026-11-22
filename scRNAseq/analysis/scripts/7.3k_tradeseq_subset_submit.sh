#!/bin/bash
# Step 7.3j + 7.3k: alpha / beta / EC driver table and gene set (7.3j, seconds, run here), then the D15 subset tradeSeq:
# 16 fit chunks (array) and the combine step after all of them. Run from the repo root (scRNAseq/):
#   bash analysis/scripts/7.3k_tradeseq_subset_submit.sh
set -euo pipefail
source /pl/active/ZhuLab_BDC/GitHub/ISL1_project/scRNAseq/study_config/config.sh
export OUT_DIR REPO_DIR
LD_LIBRARY_PATH="$(dirname "${TOOLS}")/lib" "${TOOLS}/python" "${REPO_DIR}/analysis/scripts/7.3j_driver_table.py"
fit=$(sbatch --parsable "${REPO_DIR}/analysis/scripts/7.3k_tradeseq_subset_fit.sb")
comb=$(sbatch --parsable --dependency=afterok:${fit} "${REPO_DIR}/analysis/scripts/7.3k_tradeseq_subset_combine.sb")
echo "7.3k fit ${fit} (16 chunks), combine ${comb}"
