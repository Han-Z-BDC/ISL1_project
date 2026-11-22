#!/bin/bash
# Step 7.3m: full D15 tradeSeq (beta / alpha / EC), prep -> fit array (one core per chunk) -> combine.
# Run from the repo root (scRNAseq/):  bash analysis/scripts/7.3m_tradeseq_submit.sh [test]
#   test: prep + chunk 0 only (to check time and memory before the full array)
set -euo pipefail
S=/pl/active/ZhuLab_BDC/GitHub/ISL1_project/scRNAseq/analysis/scripts
prep=$(sbatch --parsable "${S}/7.3m_tradeseq_prep.sb")
if [[ "${1:-}" == "test" ]]; then
  fit=$(sbatch --parsable --array=0 --dependency=afterok:${prep} "${S}/7.3m_tradeseq_fit.sb")
  echo "7.3m test: prep ${prep}, chunk 0 ${fit}"; exit 0
fi
fit=$(sbatch --parsable --dependency=afterok:${prep} "${S}/7.3m_tradeseq_fit.sb")
comb=$(sbatch --parsable --dependency=afterok:${fit} "${S}/7.3m_tradeseq_combine.sb")
echo "7.3m: prep ${prep}, fit ${fit}, combine ${comb}"
