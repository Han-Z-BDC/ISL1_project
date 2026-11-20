#!/bin/bash
# Step 7.4: submit the D15 + D22 endocrine trajectory chain (user, 2026-11-19). Run from the repo root:
#   bash analysis/scripts/7.4_submit.sh [subset] [velocyto job ids, comma-separated] [stop_after]
# subset      a 7.2 subset name (e.g. alpha_beta_EC; its 7.2 run must have finished); "" or "-" = the full integration
# velocyto    7.3b velocyto jobs (D15 and/or D22) still running; scVelo waits for them. "" or "-" = submit 7.3b D22
#             now (the D15 looms must then already exist)
# stop_after  "slingshot": submit prep + Slingshot only (to look at the lineages before tradeSeq / CellRank)
# amem jobs need --mem >= 240G (QOS minimum).
#   prep -> 7.4a Slingshot -> 7.4e tradeSeq fit (8 chunks) -> 7.4e combine
#   7.3b velocyto + prep -> 7.4c scVelo -> (and 7.4a) 7.4d CellRank
set -euo pipefail
S=analysis/scripts
SUB="${1:-}"; [ "$SUB" = "-" ] && SUB=""
VELO="${2:-}"; [ "$VELO" = "-" ] && VELO=""
STOP="${3:-}"
EXP="--export=ALL${SUB:+,INT_SUBSET=${SUB}}"
K="--kill-on-invalid-dep=yes"
prep=$(sbatch --parsable $EXP "$S/7.4_prep.sb")
sling=$(sbatch --parsable $EXP --dependency=afterok:${prep} $K "$S/7.4a_slingshot.sb")
echo "${SUB:-full}: prep ${prep} | slingshot ${sling}"
[ "$STOP" = "slingshot" ] && exit 0
fit=$(sbatch --parsable $EXP --dependency=afterok:${sling} $K "$S/7.4e_tradeseq_fit.sb")
comb=$(sbatch --parsable $EXP --dependency=afterok:${fit} $K "$S/7.4e_tradeseq_combine.sb")
[ -z "$VELO" ] && VELO=$(sbatch --parsable --job-name=isl1_int_velocyto "$S/7.3b_velocyto.sb" D22)   # isl1_int_ = D15+D22 run
scv=$(sbatch --parsable $EXP --dependency=afterok:${prep}:${VELO//,/:} $K "$S/7.4c_scvelo.sb")
cr=$(sbatch --parsable $EXP --dependency=afterok:${scv}:${sling} $K "$S/7.4d_cellrank.sb")
echo "tradeSeq fit ${fit} (array 0-7) -> combine ${comb} | velocyto ${VELO} -> scVelo ${scv} -> CellRank ${cr}"
