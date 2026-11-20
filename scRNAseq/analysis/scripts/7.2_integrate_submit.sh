#!/bin/bash
# Step 7.2: D15 + D22 endocrine integration on stage. Submits 7.2a -> 7.2b -> 7.2c (GPU) -> 7.2d in order.
# Run from the repo root (scRNAseq/): bash analysis/scripts/7.2_integrate_submit.sh [subset]
# With a subset name, every step gets INT_SUBSET=<subset>: cells limited to the labels in
# study_config/annotation/integration_subset_<subset>.txt, outputs in day_analysis/integration_D15_D22_endocrine_<subset>/.
set -euo pipefail
S=analysis/scripts
SUB="${1:-}"
if [ -n "$SUB" ]; then
  test -s "study_config/annotation/integration_subset_${SUB}.txt"
  EXP="--export=ALL,INT_SUBSET=${SUB}"
else
  EXP="--export=ALL"
fi
a=$(sbatch --parsable $EXP $S/7.2a_integrate_prep.sb)
b=$(sbatch --parsable $EXP --dependency=afterok:$a --kill-on-invalid-dep=yes $S/7.2b_integrate_sct.sb)
c=$(sbatch --parsable $EXP --dependency=afterok:$b --kill-on-invalid-dep=yes $S/7.2c_integrate_scvi.sb)
d=$(sbatch --parsable $EXP --dependency=afterok:$c --kill-on-invalid-dep=yes $S/7.2d_integrate_report.sb)
echo "7.2a $a  7.2b $b  7.2c $c  7.2d $d  ${SUB:+(subset $SUB)}"
