#!/bin/bash
# Step 7.2: D15 + D22 endocrine integration on stage. Submits 7.2a -> 7.2b -> 7.2c (GPU) -> 7.2d in order.
# Run from the repo root (scRNAseq/): bash analysis/scripts/7.2_integrate_submit.sh
set -euo pipefail
S=analysis/scripts
a=$(sbatch --parsable $S/7.2a_integrate_prep.sb)
b=$(sbatch --parsable --dependency=afterok:$a --kill-on-invalid-dep=yes $S/7.2b_integrate_sct.sb)
c=$(sbatch --parsable --dependency=afterok:$b --kill-on-invalid-dep=yes $S/7.2c_integrate_scvi.sb)
d=$(sbatch --parsable --dependency=afterok:$c --kill-on-invalid-dep=yes $S/7.2d_integrate_report.sb)
echo "7.2a $a  7.2b $b  7.2c $c  7.2d $d"
