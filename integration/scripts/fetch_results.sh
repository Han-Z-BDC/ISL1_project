#!/bin/bash
# Run on the Mac from integration/: copy the small result tables from Alpine into results/ for build_report.py.
#   bash scripts/fetch_results.sh
set -euo pipefail
cd "$(dirname "$0")/.."
REMOTE=hzhu10@xsede.org@login.rc.colorado.edu
SRC=/pl/active/ZhuLab_BDC/GitHub/ISL1_project/integration/results
for d in 01_great 02_direct_targets 03_binding_change 04_candidates 05_binding_signal; do
  mkdir -p "results/${d}"
done
scp -q -i ~/.ssh/id_rsa "${REMOTE}:${SRC}/01_great/summary.tsv" results/01_great/
scp -q -i ~/.ssh/id_rsa "${REMOTE}:${SRC}/02_direct_targets/*.tsv" results/02_direct_targets/
scp -q -i ~/.ssh/id_rsa "${REMOTE}:${SRC}/03_binding_change/*.tsv" results/03_binding_change/
scp -q -i ~/.ssh/id_rsa "${REMOTE}:${SRC}/04_candidates/*.tsv" results/04_candidates/
scp -q -i ~/.ssh/id_rsa "${REMOTE}:${SRC}/05_binding_signal/group_signal.tsv" "${REMOTE}:${SRC}/05_binding_signal/wilcoxon.tsv" results/05_binding_signal/
ls -l results/*/
