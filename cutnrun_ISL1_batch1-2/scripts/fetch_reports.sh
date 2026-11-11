#!/bin/bash
# Run on the Mac review copy, not on Alpine: pulls code (git) and only the HTML reports
# (MultiQC etc.) from scratch into reports/ (gitignored). No BAMs, bigwigs or peak files.
#   bash scripts/fetch_reports.sh
set -euo pipefail
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
HOST=hzhu10@xsede.org@login.rc.colorado.edu
SCRATCH=/scratch/alpine/hzhu10@xsede.org/ISL1_cutNrun

git -C "$R" pull --ff-only   # remote "alpine" = the repo on PetaLibrary
mkdir -p "$R/reports"
rsync -av -e "ssh -i ${HOME}/.ssh/id_rsa" --prune-empty-dirs \
  --include='*/' --include='*.html' --exclude='*' \
  "${HOST}:${SCRATCH}/results/" "$R/reports/"
