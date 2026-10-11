#!/bin/bash
# Copy this to config.sh and fill in your own paths. config.sh is gitignored
# on purpose -- it's specific to your account and your project, same as in
# cutnrun-pipeline.

# Root of your CellRanger outputs, one subfolder per sample (e.g.
# ${IN_DIR}/<sample_id>/outs/raw_feature_bc_matrix.h5)
IN_DIR=/projects/$USER/your_project/cellranger_outs

# Where this pipeline writes everything it produces
OUT_DIR=/projects/$USER/your_project/sc_pipeline_out

# Sample sheet -- see samples.example.tsv
SAMPLE_SHEET=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/samples.tsv

# Per-step config files (parameters, not paths -- these are versioned, unlike
# this file)
CONFIG_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/configs

# bin/ of the conda env this repo's environment.yml builds -- find it with
# `conda env list` after `conda env create -f environment.yml`
TOOLS=/projects/$USER/software/anaconda/envs/sc-pipeline/bin

# ATAC QC only: a blacklist BED of repeat/artifact regions, excluded from
# both AMULET and SnapATAC2's own feature selection.
BLACKLIST_BED=/projects/$USER/genome/GRCh38_unified_blacklist.bed
