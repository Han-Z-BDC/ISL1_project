#!/bin/bash
# Config for the ISL1 scRNA-seq (10x OCM) project -- single source of every path.
# Sourced by scripts/*.sb (Cell Ranger) and, through the pipeline/config.sh symlink
# (made by study_config/setup_links.sh), by every single-cell-pipeline step in pipeline/.

set -euo pipefail

# --- Paths -------------------------------------------------------------
# Admera project 260924-0015 (delivered 2026-10-08); this delivery also contains
# ISL1 CUT&RUN fastqs in the same folder -- only the *_scRNA_* files are used here.
RAW_DIR=/pl/active/ZhuLab_BDC/raw_data/ISL1_scRNAseq/7535-0002-10082026_215649
REPO_DIR=/pl/active/ZhuLab_BDC/GitHub/ISL1_scRNAseq   # code + configs only, no outputs
# All outputs (logs + results) on scratch -- purged after 90 days, copy keepers to /pl/active.
SCRATCH_DIR=/scratch/alpine/hzhu10@xsede.org/ISL1_scRNAseq
LOGS="${SCRATCH_DIR}/logs"
CELLRANGER_OUT="${SCRATCH_DIR}/cellranger"

# Same reference as the diff-multiome and earlier cellranger count runs, so gene
# IDs/names match across datasets. Cell Ranger accepts the ARC reference for GEX.
REFERENCE=/pl/active/ZhuLab_BDC/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

# --- Tools ---------------------------------------------------------------
# OCM (on-chip multiplexing) needs Cell Ranger >= 9.0.
CELLRANGER=/projects/hzhu10@xsede.org/software/cellranger-9.0.1/bin/cellranger

# --- Samples ---------------------------------------------------------------
# One OCM channel (GEM well) per entry: <run_id>, multi config in configs/<run_id>.multi.csv.
# The fastq_id inside each config is the fastq prefix before _S#_L00#_R#_001.fastq.gz.
RUNS=(
  HZ_Diff34_D40_scRNA
)

# --- single-cell-pipeline (pipeline/, git subtree) -----------------------
# Variable names below are the ones pipeline/scripts/* expect (see pipeline/config.example.sh).
# Cell Ranger multi writes per_sample_outs/<s>/count/sample_{raw,filtered}_feature_bc_matrix.h5;
# the pipeline expects ${IN_DIR}/<s>/outs/{raw,filtered}_feature_bc_matrix.h5, so
# setup_links.sh builds that layout under IN_DIR as symlinks (no copies).
CELLRANGER_SAMPLES="${CELLRANGER_OUT}/HZ_Diff34_D40_scRNA/outs/per_sample_outs"
IN_DIR="${SCRATCH_DIR}/cellranger_std"
OUT_DIR="${SCRATCH_DIR}/out"                              # all pipeline outputs
SAMPLE_SHEET="${REPO_DIR}/study_config/samples.tsv"
CONFIG_DIR="${REPO_DIR}/pipeline/configs"
TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/single.cell/bin   # same env as diff-multiome
