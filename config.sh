#!/bin/bash
# Config for the ISL1 scRNA-seq (10x OCM) pipeline.
# Source this from every sbatch script: `source /pl/active/ZhuLab_BDC/GitHub/ISL1_scRNAseq/config.sh`

set -euo pipefail

# --- Paths -------------------------------------------------------------
# Admera project 260924-0015 (delivered 2026-10-08); this delivery also contains
# ISL1 CUT&RUN fastqs in the same folder -- only the *_scRNA_* files are used here.
RAW_DIR=/pl/active/ZhuLab_BDC/raw_data/ISL1_scRNAseq/7535-0002-10082026_215649
OUT_DIR=/pl/active/ZhuLab_BDC/GitHub/ISL1_scRNAseq
LOGS="${OUT_DIR}/logs"
CELLRANGER_OUT="${OUT_DIR}/results/cellranger"

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
