#!/bin/bash
# Config for the ISL1 KO time-course bulk RNA-seq (run YJZ5LS).
# Sourced by every script in pipeline/.

set -euo pipefail

# --- Paths -------------------------------------------------------------
RAW_DIR=/pl/active/ZhuLab_BDC/raw_data/ISL1_bulk_RNA
OUT_DIR=/pl/active/ZhuLab_BDC/GitHub/ISL1_project/bulkRNAseq_ISL1_KO_timecourse
LOGS="${OUT_DIR}/logs"

# Same references as bulkRNAseq/ (see ../bulkRNAseq/docs/star_index_notes.md).
GENOME_DIR=/pl/active/ZhuLab_BDC/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
STAR_INDEX="${GENOME_DIR}/STARindex_2.7.5c"
GTF="${GENOME_DIR}/genes/genes.gtf"

# --- Tools ---------------------------------------------------------------
TOOLS_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/bulkRNAseq
TOOLS="${TOOLS_ENV}/bin"
R_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/r_env

# --- Samples ---------------------------------------------------------------
# Single-end, ~94 bp. Fastq basenames (no .fastq.gz). D9 samples (1, 2, 9, 10)
# are not part of this analysis. Metadata: sample_sheet.csv.
SAMPLES=(
  YJZ5LS_3_Ctr-H3_D15_rep1
  YJZ5LS_4_Ctr-H3_D15_rep2
  YJZ5LS_5_Ctr-H3_D18_rep1
  YJZ5LS_6_Ctr-H3_D18_rep2
  YJZ5LS_7_Ctr-H3_D22_rep1
  YJZ5LS_8_Ctr-H3_D22_rep2
  YJZ5LS_11_ISL1-KO_D15_rep1
  YJZ5LS_12_ISL1-KO_D15_rep2
  YJZ5LS_13_ISL1-KO_D18_rep1
  YJZ5LS_14_ISL1-KO_D18_rep2
  YJZ5LS_15_ISL1-KO_D22_rep1
  YJZ5LS_16_ISL1-KO_D22_rep2
)
