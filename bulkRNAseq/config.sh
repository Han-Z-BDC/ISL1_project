#!/bin/bash
# Config for the ISL1 rescue bulk RNA-seq pipeline.
# Source this from every sbatch script: `source "$(dirname "$0")/../config.sh"`

set -euo pipefail

# --- Paths -------------------------------------------------------------
RAW_DIR=/pl/active/ZhuLab_BDC/raw_data/ISL1_bulk_RNA
OUT_DIR=/pl/active/ZhuLab_BDC/GitHub/ISL1_bulkRNAseq
LOGS="${OUT_DIR}/logs"

# STAR reference: lab-built index, NOT the Cell Ranger ARC-bundled one
# (see docs/star_index_notes.md for why).
GENOME_DIR=/pl/active/ZhuLab_BDC/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
STAR_INDEX="${GENOME_DIR}/STARindex_2.7.5c"
GTF="${GENOME_DIR}/genes/genes.gtf"
RSEM_REF="${GENOME_DIR}/rsem_ref/GRCh38_2020A"  # built by scripts/rsem_prepare_reference.sb

# --- Tools ---------------------------------------------------------------
# Dedicated conda env for this pipeline (built via ../environment.yml):
#   micromamba create -p <TOOLS_ENV> -f environment.yml
TOOLS_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/bulkRNAseq
TOOLS="${TOOLS_ENV}/bin"

# --- Samples ---------------------------------------------------------------
# Single-end. Naming: C7S59W_<idx>_<idx>_<day>_<genotype>_<construct>_rep_<n>.fastq.gz
SAMPLES=(
  C7S59W_1_1_D12_Ctr_GFP_rep_1
  C7S59W_2_2_D12_Ctr_GFP_rep_2
  C7S59W_3_3_D12_KO_GFP_rep_1
  C7S59W_4_4_D12_KO_GFP_rep_2
  C7S59W_5_5_D12_KO_ISL1_rep_1
  C7S59W_6_6_D12_KO_ISL1_rep_2
  C7S59W_7_7_D15_Ctr_GFP_rep_1
  C7S59W_8_8_D15_Ctr_GFP_rep_2
  C7S59W_9_9_D15_KO_GFP_rep_1
  C7S59W_10_10_D15_KO_GFP_rep_2
  C7S59W_11_11_D15_KO_ISL1_rep_1
  C7S59W_12_12_D15_KO_ISL1_rep_2
)
