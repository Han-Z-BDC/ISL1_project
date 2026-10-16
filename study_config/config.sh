#!/bin/bash
# ISL1 CUT&RUN study config. setup_links.sh symlinks pipeline/config.sh -> this file, so the
# vendored pipeline/ (cutnrun-pipeline) reads it unchanged. Variable names follow
# pipeline/config.example.sh; the extra ones at the top are study-only.

set -euo pipefail

STUDY_DIR=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)
PIPELINE_DIR=${STUDY_DIR}/pipeline

# --- Study-only paths ----------------------------------------------------
# Two Admera deliveries (both verified against md5.txt, see scripts/00_md5_check.sb):
#   batch 1  7535-0001-08202026_144755  rep1 of D12 Ctr/KO, D15 Ctr/KO, D22 Ctr  (10 libraries)
#   batch 2  7535-0002-10082026_215649  later reps (16 libraries); moved here from the
#            scRNA delivery folder 2026-10-16
RAW_BASE=/pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run
RAW_BATCHES=(7535-0001-08202026_144755 7535-0002-10082026_215649)
# All outputs (logs + results) on scratch -- purged after 90 days, copy keepers to /pl/active.
SCRATCH=/scratch/alpine/hzhu10@xsede.org/ISL1_cutNrun
LOGS=${SCRATCH}/logs

# --- Pipeline inputs (names as in pipeline/config.example.sh) ------------
# The pipeline reads one RAW_DIR; setup_links.sh fills it with symlinks to both batches.
RAW_DIR=${SCRATCH}/raw_links
OUT_DIR=${SCRATCH}/results
SAMPLE_SHEET=${STUDY_DIR}/study_config/samples.tsv

# hg38 (UCSC) bowtie2 index and the GRCh38 unified blacklist already on PetaLibrary.
BOWTIE2_INDEX=/pl/active/ZhuLab_BDC/genome/hg38_bt2/hg38
BLACKLIST=/pl/active/ZhuLab_BDC/genome/GRCh38_unified_blacklist.bed

# Env built from pipeline/environment.yml (not yet built as of 2026-10-16: no existing env
# has the full tool set -- cutrun_tools lacks fastqc/multiqc/cutadapt/bowtie2).
TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/cutnrun/bin

MAPQ=20
MACS_QVALUE=0.01

# Each ISL1 library has its own IgG from the same chromatin prep -> matched.
IGG_METHOD=matched
POOL_TARGET_FRAGMENTS=4000000
POOL_IGG_SAMPLES=()
POOL_IGG_SEEDS=()
