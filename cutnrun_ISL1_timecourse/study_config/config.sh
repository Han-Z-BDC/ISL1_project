#!/bin/bash
# ISL1 CUT&RUN time course (control, D15 / D18 / D22, three replicates per day): study config.
# setup_links.sh symlinks pipeline/config.sh -> this file, so the vendored pipeline/
# (cutnrun-pipeline @1afe307) reads it unchanged. Variable names follow pipeline/config.example.sh;
# the extra ones are study-only.

set -euo pipefail

STUDY_DIR=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)
PIPELINE_DIR=${STUDY_DIR}/pipeline

# --- Study-only paths ----------------------------------------------------
# Admera deliveries (raw, untouched): batch 1 (2026-08), batch 2 (2026-10), batch 3 (2026-11).
RAW_BASE=/pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run
TENX_DELIVERY=7535-0003-11022026_200951     # 10x reads leaked into this delivery; filtered (step 01)
SCRATCH=/scratch/alpine/hzhu10@xsede.org/ISL1_cutnrun_timecourse
LOGS=${SCRATCH}/logs
METADATA=${STUDY_DIR}/study_config/sample_metadata.tsv
RAW_LINKS=${SCRATCH}/raw_links              # links to the batch 3 fastqs (input of step 01)
# 10x 3' GEX libraries sequenced on the batch 3 lane: their cell barcodes are a QC column of step 01.
TENX_FASTQ_DIR=/pl/active/ZhuLab_BDC/raw_data/ISL1_scRNAseq/7535-0003-11022026_200951
# Analysis jobs: amem, 12 CPUs, 240 GB (partition minimum is 240 GB).
AMEM=(-p amem --qos=mem-normal --cpus-per-task=12 --mem=240G)

# --- Pipeline inputs (names as in pipeline/config.example.sh) ------------
# One fastq folder: 10x-filtered batch 3 fastqs (step 01) and links to the batch 1/2 fastqs.
RAW_DIR=${SCRATCH}/fastq
SAMPLE_SHEET=${STUDY_DIR}/study_config/samples.tsv
OUT_DIR=${SCRATCH}/results
IGG_METHOD=matched                           # each library against its own deduplicated IgG

# hg38 (UCSC) bowtie2 index and the GRCh38 unified blacklist on PetaLibrary.
BOWTIE2_INDEX=/pl/active/ZhuLab_BDC/genome/hg38_bt2/hg38
BLACKLIST=/pl/active/ZhuLab_BDC/genome/GRCh38_unified_blacklist.bed

# Environments: cutnrun (pipeline/environment.yml; also deepTools, bedtools, python), featureCounts,
# DiffBind / DESeq2 / GreyListChIP.
TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/cutnrun/bin
FEATURECOUNTS=/projects/hzhu10@xsede.org/software/anaconda/envs/isl1_downstream/bin/featureCounts
DIFFBIND_TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/diffbind/bin

MAPQ=20
MACS_QVALUE=0.01
SITE_HALF_WIDTH=200                          # differential sites: summit +/- 200 bp = 401 bp
TRACK_BIN=10                                 # bigWig bin (bp)
