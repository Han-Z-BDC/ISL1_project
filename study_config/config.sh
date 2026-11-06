#!/bin/bash
# ISL1 KO / H3 control CUT&RUN time course (Admera delivery 7535-0003, 2026-11-02): study config.
# setup_links.sh symlinks pipeline/config.sh -> this file, so the vendored pipeline/
# (cutnrun-pipeline) reads it unchanged. Variable names follow pipeline/config.example.sh; the
# extra ones at the top are study-only.

set -euo pipefail

STUDY_DIR=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)
PIPELINE_DIR=${STUDY_DIR}/pipeline

# --- Study-only paths ----------------------------------------------------
# One Admera delivery, verified against its md5.txt after download (md5_check.txt: 130/130 OK);
# the 124 CUT&RUN fastqs were moved here from the scRNA folder on 2026-11-05.
RAW_BASE=/pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run
RAW_BATCHES=(7535-0003-11022026_200951)
# All outputs (logs + results) on scratch -- purged after 90 days, copy keepers to /pl/active.
SCRATCH=/scratch/alpine/hzhu10@xsede.org/ISL1_KO_histone_cutNrun
LOGS=${SCRATCH}/logs
METADATA=${STUDY_DIR}/study_config/sample_metadata.tsv

# --- Pipeline inputs (names as in pipeline/config.example.sh) ------------
RAW_DIR=${SCRATCH}/raw_links
SAMPLE_SHEET=${STUDY_DIR}/study_config/samples.tsv

# IgG background method for the pipeline's steps 02/03 (matched = each sample's own block IgG).
IGG_METHOD=${IGG_METHOD:-matched}
case "${IGG_METHOD}" in
  matched) OUT_DIR=${SCRATCH}/results ;;
  pooled)  OUT_DIR=${SCRATCH}/results_pooled ;;
  *) echo "ERROR: IGG_METHOD must be matched or pooled, got '${IGG_METHOD}'" >&2; exit 1 ;;
esac

# hg38 (UCSC) bowtie2 index and the GRCh38 unified blacklist already on PetaLibrary.
BOWTIE2_INDEX=/pl/active/ZhuLab_BDC/genome/hg38_bt2/hg38
BLACKLIST=/pl/active/ZhuLab_BDC/genome/GRCh38_unified_blacklist.bed

# Same env as ISL1_cutNrun, built from pipeline/environment.yml (cutnrun-pipeline @1afe307).
TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/cutnrun/bin

MAPQ=20
MACS_QVALUE=0.01

# Reproducible peaks (scripts/04-05, ENCODE design). Relaxed calls: p 0.01, top N peaks
# (ENCODE: 300k for TFs, 500k for histone marks). IDR env as built for ISL1_cutNrun.
REPRO_LINES="H3CTR ISL1KO"       # D15 CtrGFP / KOISL1 set aside (2026-11-06)
REPRO_ANTIBODIES="ISL1 H3K27ac H3K27me3"
RELAXED_P=0.01
CAP_TF=300000
CAP_HISTONE=500000
IDR_TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/idr/bin
IDR_THRESHOLD=0.05

# Pooled IgG (pipeline 02, IGG_METHOD=pooled): not used yet; set POOL_TARGET_FRAGMENTS from the
# alignment metrics before using it.
mapfile -t POOL_IGG_SAMPLES < <(awk -F'\t' 'NR>1 && $3=="IgG" {print $2}' "${SAMPLE_SHEET}")
POOL_IGG_SEEDS=($(seq 101 $((100 + ${#POOL_IGG_SAMPLES[@]}))))
POOL_TARGET_FRAGMENTS=
if [[ "${IGG_METHOD}" == pooled ]]; then
  [[ -n "${POOL_TARGET_FRAGMENTS}" ]] || { echo "ERROR: POOL_TARGET_FRAGMENTS not set in study_config/config.sh" >&2; exit 1; }
fi
