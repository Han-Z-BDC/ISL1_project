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
SAMPLE_SHEET=${STUDY_DIR}/study_config/samples.tsv

# IgG background: every ISL1 sample is called both ways and the peak sets compared
# (scripts/03b_compare_igg_methods.sb); a sample whose peaks change a lot flags its own IgG.
# Pick the method per job, default matched:
#   cd pipeline && sbatch --export=ALL,IGG_METHOD=pooled scripts/03_call_macs3.sb
# Each method gets its own OUT_DIR so igg_background/ and macs3/ don't overwrite each other;
# results_pooled/bam is a symlink to results/bam (setup_links.sh), so alignment runs once.
IGG_METHOD=${IGG_METHOD:-matched}
case "${IGG_METHOD}" in
  matched) OUT_DIR=${SCRATCH}/results ;;
  pooled)  OUT_DIR=${SCRATCH}/results_pooled ;;
  *) echo "ERROR: IGG_METHOD must be matched or pooled, got '${IGG_METHOD}'" >&2; exit 1 ;;
esac

# hg38 (UCSC) bowtie2 index and the GRCh38 unified blacklist already on PetaLibrary.
BOWTIE2_INDEX=/pl/active/ZhuLab_BDC/genome/hg38_bt2/hg38
BLACKLIST=/pl/active/ZhuLab_BDC/genome/GRCh38_unified_blacklist.bed

# Env built from pipeline/environment.yml by scripts/00b_build_env.sb (job 34348808);
# resolved versions in study_config/cutnrun_env.export.yml.
TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/cutnrun/bin

MAPQ=20
MACS_QVALUE=0.01

# Pool = all 13 IgGs (every row with role IgG), each deduplicated and downsampled to the same
# POOL_TARGET_FRAGMENTS, one fixed seed each. The target must be <= the smallest IgG's
# deduplicated fragment count, and strictly below it (a fraction of exactly 1 breaks
# samtools -s SEED.FRAC). Smallest after 01_align (job 34354886): D12_KO_rep1_IgG,
# 9,350,997 deduplicated fragments (flagstat: (total - duplicates) / 2). Pool = 13 x 9 M.
mapfile -t POOL_IGG_SAMPLES < <(awk -F'\t' 'NR>1 && $3=="IgG" {print $2}' "${SAMPLE_SHEET}")
POOL_IGG_SEEDS=($(seq 101 $((100 + ${#POOL_IGG_SAMPLES[@]}))))
POOL_TARGET_FRAGMENTS=9000000
# Final peak calls (study-only, scripts/02c + 03c): one IgG background per biological group
# (day x genotype, from sample_metadata.tsv) = that group's deduplicated IgGs merged, used for
# every ISL1 replicate in the group. D12_KO_rep1_IgG is excluded: it is the one sample flagged
# by the pre-specified matched-vs-pooled check (28 peaks vs own IgG, 10,803 vs pool; job 34388152)
# and the shallowest, least complex IgG (9.35 M unique fragments, 55% duplicates).
METADATA=${STUDY_DIR}/study_config/sample_metadata.tsv
GROUP_OUT=${SCRATCH}/results_group
GROUP_IGG_EXCLUDE=(D12_KO_rep1_IgG)

if [[ "${IGG_METHOD}" == pooled ]]; then
  [[ -n "${POOL_TARGET_FRAGMENTS}" ]] || { echo "ERROR: POOL_TARGET_FRAGMENTS not set in study_config/config.sh" >&2; exit 1; }
  (( ${#POOL_IGG_SAMPLES[@]} >= 2 )) || { echo "ERROR: need >=2 IgG rows in ${SAMPLE_SHEET} to pool" >&2; exit 1; }
fi
