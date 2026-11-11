#!/bin/bash
# Copy this file to config.sh (in the same directory) and edit every path below
# for your own project. config.sh is gitignored on purpose -- it's yours, it
# will contain paths specific to your account, and it should never be committed.
#
#   cp config.example.sh config.sh
#   nano config.sh

set -euo pipefail

PIPELINE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ---- Your data --------------------------------------------------------
RAW_DIR=/path/to/your/raw_fastqs
OUT_DIR=/path/to/your/output_directory
SAMPLE_SHEET=${PIPELINE_DIR}/samples.tsv

# ---- Reference genome ---------------------------------------------------
# A pre-built bowtie2 index (prefix, no .bt2 extension) and an ENCODE-style
# blacklist BED for your genome build. Ask your lab where these already live
# before building your own -- it's likely already on /pl/active/ somewhere.
BOWTIE2_INDEX=/path/to/genome/bowtie2_index/genome
BLACKLIST=/path/to/genome/blacklist.bed

# ---- Software -------------------------------------------------------------
# Built from environment.yml: module load miniforge; conda env create -f environment.yml
# Point this at that environment's bin/ directory.
TOOLS=/projects/$USER/software/anaconda/envs/cutnrun/bin

# ---- Alignment / filtering --------------------------------------------
MAPQ=20
MACS_QVALUE=0.01

# ---- IgG background method ----------------------------------------------
# Pick ONE. Either way, whatever IgG BAM ends up as the MACS3 control is
# deduplicated first -- only the target (TF) sample keeps duplicates.
#
#   matched  Each target sample uses its own paired IgG control (the
#            control_sample column in samples.tsv), deduplicated.
#            Simple, no extra samples needed, works fine with a single
#            good-quality IgG replicate per target sample.
#
#   pooled   Several IgG replicates are each deduplicated, downsampled to a
#            common fragment depth, and merged into one shared background
#            used for every target sample. More robust than any single IgG
#            replicate, and useful if you have (or suspect) a weak/noisy
#            individual IgG control -- but needs >=2 usable IgG replicates.
IGG_METHOD=matched

# Only read when IGG_METHOD=pooled:
POOL_TARGET_FRAGMENTS=4000000
POOL_IGG_SAMPLES=(IgG_sample_1 IgG_sample_2 IgG_sample_3)
POOL_IGG_SEEDS=(101 102 103)   # one fixed downsampling seed per sample, same order
