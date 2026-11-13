#!/bin/bash
# Config for the ISL1 scRNA-seq (10x OCM) project -- single source of every path.
# Sourced by scripts/*.sb (Cell Ranger) and, through the pipeline/config.sh symlink
# (made by study_config/setup_links.sh), by every single-cell-pipeline step in pipeline/.

set -euo pipefail

# --- Paths -------------------------------------------------------------
# Admera project 260924-0015 (delivered 2026-10-08). The delivery also had 32 ISL1 CUT&RUN
# fastqs; moved 2026-10-16 to /pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run/7535-0002-10082026_215649
# (same folder name). Here: the 2 scRNA fastqs, md5.txt (their 2 lines) and
# md5.admera_full_delivery.txt (Admera's original, all 34 files).
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
# Admera project 261016-0008 (delivered 2026-11-02, raw_data/ISL1_scRNAseq/7535-0003-11022026_200951):
# extra sequencing of the same D40 library (S268, L003) plus two new OCM channels, Diff37 D15
# and D22. All channels: OB1 = ctr rep1, OB2 = KO rep1, OB3 = ctr rep2, OB4 = KO rep2.
# HZ_Diff34_D40_scRNA = first D40 sequencing only (the run the current out/ analysis uses);
# HZ_Diff34_D40_scRNA_combined = same library, both sequencing runs (two [libraries] rows).
RUNS=(
  HZ_Diff34_D40_scRNA_combined
  HZ_Diff37_D15_scRNA
  HZ_Diff37_D22_scRNA
)

# --- single-cell-pipeline (pipeline/, git subtree) -----------------------
# Variable names below are the ones pipeline/scripts/* expect (see pipeline/config.example.sh).
# Cell Ranger multi writes per_sample_outs/<s>/count/sample_{raw,filtered}_feature_bc_matrix.h5;
# the pipeline expects ${IN_DIR}/<s>/outs/{raw,filtered}_feature_bc_matrix.h5, so
# setup_links.sh builds that layout under IN_DIR as symlinks (no copies).
# 2026-11-08: all 12 samples (D40 from the combined run, plus D15 + D22). setup_links.sh looks
# each sample up in these per_sample_outs dirs in order. The D40-only analysis (first D40
# sequencing, run HZ_Diff34_D40_scRNA) is left as it was in cellranger_std/ and out/.
CELLRANGER_SAMPLE_DIRS=(
  "${CELLRANGER_OUT}/HZ_Diff34_D40_scRNA_combined/outs/per_sample_outs"
  "${CELLRANGER_OUT}/HZ_Diff37_D15_scRNA/outs/per_sample_outs"
  "${CELLRANGER_OUT}/HZ_Diff37_D22_scRNA/outs/per_sample_outs"
)
IN_DIR="${SCRATCH_DIR}/cellranger_std_all"
OUT_DIR="${SCRATCH_DIR}/out_all"                          # all pipeline outputs
SAMPLE_SHEET="${REPO_DIR}/study_config/samples.tsv"
CONFIG_DIR="${REPO_DIR}/pipeline/configs"
TOOLS=/projects/hzhu10@xsede.org/software/anaconda/envs/single.cell/bin   # same env as diff-multiome
# single.cell needs these when its binaries are called by full path (no conda activate in
# jobs): otherwise pandas loads the system libstdc++ (GLIBCXX_3.4.29 not found), and the
# rpy2-based steps (SoupX, scDblFinder, SCTransform) need R_HOME. Same as diff-multiome.
_ENV_LIB=/projects/hzhu10@xsede.org/software/anaconda/envs/single.cell/lib
export LD_LIBRARY_PATH="${_ENV_LIB}:${LD_LIBRARY_PATH:-}"
export R_HOME="${_ENV_LIB}/R"

# --- Cross-sample merge + integration (pipeline 5.2-5.4) ------------------
# The pipeline parses each sample name as {diff_batch}_{stage}_{tech_batch}. Mapping here
# (rename_manifest.csv): diff_batch = Diff34 / Diff37 (differentiation), "stage" = <day>-ctr /
# <day>-KO (the biology to PRESERVE -- Harmony never corrects on it), tech_batch = rep1 / rep2.
# Each day is one OCM channel (4 samples), so channel = day. Diff34 is only D40 and Diff37 only
# D15 + D22, so diff_batch is confounded with day. rep1/rep2 are biological replicates that
# are NOT paired (not across days, not ctr vs KO within a day), so tech_batch is never a batch
# variable: each of the 12 samples is one replicate (pseudobulk unit for DE/composition).
# Integration (user, 2026-11-08) is run three ways, Harmony (analysis 5.4c) and scVI (5.5c):
# none / diff_batch only / per sample -- compared in 5.4b. The pipeline's 5.4/5.5 A/B strategies
# (diff_batch + tech_batch) are not used for this study. 6.1 clusters on log unintegrated.
INTEGRATION_SAMPLES_FILE="${REPO_DIR}/study_config/integration/samples_for_integration.txt"
INTEGRATION_SAMPLE_MAP_CSV="${REPO_DIR}/study_config/rename_manifest.csv"
HARMONY_THETA_DIFF_BATCH=3.0
HARMONY_THETA_TECH_BATCH=0.5
