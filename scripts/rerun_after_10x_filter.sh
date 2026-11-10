#!/bin/bash
# Re-run of the settled analysis on the 10x-filtered fastqs (agreed 2026-11-10). Run from the repo
# root on a login node after step 00f has finished (or pass its job id) and after the old results
# were moved aside:  mv ${SCRATCH}/results ${SCRATCH}/results_with10x
#   bash scripts/rerun_after_10x_filter.sh [<00f job id>]
# Every job runs on amem (12 cores, 240 GB); parameters as in the original runs:
#   01 align, 02 IgG background (pipeline/, unchanged), 03 peaks, 04 relaxed + pseudoreps,
#   05 reproducible sets, 13 / 14 earlier control ISL1 (D18 rep4 is a batch 3 library),
#   15 three-replicate sets (ctr_trio), 18 kept ISL1 sets, 16 / 17 / 19 / 21 histone differential
#   (H3K27ac: full union, drop-one D18 H3CTR rep2, 500 bp windows; H3K27me3: 3 kb domains,
#   drop-one D22 H3CTR rep2, 5 kb windows, broad DiffBind calls), 23a-c tracks, 24 comparison.
# Not re-run (exploratory; old outputs stay in results_with10x): 06-12, 20, 22, ctr_multi.
set -euo pipefail
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd); cd "${R}"
IGG_METHOD=matched
source study_config/config.sh
[[ -d "${SCRATCH}/results_with10x" && ! -e "${OUT_DIR}/bam" ]] || { echo "ERROR: move ${OUT_DIR} to results_with10x first" >&2; exit 1; }
AM=(-p amem --qos=mem-normal --cpus-per-task=12 --mem=240G)
dep () { local d=""; for j in "$@"; do [[ -n "${j}" ]] && d+=":${j}"; done; [[ -n "${d}" ]] && echo "--dependency=afterok${d}"; }
sb () { sbatch --parsable "${AM[@]}" "$@"; }

F=${1:-}
cd pipeline
A=$(sb --array=0-61%20 --time=12:00:00 $(dep "${F}") scripts/01_align.sb)
G=$(sb --time=12:00:00 $(dep "${A}") scripts/02_build_igg_background.sb)
cd ..
P3=$(sb --array=0-43 --time=12:00:00 $(dep "${G}") scripts/03_call_peaks.sb)
P4=$(sb --array=0-39 --time=12:00:00 $(dep "${G}") scripts/04_relaxed_pseudoreps.sb)
P5=$(sb --array=0-19 --time=12:00:00 $(dep "${P3}" "${P4}") scripts/05_reproducible_peaks.sb)
C13=$(sb --array=0-6 $(dep "${G}") scripts/13_ctr_matched_calls.sb)
C14=$(sb --array=0-6 $(dep "${G}") scripts/14_ctr_relaxed.sb)
C15=$(SUBSET="D15:1,3,4 D18:1,3,4 D22:2,4,5" OUT_NAME=ctr_trio sb --array=0-2 --export=ALL $(dep "${P3}" "${P4}" "${C13}" "${C14}") scripts/15_ctr_multirep_sets.sb)
C18=$(sb --job-name=isl1_final --time=02:00:00 --output="${LOGS}/isl1_final.%j.out" $(dep "${C15}") --wrap "bash ${R}/scripts/18_isl1_ctr_final_peaks.sh")
H16a=$(MARK=H3K27ac sb --export=ALL $(dep "${P5}") scripts/16_histone_diff_counts.sb)
H16m=$(MARK=H3K27me3 MERGE_GAP=3000 sb --export=ALL $(dep "${P5}") scripts/16_histone_diff_counts.sb)
H17a=$(MARK=H3K27ac DROP_SAMPLE=D18_H3CTR_rep2_H3K27ac sb --export=ALL $(dep "${H16a}") scripts/17_histone_deseq2.sb)
H17m=$(MARK=H3K27me3 DROP_SAMPLE=D22_H3CTR_rep2_H3K27me3 sb --export=ALL $(dep "${H16m}") scripts/17_histone_deseq2.sb)
H19a=$(MARK=H3K27ac WIN=500 sb --export=ALL $(dep "${H16a}") scripts/19_histone_diff_sensitive.sb)
H19m=$(MARK=H3K27me3 WIN=5000 sb --export=ALL $(dep "${H16m}") scripts/19_histone_diff_sensitive.sb)
H21a=$(MARK=H3K27ac sb --export=ALL $(dep "${H17a}" "${H19a}") scripts/21_histone_diffbind.sb)
H21m=$(MARK=H3K27me3 DB_PEAK_DIR=macs3_broad DB_PEAK_SUFFIX=_peaks.filtered.broadPeak DB_PEAK_FORMAT=bed DB_SUMMITS=0 \
  sb --export=ALL $(dep "${H17m}" "${H19m}") scripts/21_histone_diffbind.sb)
T1=$(sb $(dep "${C18}" "${H16a}" "${H16m}") scripts/23a_tracks_prepare.sb)
T2=$(sb --array=0-69%20 $(dep "${T1}") scripts/23b_tracks_bigwig.sb)
T3=$(sb $(dep "${T2}" "${H21a}" "${H21m}") scripts/23c_tracks_finalize.sb)
X=$(sb $(dep "${T1}" "${H21a}" "${H21m}") scripts/24_tenx_before_after.sb)
cat <<EOF
01 align ${A} | 02 IgG ${G} | 03 ${P3} | 04 ${P4} | 05 ${P5}
13 ${C13} | 14 ${C14} | 15 ${C15} | 18 ${C18}
16 ${H16a} ${H16m} | 17 ${H17a} ${H17m} | 19 ${H19a} ${H19m} | 21 ${H21a} ${H21m}
23a ${T1} | 23b ${T2} | 23c ${T3} | 24 ${X}
EOF
