#!/bin/bash
# Submit the whole ISL1 time-course analysis as one Slurm dependency chain (login node on Alpine,
# after bash study_config/setup_links.sh). Every job runs on amem (12 CPUs, 240 GB).
#   01 10x filter (batch 3 fastqs)  -> 02 align (pipeline 01_align.sb, 18 libraries)
#   -> 03 deduplicated matched IgG (pipeline 02)  -> 04 MACS3 per library  -> 05 consensus per day
#   -> 06 401 bp sites  -> 07 DiffBind + DESeq2 check  -> 08a-c tracks  -> 09 report + keepers
#   bash scripts/run_all.sh
set -euo pipefail
R=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd); cd "${R}"
source study_config/config.sh
[[ -L pipeline/config.sh ]] || { echo "ERROR: run study_config/setup_links.sh first" >&2; exit 1; }
dep () { local d=""; for j in "$@"; do [[ -n "${j}" ]] && d+=":${j}"; done; [[ -n "${d}" ]] && echo "--dependency=afterok${d}"; }
sb () { sbatch --parsable "${AMEM[@]}" "$@"; }

n_tenx=$(awk -F'\t' -v d="${TENX_DELIVERY}" 'NR>1 && $7==d' "${METADATA}" | wc -l)
n_lib=$(($(wc -l < "${SAMPLE_SHEET}") - 1))
n_tf=$(awk -F'\t' 'NR>1 && $3=="TF"' "${SAMPLE_SHEET}" | wc -l)

S01=$(sb --array=0-$((n_tenx - 1)) scripts/01_filter_10x.sb)
cd pipeline   # the pipeline finds config.sh and logs/ through SLURM_SUBMIT_DIR
S02=$(sb --array=0-$((n_lib - 1)) --time=12:00:00 $(dep "${S01}") scripts/01_align.sb)
S03=$(sb --time=12:00:00 $(dep "${S02}") scripts/02_build_igg_background.sb)
cd ..
S04=$(sb --array=0-$((n_tf - 1)) $(dep "${S03}") scripts/04_call_peaks.sb)
S05=$(sb $(dep "${S04}") scripts/05_consensus.sb)
S06=$(sb $(dep "${S05}") scripts/06_diff_sites.sb)
S07=$(sb $(dep "${S06}") scripts/07_diffbind.sb)
S08a=$(sb $(dep "${S07}") scripts/08a_tracks_prepare.sb)
S08b=$(sb --array=0-$((2 * n_tf - 1)) $(dep "${S08a}") scripts/08b_tracks_bigwig.sb)
S08c=$(sb $(dep "${S08b}") scripts/08c_tracks_finalize.sb)
S09=$(sb $(dep "${S07}" "${S08c}") scripts/09_report.sb)
cat <<EOF
01 filter ${S01} | 02 align ${S02} | 03 IgG ${S03} | 04 peaks ${S04} | 05 consensus ${S05}
06 sites ${S06} | 07 diffbind ${S07} | 08a ${S08a} | 08b ${S08b} | 08c ${S08c} | 09 report ${S09}
EOF
