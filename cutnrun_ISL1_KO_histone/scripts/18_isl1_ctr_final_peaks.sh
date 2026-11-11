#!/bin/bash
# Step 18: the kept control ISL1 peak sets (agreed 2026-11-09), from the three-replicate run of
# step 15 (results/ctr_trio: D15/D18 rep1 + rep3 + rep4, D22 rep2 + rep4 + rep5):
#   macs3_ge2     MACS3 q 0.01 consensus, region with peaks from >= 2 of the 3 replicates
#   idr_any_pair  pooled relaxed peaks passing IDR < IDR_THRESHOLD in at least one replicate pair
# Columns (BED, no header): chr start end name n_reps reps score n_batch3 n_earlier
#   score = max -log10 p of the supporting calls (macs3_ge2: per-replicate q 0.01 calls;
#   idr_any_pair: pooled relaxed call)
# Written to results/isl1_ctr_final/<set>/<day>_Ctr_ISL1.bed on scratch and copied to the same
# path under the repo on /pl/active (results/ is not tracked; /pl/active is not purged).
# Light (awk + cp): runs on a login node.
#   bash scripts/18_isl1_ctr_final_peaks.sh

set -euo pipefail
IGG_METHOD=matched
source /pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun/study_config/config.sh
SRC=${OUT_DIR}/ctr_trio
F=${OUT_DIR}/isl1_ctr_final
KEEP=${STUDY_DIR}/results/isl1_ctr_final

summary=${F}/summary.tsv
mkdir -p "${F}"
echo -e "set\tday\treplicates\tpeaks\tMb\tmedian_bp" > "${summary}"
for set in macs3_ge2 idr_any_pair; do
  [[ ${set} == macs3_ge2 ]] && sub=A || sub=idr
  mkdir -p "${F}/${set}"
  for day in D15 D18 D22; do
    in=${SRC}/${day}/${sub}/k2.bed
    [[ -s "${in}" ]] || { echo "ERROR: missing ${in} -- run step 15 with SUBSET first" >&2; exit 1; }
    out=${F}/${set}/${day}_Ctr_ISL1.bed
    sort -k1,1V -k2,2n "${in}" | awk -v p="${day}_Ctr_ISL1_${set}_" 'BEGIN{FS=OFS="\t"} {$4 = sprintf("%s%05d", p, NR); print}' > "${out}"
    reps=$(cut -f6 "${out}" | tr ',' '\n' | sort -nu | paste -sd, -)
    echo -e "${set}\t${day}\t${reps}\t$(wc -l < "${out}" | tr -d ' ')\t$(awk '{t+=$3-$2} END{printf "%.2f", t/1e6}' "${out}")\t$(awk '{print $3-$2}' "${out}" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')" >> "${summary}"
  done
done
mkdir -p "${KEEP}"
cp -r "${F}/." "${KEEP}/"
column -t "${summary}"
echo "copied to ${KEEP}"
