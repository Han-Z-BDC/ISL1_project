#!/bin/bash
#SBATCH --job-name=int12_abc_bedpe
#SBATCH --account=amc-general
#SBATCH --partition=acpu
#SBATCH --qos=cpu-testing
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/pl/active/ZhuLab_BDC/GitHub/ISL1_project/integration/logs/%x.%j.out
# Step 12: ABC enhancer -> gene links (bedpe) of the TFs in the TF summary (user, 2026-11-25: "subset the bedpe only
# containing enhancers linking to those genes"; groups ENP2, ENP3, early_SC_alpha = the ABC groups of CUT&RUN step 13).
# Input: diff_multiome ABC <group>/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.bedpe.gz
#   (enhancer chr/start/end, TSS chr/pos/pos, GENE|enhancer, ABC score, ., .).
# Kept: links whose gene (name before "|") is one of the 15 TFs; the gene's own promoter (element containing its TSS)
#   dropped, so only enhancer -> gene links remain.
# Output in <step 12 out dir>/abc_bedpe_TFs/: <group>_TF_links.bedpe.gz (same columns) and all_groups_TF_links.bedpe.gz
#   (the same plus the group in column 11).
# Submit from integration/:  sbatch scripts/12_abc_bedpe_TFs.sh
set -euo pipefail
source /pl/active/ZhuLab_BDC/GitHub/ISL1_project/integration/config.sh
ABC_DIR=/pl/active/ZhuLab_BDC/Han_analysis/diff_multiome_2026/out/abc/results   # as CUT&RUN step 13
ABC_GROUPS="ENP2 ENP3 early_SC_alpha"
TFS="SIM1,PDX1,PBX3,FOXO1,TBX3,PAX6,ETV1,ARX,HHEX,CREB3L2,FEV,ASCL2,LMX1A,MNX1,HEY2"
O=${RESULTS}/12_trajectory_binding/abc_bedpe_TFs; mkdir -p "${O}"
: > "${O}/all_groups_TF_links.bedpe"
for g in ${ABC_GROUPS}; do
  zcat "${ABC_DIR}/${g}/Predictions/EnhancerPredictionsFull_threshold0.02_self_promoter.bedpe.gz" |
    awk -F'\t' -v OFS='\t' -v tfs="${TFS}" 'BEGIN { n = split(tfs, a, ","); for (i = 1; i <= n; i++) k[a[i]] = 1 }
      { split($7, p, "|") } (p[1] in k) && !($1 == $4 && $2 <= $5 && $5 <= $3)' |
    sort -k1,1 -k2,2n -k3,3n -k7,7 > "${O}/${g}_TF_links.bedpe"
  awk -v OFS='\t' -v g="${g}" '{ print $0, g }' "${O}/${g}_TF_links.bedpe" >> "${O}/all_groups_TF_links.bedpe"
done
sort -k1,1 -k2,2n -k3,3n -k7,7 -k11,11 -o "${O}/all_groups_TF_links.bedpe" "${O}/all_groups_TF_links.bedpe"
echo "links per gene and group:"
awk -F'\t' '{ split($7, p, "|"); n[p[1] "\t" $11]++ } END { for (x in n) print x "\t" n[x] }' \
  "${O}/all_groups_TF_links.bedpe" | sort -k1,1 -k2,2
for f in "${O}"/*.bedpe; do gzip -f "${f}"; done
ls -la "${O}"
