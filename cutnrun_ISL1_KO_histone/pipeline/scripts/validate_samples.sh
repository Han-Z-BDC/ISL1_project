#!/bin/bash
# Sanity-checks samples.tsv and confirms every FASTQ it references actually
# exists before any pipeline step burns compute time on a typo.
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_DIR=$(dirname "${SCRIPT_DIR}")
source "${PIPELINE_DIR}/config.sh"

[[ -s "${SAMPLE_SHEET}" ]] || { echo "ERROR: Missing ${SAMPLE_SHEET}" >&2; exit 1; }

awk -F '\t' '
NR==1 {
  if ($1!="library_id"||$2!="sample_name"||$3!="role"||$4!="control_sample"||$5!="replicate")
    bad("invalid header")
  next
}
function bad(x){ print "ERROR: " x > "/dev/stderr"; err=1 }
{
  if (NF!=5) bad("row " NR " must have 5 fields")
  if ($2=="TBD"||$3=="TBD"||$4=="TBD") bad("unresolved TBD on row " NR)
  if ($3!="TF" && $3!="IgG") bad("invalid role on row " NR ": " $3)
  if ($3=="IgG" && $4!="NA") bad("IgG control must be NA on row " NR)
  if (lib[$1]++) bad("duplicate library " $1)
  if (samp[$2]++) bad("duplicate sample " $2)
  role[$2]=$3; ctl[$2]=$4
}
END {
  for (s in ctl)
    if (role[s]=="TF" && (!(ctl[s] in role) || role[ctl[s]]!="IgG"))
      bad("invalid IgG control for " s)
  exit err
}' "${SAMPLE_SHEET}"

while IFS=$'\t' read -r library sample role control replicate; do
  [[ "${library}" == library_id ]] && continue
  for mate in R1 R2; do
    f="${RAW_DIR}/${library}_${mate}_001.fastq.gz"
    [[ -s "${f}" ]] || { echo "ERROR: Missing ${f}" >&2; exit 1; }
  done
done < "${SAMPLE_SHEET}"

echo "Sample sheet and FASTQs validated."
