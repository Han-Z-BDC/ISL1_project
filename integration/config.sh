# integration/config.sh -- input and output paths for the ISL1 CUT&RUN x bulk RNA-seq integration.
# Sourced by every script in scripts/.

PROJECT=/pl/active/ZhuLab_BDC/GitHub/ISL1_project
INT_DIR=${PROJECT}/integration
RESULTS=${INT_DIR}/results
LOGS=${INT_DIR}/logs

# ---- CUT&RUN: ISL1, control line, D15 / D18 / D22 (cutnrun_ISL1_timecourse, on scratch) ------------
# Step 06 sites, step 07 IgG greylist and DiffBind (the 9 libraries, the same counts as the exploratory
# cutnrun_ISL1_KO_histone step 27).
CNR_RESULTS=/scratch/alpine/hzhu10@xsede.org/ISL1_cutnrun_timecourse/results
CNR_DIFF=${CNR_RESULTS}/diff
CNR_SITES=${CNR_DIFF}/sites_used.bed          # 401 bp sites after the IgG greylist
CNR_COUNTS=${CNR_DIFF}/counts_diffbind.tsv    # DiffBind reads per site, 9 libraries (steps 03, 05, 08, 09, 10)
CNR_NORM=${CNR_DIFF}/norm_factors.tsv         # size factors; the integration uses method diffbind, background
CNR_SAMPLES=${CNR_DIFF}/samples.tsv           # library -> day, batch
CNR_TABLES=${CNR_DIFF}/tables                 # per-site differential binding tables
CHROM_SIZES=${CNR_RESULTS}/genome.txt
CNR_OVERLAP=${CNR_RESULTS}/peak_overlap      # step 12: overlap classes of the per-day peak calls (step 09)

# ---- Bulk RNA-seq: Ctr vs ISL1 KO, D15 / D18 / D22 (run YJZ5LS) ---------------------------------
BULK_DE=${PROJECT}/bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2
GTF=/pl/active/ZhuLab_BDC/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf

# ---- Islet cell-type marker genes (step 06, report section 06)
MARKERS=${INT_DIR}/markers.tsv

# ---- Environments -------------------------------------------------------------------------------
PY=/projects/hzhu10@xsede.org/software/anaconda/envs/bulkRNAseq/bin/python
R_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/r_env
