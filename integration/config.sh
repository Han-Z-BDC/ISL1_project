# integration/config.sh -- input and output paths for the ISL1 CUT&RUN x bulk RNA-seq integration.
# Sourced by every script in scripts/.

PROJECT=/pl/active/ZhuLab_BDC/GitHub/ISL1_project
INT_DIR=${PROJECT}/integration
RESULTS=${INT_DIR}/results
LOGS=${INT_DIR}/logs

# ---- CUT&RUN: ISL1, control line, D15 / D18 / D22 -----------------------------------------------
# Currently from cutnrun_ISL1_KO_histone (step 26 sites, step 27 DiffBind, on scratch). Switch this
# block to cutnrun_ISL1_timecourse once the unified CUT&RUN data are there.
CNR_DIFF=/scratch/alpine/hzhu10@xsede.org/ISL1_KO_histone_cutNrun/results/isl1_diff
CNR_SITES=${CNR_DIFF}/sites_used.bed          # 401 bp sites after the IgG greylist
CNR_COUNTS=${CNR_DIFF}/counts_sites.tsv       # reads per site, 9 libraries
CNR_NORM=${CNR_DIFF}/norm_factors.tsv         # DiffBind / DESeq2 normalization factors
CNR_SAMPLES=${CNR_DIFF}/samples.tsv           # library -> day, batch
CNR_TABLES=${CNR_DIFF}/tables                 # per-site differential binding tables
CHROM_SIZES=/scratch/alpine/hzhu10@xsede.org/ISL1_KO_histone_cutNrun/results/isl1_consensus/genome.txt
# DiffBind counts and background size factors of the 9 libraries, exported by step 00 (used by steps 03, 05, 08)
DB_COUNTS=${RESULTS}/00_diffbind_counts

# ---- Bulk RNA-seq: Ctr vs ISL1 KO, D15 / D18 / D22 (run YJZ5LS) ---------------------------------
BULK_DE=${PROJECT}/bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2
GTF=/pl/active/ZhuLab_BDC/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf

# ---- Islet cell-type marker genes (step 06, report section 06)
MARKERS=${INT_DIR}/markers.tsv

# ---- Environments -------------------------------------------------------------------------------
PY=/projects/hzhu10@xsede.org/software/anaconda/envs/bulkRNAseq/bin/python
R_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/r_env
DIFFBIND_ENV=/projects/hzhu10@xsede.org/software/anaconda/envs/diffbind
