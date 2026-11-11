# ISL1 KO time-course bulk RNA-seq (run YJZ5LS)

Bulk RNA-seq of stem-cell islet differentiation, control (Ctr-H3) vs ISL1 KO, at D15, D18 and D22,
2 replicates each: a full 2 genotype x 3 day design, 12 samples. The run also has D9 samples
(YJZ5LS 1, 2, 9, 10); they are not part of this analysis.

| Day | Ctr | KO |
|---|---|---|
| D15 | YJZ5LS_3, _4 | YJZ5LS_11, _12 |
| D18 | YJZ5LS_5, _6 | YJZ5LS_13, _14 |
| D22 | YJZ5LS_7, _8 | YJZ5LS_15, _16 |

- Raw fastqs: `/pl/active/ZhuLab_BDC/raw_data/ISL1_bulk_RNA/YJZ5LS_*.fastq.gz` (single-end, ~94 bp).
- Samples: `config.sh` (`SAMPLES`); metadata `sample_sheet.csv` (`sample, day, genotype, replicate`;
  genotype `Ctr` = Ctr-H3, `KO` = ISL1-KO).
- Outputs (not in git): `results/`, `qc/`, `logs/` in this folder on Alpine.

## Pipeline

`pipeline/` is a copy of the processing scripts in `../bulkRNAseq/scripts/` (ISL1_project commit `290c3e1`),
changed only to use this folder's `config.sh`, output and log paths. Same references (GRCh38 2020-A,
`STARindex_2.7.5c`, RSEM `GRCh38_2020A`) and conda env (`bulkRNAseq`). All steps run with `sbatch` on acpu.

| Step | Script | Output |
|---|---|---|
| FastQC + MultiQC on raw reads | `pipeline/00_fastqc.sb` | `qc/raw_fastqc/` |
| fastp trimming (poly-G trimmed) | `pipeline/01_fastp.sb` | `results/trimmed/`, `qc/fastp/` |
| STAR alignment (array 0-11) | `pipeline/02_star_align.sb` | `results/star/<sample>/` |
| RSEM quantification (array 0-11) | `pipeline/03_rsem_quant.sb` | `results/rsem/` |
| BAM QC: flagstat, idxstats, rseqc (array 0-11) | `pipeline/04_qc_bam.sb` | `qc/bam/` |
| MultiQC over all QC | `pipeline/04b_multiqc_full.sb` | `qc/multiqc_full_report.html` |
| CPM bigwigs (array 0-11) | `pipeline/05_bigwig.sb` | `results/bigwig/` |

`pipeline/rsem_prepare_reference.sb` builds the shared RSEM reference; it already exists and the script skips it.

`03_rsem_quant.sb` uses `--strandedness forward`, which was confirmed on the earlier run (C7S59W). For this
run, the rseqc `infer_experiment` output from `04_qc_bam.sb` must confirm it before the RSEM counts are used.

Differential expression (DESeq2, with limma-voom as cross-check) is to be added once counts exist.
