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
- Library: Plasmidsaurus bulk RNA-seq, a 3'-end counting protocol (reads come mostly from the last ~400 nt of each
  transcript), forward-stranded. Fastqs are delivered adapter-trimmed and not deduplicated; each read name ends with a
  14-nt UMI (`..._GATACTCTTACGGC`).
- Samples: `config.sh` (`SAMPLES`); metadata `sample_sheet.csv` (`sample, day, genotype, replicate`;
  genotype `Ctr` = Ctr-H3, `KO` = ISL1-KO).
- Outputs (not in git): `results/`, `qc/`, `logs/` in this folder on Alpine.

## Pipeline

`pipeline/` started as a copy of the processing scripts in `../bulkRNAseq/scripts/` (ISL1_project commit `290c3e1`)
and keeps their references (GRCh38 2020-A, `STARindex_2.7.5c`, `genes.gtf`) and conda env (`bulkRNAseq`, now with
`umi_tools`). Quantification differs, because these are 3'-end UMI libraries: STAR BAMs are UMI-deduplicated and
genes are counted with featureCounts (one count per molecule, no length correction) instead of RSEM.
All steps run with `sbatch` on acpu, in this order:

| Step | Script | Output |
|---|---|---|
| FastQC + MultiQC on raw reads | `pipeline/00_fastqc.sb` | `qc/raw_fastqc/` |
| fastp trimming (poly-G trimmed) | `pipeline/01_fastp.sb` | `results/trimmed/`, `qc/fastp/` |
| STAR alignment (array 0-11) | `pipeline/02_star_align.sb` | `results/star/<sample>/` |
| UMI deduplication, `umi_tools dedup` (array 0-11) | `pipeline/03_umi_dedup.sb` | `results/dedup/`, `qc/umi_dedup/` |
| Gene counts, `featureCounts -s 1` | `pipeline/03b_featurecounts.sb` | `results/counts/featureCounts.gene.txt` |
| BAM QC: flagstat, idxstats (STAR BAM); rseqc infer_experiment (dedup BAM) (array 0-11) | `pipeline/04_qc_bam.sb` | `qc/bam/` |
| MultiQC over all QC | `pipeline/04b_multiqc_full.sb` | `qc/multiqc_full_report.html` |
| CPM bigwigs from dedup BAMs (array 0-11) | `pipeline/05_bigwig.sb` | `results/bigwig/` |

UMI deduplication reads the UMI from the read name (`--extract-umi-method=read_id --umi-separator=_`, directional
method). featureCounts counts uniquely mapped reads on exons, summarized by `gene_id`, forward strand. The rseqc
`infer_experiment` output from `04_qc_bam.sb` must confirm forward strandedness before the counts are used.
There is no gene-body coverage step: in 3'-end libraries reads pile up at the 3' end by design, so it cannot show
RNA degradation.

## Differential expression (DESeq2)

`pipeline/06_deseq2_timecourse.sb` writes `results/counts/gene_annotation.tsv` (gene_id, gene_name, gene_type from the
GTF) and runs `pipeline/06_deseq2_timecourse.R` in `r_env` (DESeq2 1.46, apeglm 1.28). Output: `results/DE_DESeq2/`.

- Model: `~ genotype + day + genotype:day`, fitted as the equivalent 6-level `group` (Ctr_D15 ... KO_D22). Day is a
  factor. No replicate term: rep1 and rep2 are not differentiation batches. All 12 samples in one fit.
- Genes kept: >= 10 counts in >= 2 samples.
- Wald tests (adjusted p from the MLE fold change; `log2FC_shrunk` from apeglm for ranking and plots):
  `KO_vs_Ctr_D15/D18/D22.tsv`; day effects within each genotype, `<Ctr|KO>_D18_vs_D15`, `_D22_vs_D15`, `_D22_vs_D18`.
- Likelihood-ratio tests: `LRT_interaction.tsv` (full vs `~ genotype + day`: the KO effect changes with day) and
  `LRT_genotype.tsv` (full vs `~ day`: a KO effect at any day).
- Also: `summary.tsv` (genes at padj < 0.05 per test), `normalized_counts.tsv`, `vst.tsv`, `pca_top500.tsv`,
  `pca_variance.tsv`, `sample_correlation_vst.tsv`, `size_factors.tsv`, `sessionInfo.txt`.

limma-voom as a cross-check is to be added.
