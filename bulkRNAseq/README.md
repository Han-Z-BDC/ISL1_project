# ISL1 rescue bulk RNA-seq

Bulk RNA-seq for the ISL1 rescue stem-cell islet differentiation experiment — same
design as the companion `cutnrun-pipeline` project (Ctr-GFP vs ISL1-KO-GFP vs
ISL1-KO-rescue, Day 12 / Day 15, n=2 replicates/condition; see that repo / lab notes for
the actual-harvest-day caveat on "D12"/"D15" labels).

12 single-end samples (`raw_data/ISL1_bulk_RNA/`, symlinked from
`/pl/active/ZhuLab_BDC/raw_data/ISL1_bulk_RNA`), metadata in `scripts/sample_sheet.csv`:

D12: Ctr-GFP (x2), KO-GFP (x2), KO-ISL1-rescue (x2)
D15: Ctr-GFP (x2), KO-GFP (x2), KO-ISL1-rescue (x2)

Design is **not** a full genotype x construct factorial — Ctr cells only ever get GFP,
only KO cells get the ISL1-rescue construct — so downstream DE uses one 6-level `group`
factor (day x group3), not separate main effects. See `scripts/06_*`/`07_*` for details.

## Environment

Two dedicated conda envs, per lab convention (one per pipeline, built from a checked-in
`environment.yml`/package list rather than reusing an unrelated env):

- **`bulkRNAseq`** (this repo's `environment.yml`): fastqc, multiqc, fastp,
  STAR 2.7.5c, bowtie2, samtools, RSEM, featureCounts (subread, cross-check only),
  rseqc, deeptools.
  ```
  micromamba create -p /projects/$USER/software/anaconda/envs/bulkRNAseq -f environment.yml
  ```
- **`r_env`** (lab's existing shared R env) — DESeq2 + apeglm added to it for this
  project, alongside the edgeR/limma already there. Differential expression is run
  through *both* engines and cross-checked (see step 8) rather than committing to one.
  Pinned package list: `envs/r_env.yml`
  ```
  micromamba create -p /projects/$USER/software/anaconda/envs/r_env -f envs/r_env.yml
  ```

## Pipeline

1. `scripts/00_fastqc.sb` — FastQC + MultiQC on raw fastqs (`sbatch`, `acpu` partition)
2. `scripts/01_fastp.sb` — adapter/quality trim
3. `scripts/02_star_align.sb` — STAR alignment to `STARindex_2.7.5c` (see
   `docs/star_index_notes.md` for why this index and not the Cell Ranger ARC-bundled
   `star/` one), ENCODE-standard STAR options, `--quantMode TranscriptomeSAM` for RSEM
4. `scripts/03_rsem_quant.sb` — gene/isoform-level quantification via RSEM (ENCODE's
   standard quantification tool downstream of STAR); featureCounts available as a
   quick sanity-check but RSEM output is what steps 6-8 consume
5. `scripts/04_qc_bam.sb` — `samtools flagstat`/`idxstats`, `rseqc` (strandedness,
   gene-body coverage), second MultiQC pass over the full pipeline
6. `scripts/05_bigwig.sb` — `deeptools bamCoverage` for genome-browser tracks
7. `scripts/06_differential_expression.R` — DE via **limma-voom** (`r_env`), reads
   `results/rsem/*.genes.results` → `results/DE/`
8. `scripts/07_differential_expression_DESeq2.R` — the same design/contrasts via
   **DESeq2** (`r_env`) → `results/DE_DESeq2/`
9. `scripts/08_compare_DE_methods.R` — merges both methods' results per contrast:
   logFC correlation, significant-gene overlap (Venn), confirmation rate →
   `results/DE_comparison/`

Steps 1 (FastQC) and 7-9 (R, on small count tables) are light and run directly via
`sbatch`/login node. STAR alignment (step 3) and RSEM (step 4) are the heavy steps —
`sbatch`, sized to `acpu`'s 3.75GB/cpu cap unless a sample needs more, in which case use
`amem`. Any conda env *build/install* step (not analysis) goes through an `acompile`
node, never the plain login node — package linking/installation is CPU-heavy enough to
get silently SIGKILLed by the login node's watchdog otherwise.

Steps 1-2 are written; 2-5 (STAR/RSEM/QC/bigwig `sbatch` scripts) are still to be
written once FastQC results are in.
