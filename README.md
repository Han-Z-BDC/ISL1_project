# ISL1 rescue bulk RNA-seq

Bulk RNA-seq for the ISL1 rescue stem-cell islet differentiation experiment — same
design as the companion `cutnrun-pipeline` project (Ctr-GFP vs ISL1-KO-GFP vs
ISL1-KO-rescue, Day 12 / Day 15, n=2 replicates/condition; see that repo / lab notes for
the actual-harvest-day caveat on "D12"/"D15" labels).

12 single-end samples (`raw_data/ISL1_bulk_RNA/`, symlinked from
`/pl/active/ZhuLab_BDC/raw_data/ISL1_bulk_RNA`):

D12: Ctr-GFP (x2), KO-GFP (x2), KO-ISL1-rescue (x2)
D15: Ctr-GFP (x2), KO-GFP (x2), KO-ISL1-rescue (x2)

## Environment

Dedicated conda env, built via `environment.yml`:

```
micromamba create -p /projects/$USER/software/anaconda/envs/bulkRNAseq -f environment.yml
```

fastqc, multiqc, fastp, STAR 2.7.5c, bowtie2, samtools, subread (featureCounts), rseqc,
deeptools. Differential expression is done separately in the lab's existing `r_env`
(edgeR/limma), not duplicated here.

## Pipeline

1. `scripts/00_fastqc.sb` — FastQC + MultiQC on raw fastqs (`sbatch`, `acpu` partition)
2. `scripts/01_fastp.sb` — adapter/quality trim
3. `scripts/02_star_align.sb` — STAR alignment to `STARindex_2.7.5c` (see
   `docs/star_index_notes.md` for why this index and not the Cell Ranger ARC-bundled
   `star/` one), `--quantMode GeneCounts`
4. `scripts/03_featurecounts.sb` — gene-level counts via `featureCounts` from the STAR
   BAMs
5. `scripts/04_qc_bam.sb` — `samtools flagstat`/`idxstats`, `rseqc` (strandedness,
   gene-body coverage), second MultiQC pass over the full pipeline
6. `scripts/05_bigwig.sb` — `deeptools bamCoverage` for genome-browser tracks
7. Differential expression (edgeR, `r_env`) — separate script/notebook, not yet written

Steps 1 (FastQC) is light and runs directly via `sbatch` on `acpu`. STAR alignment (step
3) is the heavy step — also `sbatch`, sized to `acpu`'s 3.75GB/cpu cap unless a sample
needs more, in which case use `amem`.
