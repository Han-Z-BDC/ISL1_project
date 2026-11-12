# ISL1 CUT&RUN x bulk RNA-seq integration

**Question.** Which genes does ISL1 regulate directly, and how does that regulation change over islet
differentiation (D15, D18, D22)? The analysis is RNA-led: the time-dependent KO effect from bulk RNA-seq defines the
gene groups, and ISL1 CUT&RUN tells which genes in each group are bound and whether binding near them changes over
time.

## Inputs

All paths are in `config.sh`.

| Data | Source | Used |
|---|---|---|
| Bulk RNA-seq, Ctr vs ISL1 KO, D15 / D18 / D22, 2 replicates (run YJZ5LS) | `bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2/` | `LRT_interaction.tsv` (KO effect changes with day), `LRT_genotype.tsv` (KO effect at any day), `KO_vs_Ctr_D15/D18/D22.tsv`, `timecourse_clusters_LRT_interaction.tsv` |
| ISL1 CUT&RUN, control line, D15 / D18 / D22, 3 libraries per day | `cutnrun_ISL1_KO_histone` step 26 (sites) and step 27 (DiffBind), on scratch | 14,431 sites of 401 bp, site counts, normalization factors, per-site differential tables |
| Genes | 10x GRCh38 2020-A `genes.gtf`, the annotation of the bulk counts | gene IDs, TSSs |

The CUT&RUN block of `config.sh` moves to `cutnrun_ISL1_timecourse/` once the unified CUT&RUN data are there.

## Definitions

- **KO-effect groups (bulk).** Changes over time: the 6 k-means clusters of the 2,231 genes with
  `LRT_interaction` padj < 0.05, from the bulk report. Stays: `LRT_genotype` padj < 0.05, `LRT_interaction`
  padj >= 0.05, and the same sign of the shrunk KO log2FC at D15, D18 and D22. Other KO: `LRT_genotype`
  padj < 0.05 with a sign change but no significant interaction. No KO effect: `LRT_genotype` padj >= 0.05 (the
  reference group). A padj of NA (964 genes in `LRT_interaction`) counts as not significant.
- **Site -> gene.** GREAT basal plus extension, default settings: basal domain 5 kb upstream to 1 kb downstream of
  the TSS, extended to the nearest other basal domain, at most 1 Mb. Domains are built from every gene in the GTF on
  chr1-22/X/Y; tests use the 16,563 genes tested in the bulk analysis. A site belongs to every gene whose domain
  contains its centre. Check: nearest TSS within 50 kb.
- **Direct target.** A gene with at least one ISL1 site (any day) in its domain. Enrichment per KO-effect group
  against the no-KO-effect genes: Fisher's exact test, and logistic regression `bound ~ group +
  log10(domain length)`, since longer domains collect more sites by chance (BH over groups).
- **Binding classes (step 27, background normalization).** Gained: FDR < 0.05 up in D18 vs D15 or D22 vs D15.
  Stable-bound: every other site. The sites lost over time (about 10) are listed, not tested.
- **Binding change per KO-effect group.** Sites are members of a set and the 9 libraries are the replicates: limma
  voom on the site counts (step 27 background normalization, `~ batch + day`), `camera` for each group against all
  other sites, contrasts D18 vs D15 and D22 vs D15. Figure: cumulative distribution of the site log2FC per group
  against all sites.

## Steps

Run on Alpine from this folder with `sbatch`; outputs in `results/` (not in git).

| Step | Script | Output |
|---|---|---|
| 01 GREAT domains, site -> gene | `scripts/01_great_site_to_gene.sb` (acpu, 2 CPUs, 7 GB) | `results/01_great/`: `gene_regions.bed`, `site_gene.tsv`, `summary.tsv` |
| 02 Direct targets per KO-effect group | `scripts/02_direct_targets.sb` (acpu, 2 CPUs, 7 GB) | `results/02_direct_targets/`: `gene_table.tsv` (group, KO log2FC per day, sites per gene), `enrichment.tsv` |

Step 02 reads `timecourse_clusters_LRT_interaction.tsv`, which the bulk report builder writes on the Mac; it is
copied to `bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2/` on Alpine.

## Results

### 01 GREAT domains, site -> gene

- 36,559 genes in the GTF; 16,546 of the 16,563 bulk-tested genes have a domain (the other 17 are 13
  mitochondrial genes and 4 on unplaced contigs). Domains average 180 kb; none is over 2 Mb.
- 14,430 of the 14,431 sites link to at least one gene: 3,580 to one gene, 7,639 to two, 3,211 to three or more
  (basal domains overlap). 29,103 site-gene links, 8,790 of them in a basal domain; median distance to the TSS
  23.4 kb.
- **10,056 of the tested genes (60.8%) have at least one site**, a median of 1. Bound-gene fractions of the
  KO-effect groups are therefore compared with this background, not read on their own.

### 02 Direct targets per KO-effect group

Bound = at least one site in the GREAT domain. Adjusted OR and FDR from the logistic regression with domain
length; all groups compared with no KO effect.

| Group | KO effect (median shrunk log2FC D15 / D18 / D22) | Genes | Bound | % bound | Adjusted OR | FDR | % bound, nearest TSS <= 50 kb (FDR) |
|---|---|---|---|---|---|---|---|
| No KO effect | | 10,605 | 6,319 | 59.6 | 1 | | 38.2 |
| Cluster 1 | down, strongest at D18 (-0.24 / -1.02 / -0.52) | 759 | 463 | 61.0 | 0.91 | 0.40 | 38.3 (0.79) |
| **Cluster 2** | up early, down by D22 (+0.37 / +0.25 / -0.34) | 383 | 286 | **74.7** | **1.72** | **7.9e-5** | **48.8 (0.0027)** |
| Cluster 3 | up, emerging from D18 (-0.03 / +0.53 / +0.72) | 358 | 223 | 62.3 | 0.96 | 0.74 | 39.1 (0.91) |
| Cluster 4 | down, emerging and growing (-0.12 / -0.92 / -2.07) | 357 | 255 | 71.4 | 1.24 | 0.33 | 44.8 (0.26) |
| Cluster 5 | up, late (-0.16 / -0.23 / +0.48) | 193 | 130 | 67.4 | 1.22 | 0.40 | 46.6 (0.15) |
| Cluster 6 | like cluster 2 (+0.38 / +0.46 / -0.33) | 175 | 113 | 64.6 | 1.08 | 0.73 | 40.0 (0.91) |
| Stays | one sign at all days | 3,147 | 1,904 | 60.5 | 0.95 | 0.40 | 40.4 (0.26) |
| Other KO | sign change, no interaction | 569 | 363 | 63.8 | 1.09 | 0.47 | 40.2 (0.79) |

- Only **cluster 2** is enriched for ISL1 binding after the domain-length adjustment, under both assignments.
- **Cluster 4** is bound more often (71.4%; Fisher OR 1.70, p 6e-6), but its genes have the longest domains
  (median 200 kb); with domain length in the model the excess is not significant.
- The other clusters and the genes with a lasting KO effect are bound no more often than genes without a KO effect.
  With 60% of genes bound, presence of a site alone separates the groups little.
