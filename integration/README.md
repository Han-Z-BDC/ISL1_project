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
  other sites, contrasts D18 vs D15, D22 vs D18 and D22 vs D15. Figure: cumulative distribution of the site log2FC per group
  against all sites.
- **Candidate time-dependent direct target.** A gene with a site in its domain whose binding change matches the
  timing of its KO effect (the sign of the KO effect is not used): cluster 1 (strongest at D18) a site gained D18 vs
  D15; clusters 3 and 4 (emerging from D18) a site gained D18 or D22 vs D15; cluster 5 (late) a site gained D22 vs
  D15; clusters 2 and 6 (fading) a site lost. Each candidate also records whether its own KO effect is significant
  (padj < 0.05) on the day of the matching change. Gene-level check: genes with at least one gained site per group
  against genes without a KO effect, logistic regression with domain length and the number of sites in the domain.

## Steps

Run on Alpine from this folder with `sbatch`; outputs in `results/` (not in git).

| Step | Script | Output |
|---|---|---|
| 01 GREAT domains, site -> gene | `scripts/01_great_site_to_gene.sb` (acpu, 2 CPUs, 7 GB) | `results/01_great/`: `gene_regions.bed`, `site_gene.tsv`, `summary.tsv` |
| 02 Direct targets per KO-effect group | `scripts/02_direct_targets.sb` (acpu, 2 CPUs, 7 GB) | `results/02_direct_targets/`: `gene_table.tsv` (group, KO log2FC per day, sites per gene), `enrichment.tsv` |

| 03 Binding change per KO-effect group | `scripts/03_binding_change.sb` (acpu, 2 CPUs, 7 GB) | `results/03_binding_change/`: `camera.tsv`, `site_stats.tsv` (per site: DiffBind and voom log2FC, FDR, gained), `site_groups.tsv` |
| 04 Candidate time-dependent direct targets | `scripts/04_candidates.sb` (acpu, 2 CPUs, 7 GB) | `results/04_candidates/`: `candidate_genes.tsv`, `candidates.tsv` (every matching site), `summary.tsv`, `gained_enrichment.tsv` |
| 05 Binding signal per day at each cluster's sites | `scripts/05_binding_signal.sb` (acpu, 2 CPUs, 7 GB) | `results/05_binding_signal/`: `site_signal.tsv`, `group_signal.tsv`, `wilcoxon.tsv` |

Step 02 reads `timecourse_clusters_LRT_interaction.tsv`, which the bulk report builder writes on the Mac; it is
copied to `bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2/` on Alpine.

Report (Mac, from this folder): `bash scripts/fetch_results.sh`, then `~/anaconda3/bin/python3 build_report.py`
-> `results/report_integration.html` (steps 01-05: tables, bound share per group, cumulative distributions of the
site log2FC per group, candidate lists, boxplots of the binding signal per day and cluster; section 06, built in the
report from the step 01 and 05 tables and the bulk normalized counts: binding and expression of islet cell-type
marker genes).

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

### 03 Binding change per KO-effect group

758 of the 14,431 sites are gained (FDR < 0.05 up in D18 vs D15 or D22 vs D15) and 12 lost. Site log2FC from voom
and from DiffBind agree (Pearson r 0.97 for D18 vs D15 and D22 vs D15, 0.96 for D22 vs D18). A site counts for every
group with a linked gene. Cells: median DiffBind log2FC of the group's sites (other sites), camera direction and FDR
(limma's default inter-site correlation).

| Group | Sites | D18 vs D15 | D22 vs D18 | D22 vs D15 | % gained (other sites) |
|---|---|---|---|---|---|
| No KO effect | 7,631 | 0.317 (0.345), down 0.88 | -0.006 (-0.060), **up 0.023** | 0.289 (0.229), up 0.39 | 4.1 (6.6) |
| Cluster 1 | 850 | 0.330 (0.330), down 0.89 | -0.043 (-0.032), down 0.89 | 0.235 (0.265), down 0.74 | 4.7 (5.3) |
| Cluster 2 | 629 | 0.343 (0.329), up 0.89 | -0.060 (-0.032), down 0.82 | 0.241 (0.264), down 0.94 | 7.2 (5.2) |
| Cluster 3 | 390 | 0.307 (0.330), down 0.89 | -0.038 (-0.033), up 0.94 | 0.253 (0.263), up 0.95 | 5.1 (5.3) |
| **Cluster 4** | 613 | **0.388** (0.328), up 0.23 | **-0.084** (-0.030), **down 0.023** | 0.224 (0.264), down 0.74 | **9.8** (5.1) |
| Cluster 5 | 219 | 0.301 (0.330), down 0.88 | -0.021 (-0.033), up 1.00 | 0.241 (0.263), down 0.74 | 2.7 (5.3) |
| Cluster 6 | 218 | 0.324 (0.330), up 0.89 | -0.054 (-0.032), down 0.25 | 0.213 (0.264), down 0.39 | 3.7 (5.3) |
| Stays | 3,077 | 0.324 (0.331), down 0.89 | -0.027 (-0.034), up 1.00 | 0.258 (0.264), down 0.95 | 4.9 (5.3) |
| Other KO | 625 | 0.332 (0.329), down 0.89 | -0.020 (-0.033), up 0.89 | 0.268 (0.263), up 0.95 | 4.6 (5.3) |

- The whole site set rises from D15 to D18 (median log2FC about 0.33) and stays there (D22 vs D18 about -0.03); the
  groups follow that rise.
- **Cluster 4** (KO effect down, emerging and growing) rises most from D15 to D18 (FDR 0.23) and is the one cluster
  whose sites then **fall from D18 to D22 relative to all other sites** (median -0.084 vs -0.030, FDR 0.023). Its
  binding peaks at D18, when its KO effect first appears. The genes without a KO effect move the opposite way
  (FDR 0.023), as the complement of the KO-effect clusters.
- Sites move together: the inter-site correlation estimated by camera is 0.21-0.31, far above limma's default 0.01.
  With the estimated correlation every FDR is >= 0.97 (D22 vs D18: 0.999), so the cluster 4 results are leads, not
  firm differences.
- Sites near genes without a KO effect are gained less often (4.1%) than the other sites (6.6%), which lie mostly
  near genes with a KO effect.

### 04 Candidate time-dependent direct targets

| Group | Genes | With a gained site (%) | Gained site vs no KO effect, adjusted OR (FDR) | Candidates | KO effect significant on the day of the change |
|---|---|---|---|---|---|
| No KO effect | 10,605 | 390 (3.7) | 1 | | |
| Cluster 1: down, strongest at D18 | 759 | 38 (5.0) | 0.95 (0.78) | 35 | 32 |
| Cluster 2: up early, down by D22 | 383 | 41 (10.7) | 1.39 (0.25) | 1 (lost site) | 0 |
| Cluster 3: up, emerging from D18 | 358 | 19 (5.3) | 1.07 (0.78) | 19 | 16 |
| **Cluster 4: down, emerging and growing** | 357 | **54 (15.1)** | **2.06 (3.5e-4)** | **54** | **43** |
| Cluster 5: up, late (D22) | 193 | 6 (3.1) | 0.60 (0.45) | 2 | 0 |
| Cluster 6: up early, down by D22 | 175 | 6 (3.4) | 0.47 (0.25) | 0 | 0 |
| Stays | 3,147 | 155 (4.9) | 1.12 (0.45) | | |
| Other KO | 569 | 24 (4.2) | 0.86 (0.66) | | |

- **111 candidate genes**, 91 of them with their own KO effect significant on the day of the matching binding change.
- **Cluster 4 is the one group linked to gained binding at the gene level**: 15.1% of its genes have a gained site,
  against 3.7% of genes without a KO effect (adjusted OR 2.06, FDR 3.5e-4, after domain length and number of sites).
  Its genes lose expression in the KO from D18 on, so ISL1 binding gained after D15 goes with genes ISL1 comes to
  activate. All 54 of its genes with a gained site are candidates; for 53 the best site is gained in D18 vs D15
  (17 of them also in D22 vs D15).
- Cluster 2 genes have gained sites more often (10.7%) but not after adjusting for their many sites, and binding
  near them is almost never lost (1 candidate), so the fading KO effect is not matched by lost binding.
- Islet-related genes among the candidates include the transcription factors HHEX, MEIS2, TCF7L2 and REST
  (cluster 1), NKX6-1 (cluster 3) and ETV1 (cluster 4), and the endocrine genes PCP4 and SYT4 (cluster 4).

### 05 Binding signal per day at each cluster's sites

Signal: reads per site / step 27 background size factor, mean of the 3 libraries per day. Paired Wilcoxon between
days within each group (BH).

| Group | Sites | Median signal D15 / D18 / D22 | Sites up D18 vs D15 | Sites up D22 vs D18 | Median log2 ratio D22 vs D18 (FDR) |
|---|---|---|---|---|---|
| No KO effect | 7,631 | 24.0 / 32.3 / 32.5 | 84.7% | 50.4% | 0.00 (3.3e-6) |
| Cluster 1 | 850 | 23.4 / 31.3 / 29.2 | 82.8% | 39.9% | -0.15 (9.0e-15) |
| Cluster 2 | 629 | 22.3 / 31.4 / 28.7 | 84.4% | 36.1% | -0.21 (7.0e-19) |
| Cluster 3 | 390 | 22.9 / 30.7 / 29.8 | 84.9% | 43.3% | -0.10 (1.3e-4) |
| **Cluster 4** | 613 | 22.9 / **33.7** / 27.4 | 86.6% | **27.9%** | **-0.31** (7.2e-38) |
| Cluster 5 | 219 | 23.4 / 31.9 / 30.9 | 85.4% | 44.3% | -0.06 (0.010) |
| Cluster 6 | 218 | 23.0 / 31.2 / 28.3 | 86.7% | 37.2% | -0.17 (5.5e-8) |

- Binding rises from D15 to D18 at 83-87% of the sites of every group (all FDR < 1e-12), the overall ISL1 increase.
- **From D18 to D22 the sites near the KO-effect clusters lose signal while the reference stays flat** (50.4% of
  sites up, median ratio 0.00). Cluster 4 shows both the largest rise to D18 (median 22.9 -> 33.7) and the largest
  fall to D22 (72% of its sites down). Per site, none of these D22 vs D18 changes reaches significance in step 27.
  Against all other sites, cluster 4's fall is the one that stands out (step 03 camera D22 vs D18, FDR 0.023, not
  significant with the estimated inter-site correlation). The batch-aware model gives a smaller shift than these
  day means (median DiffBind log2FC -0.084 vs -0.030 for the other sites).
- The paired tests use sites as units; with hundreds of sites per group even small shifts give tiny p values
  (the reference's D22 vs D18 FDR of 3.3e-6 comes with a median change of 0). They describe each group; whether a
  group differs from the others is the camera test of step 03.

### 06 Islet cell-type programs (report section)

Marker genes: alpha GCG, ARX, IRX1, IRX2, TTR, GC, LOXL4, DPP4, FAP, TM4SF4, CRYBA2, PTPRT, KCTD12, MAFB, POU3F4,
PLCE1, TMEM236, SMIM24, SLC7A2; beta INS, IAPP, MAFA, NKX6-1, PDX1, SLC30A8, G6PC2, HADH, PCSK1; delta SST, HHEX,
LEPR, RBP4. Binding: signal of the sites linked to them (step 05 signal); expression: mean normalized counts in
control islets.

| Sites near | Sites | Median signal D15 / D18 / D22 | Sites up D18 vs D15 | Sites up D22 vs D18 (FDR) |
|---|---|---|---|---|
| All sites | 14,431 | 23.1 / 31.7 / 29.9 | 84.9% | 43.5% (1e-126) |
| Alpha genes | 26 | 17.7 / 29.6 / 24.9 | 96.2% | 38.5% (0.039) |
| Beta genes | 16 | 19.1 / 39.7 / 35.0 | 100% | 18.8% (0.035) |
| Delta genes | 9 | 23.0 / 33.2 / 26.5 | 88.9% | 0% (0.006) |

- Cluster 4 holds 9 of the 19 alpha genes (GCG, ARX, IRX2, TTR, DPP4, LOXL4, TMEM236, PLCE1, SLC7A2), all lost in
  the KO with a growing effect (GCG -6.2 / -8.9 / -11.2, ARX -1.4 / -3.4 / -5.3), together with SST and the beta
  genes PDX1 and HADH.
- In control islets these genes keep rising from D18 to D22 (GCG 35,545 -> 203,556; ARX 485 -> 1,020), and so does
  ISL1 (373 / 1,990 / 4,425 at D15 / D18 / D22).
- Binding near the marker genes of all three cell types peaks at D18 and falls by D22, more than at all sites,
  while their expression keeps rising. The rise-then-fall is shared by the alpha, beta and delta programs, not
  specific to alpha. Among alpha genes it is carried by LOXL4, DPP4, TTR, GC and MAFB; binding near GCG keeps rising
  (48 / 77 / 84) and near ARX stays level after D18 (32 / 42 / 41).
- Bound alpha genes without a significant site (GCG, ARX, TTR) are not step 04 candidates.
- The marker sets are small (9-26 sites) and the paired tests treat sites as independent, so these are descriptive.
