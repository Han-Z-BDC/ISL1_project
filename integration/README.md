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
| ISL1 CUT&RUN, control line, D15 / D18 / D22, 3 libraries per day | `cutnrun_ISL1_KO_histone` (git tag `cutnrun-exploratory`) step 26 (sites) and step 27 (DiffBind), on scratch | 14,431 sites of 401 bp; per-site differential tables (main analysis: DiffBind, background normalization); the DiffBind object, from which step 00 exports counts and size factors |
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
- **Counts and signal.** Every step uses step 27's DiffBind counts (paired-end fragments per site, duplicates kept)
  and its background size factors (DESeq2 on 15 kb genome bins), exported by step 00. Signal per site and library =
  reads / size factor.
- **Binding classes (step 27, DiffBind, background normalization).** Gained: FDR < 0.05 up in D18 vs D15 or D22 vs
  D15. Stable-bound: every other site. The sites lost over time (12) are listed, not tested.
- **Binding change per KO-effect group.** Sites are members of a set and the 9 libraries are the replicates: limma
  voom on the DiffBind site counts (background size factors, `~ batch + day`), `camera` for each group against all
  other sites, contrasts D18 vs D15, D22 vs D18 and D22 vs D15. Figure: cumulative distribution of the DiffBind site
  log2FC per group against all sites.
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
| 00 DiffBind counts and size factors | `scripts/00_diffbind_counts.sb` (amem, 12 CPUs, 240 GB; re-reads the BAMs for the background bins) | `results/00_diffbind_counts/`: `counts_raw.tsv`, `norm_factors.tsv`, `counts_norm.tsv` |
| 01 GREAT domains, site -> gene | `scripts/01_great_site_to_gene.sb` (acpu, 2 CPUs, 7 GB) | `results/01_great/`: `gene_regions.bed`, `site_gene.tsv`, `summary.tsv` |
| 02 Direct targets per KO-effect group | `scripts/02_direct_targets.sb` (acpu, 2 CPUs, 7 GB) | `results/02_direct_targets/`: `gene_table.tsv` (group, KO log2FC per day, sites per gene), `enrichment.tsv` |
| 03 Binding change per KO-effect group | `scripts/03_binding_change.sb` (acpu, 2 CPUs, 7 GB) | `results/03_binding_change/`: `camera.tsv`, `site_stats.tsv` (per site: DiffBind and voom log2FC, FDR, gained), `site_groups.tsv` |
| 04 Candidate time-dependent direct targets | `scripts/04_candidates.sb` (acpu, 2 CPUs, 7 GB) | `results/04_candidates/`: `candidate_genes.tsv`, `candidates.tsv` (every matching site), `summary.tsv`, `gained_enrichment.tsv` |
| 05 Binding signal per day at each cluster's sites | `scripts/05_binding_signal.sb` (acpu, 2 CPUs, 7 GB) | `results/05_binding_signal/`: `site_signal.tsv`, `group_signal.tsv`, `wilcoxon.tsv` |
| 06 Alpha-leaning vs beta-leaning binding per day, libraries as replicates | `scripts/06_cluster_contrast.sb` (acpu, 2 CPUs, 7 GB) | `results/06_cluster_contrast/`: `library_values.tsv`, `tests.tsv`, `set_sizes.tsv` |
| 07 Target genes of the differentially bound sites and their KO effect | `scripts/07_diff_site_targets.sb` (acpu, 2 CPUs, 7 GB) | `results/07_diff_site_targets/`: `diff_site_genes.tsv` (every differential site and target gene), `gene_classes.tsv`, `summary.tsv`, `tests.tsv` |

Step 02 reads `timecourse_clusters_LRT_interaction.tsv`, which the bulk report builder writes on the Mac; it is
copied to `bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2/` on Alpine.

Report (Mac, from this folder): `bash scripts/fetch_results.sh`, then `~/anaconda3/bin/python3 build_report.py`
-> `results/report_integration.html` (steps 01-05: tables, bound share per group, cumulative distributions of the
site log2FC per group, candidate lists, boxplots of the binding signal per day and cluster; section 06: alpha-leaning vs
beta-leaning sites per day, and binding and expression of their marker genes, using the bulk normalized counts;
section 07: KO effect of the target genes of the differentially bound sites).

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
| No KO effect | 7,631 | 0.317 (0.345), down 0.85 | -0.006 (-0.060), **up 0.018** | 0.289 (0.229), up 0.40 | 4.1 (6.6) |
| Cluster 1 | 850 | 0.330 (0.330), down 0.85 | -0.043 (-0.032), down 0.91 | 0.235 (0.265), down 0.67 | 4.7 (5.3) |
| Cluster 2 | 629 | 0.343 (0.329), up 0.85 | -0.060 (-0.032), down 0.73 | 0.241 (0.264), down 0.90 | 7.2 (5.2) |
| Cluster 3 | 390 | 0.307 (0.330), down 0.85 | -0.038 (-0.033), up 0.94 | 0.253 (0.263), up 0.95 | 5.1 (5.3) |
| **Cluster 4** | 613 | **0.388** (0.328), up 0.21 | **-0.084** (-0.030), **down 0.018** | 0.224 (0.264), down 0.67 | **9.8** (5.1) |
| Cluster 5 | 219 | 0.301 (0.330), down 0.85 | -0.021 (-0.033), up 0.98 | 0.241 (0.263), down 0.67 | 2.7 (5.3) |
| Cluster 6 | 218 | 0.324 (0.330), up 0.85 | -0.054 (-0.032), down 0.24 | 0.213 (0.264), down 0.40 | 3.7 (5.3) |
| Stays | 3,077 | 0.324 (0.331), down 0.85 | -0.027 (-0.034), up 0.98 | 0.258 (0.264), down 0.95 | 4.9 (5.3) |
| Other KO | 625 | 0.332 (0.329), down 0.85 | -0.020 (-0.033), up 0.91 | 0.268 (0.263), up 0.95 | 4.6 (5.3) |

- The whole site set rises from D15 to D18 (median log2FC about 0.33) and stays there (D22 vs D18 about -0.03); the
  groups follow that rise.
- **Cluster 4** (KO effect down, emerging and growing) rises most from D15 to D18 (FDR 0.21) and is the one cluster
  whose sites then **fall from D18 to D22 relative to all other sites** (median -0.084 vs -0.030, FDR 0.018). Its
  binding peaks at D18, when its KO effect first appears. The genes without a KO effect move the opposite way
  (FDR 0.018), as the complement of the KO-effect clusters.
- Sites move together: the inter-site correlation estimated by camera is 0.19-0.28, far above limma's default 0.01.
  With the estimated correlation every FDR is >= 0.96 (D22 vs D18: 0.995), so the cluster 4 results are leads, not
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

Signal: DiffBind reads per site / background size factor (step 00), mean of the 3 libraries per day. Paired
Wilcoxon between days within each group (BH).

| Group | Sites | Median signal D15 / D18 / D22 | Sites up D18 vs D15 | Sites up D22 vs D18 | Median log2 ratio D22 vs D18 (FDR) |
|---|---|---|---|---|---|
| No KO effect | 7,631 | 23.6 / 31.8 / 32.4 | 84.9% | 51.6% | +0.02 (0.054) |
| Cluster 1 | 850 | 23.0 / 31.0 / 29.1 | 83.4% | 41.5% | -0.13 (5.7e-12) |
| Cluster 2 | 629 | 22.0 / 31.1 / 28.8 | 85.1% | 36.9% | -0.18 (1.4e-16) |
| Cluster 3 | 390 | 22.6 / 30.4 / 29.8 | 84.6% | 43.6% | -0.07 (1.4e-3) |
| **Cluster 4** | 613 | 22.5 / **33.3** / 27.6 | 87.3% | **29.2%** | **-0.28** (7.9e-35) |
| Cluster 5 | 219 | 23.1 / 31.7 / 30.7 | 87.2% | 46.6% | -0.04 (0.035) |
| Cluster 6 | 218 | 22.7 / 31.0 / 28.4 | 87.2% | 36.7% | -0.16 (5.8e-7) |

- Binding rises from D15 to D18 at 83-87% of the sites of every group (all FDR < 1e-28), the overall ISL1 increase.
- **From D18 to D22 the sites near the KO-effect clusters lose signal while the reference stays flat** (51.6% of
  sites up, median ratio +0.02). Cluster 4 shows both the largest rise to D18 (median 22.5 -> 33.3) and the largest
  fall to D22 (71% of its sites down). Per site, none of these D22 vs D18 changes reaches significance in step 27.
  Against all other sites, cluster 4's fall is the one that stands out (step 03 camera D22 vs D18, FDR 0.018, not
  significant with the estimated inter-site correlation). The batch-aware model gives a smaller shift than these
  day means (median DiffBind log2FC -0.084 vs -0.030 for the other sites).
- The paired tests use sites as units; with hundreds of sites per group even small shifts give tiny p values
  (cluster 5's D22 vs D18 FDR of 0.035 comes with a median change of -0.04). They describe each group; whether a
  group differs from the others is the camera test of step 03.

### 06 Alpha-leaning against beta-leaning sites, per day and over time

Marker genes (human islet cell-type markers): `markers.tsv`, 28 alpha and 34 beta. Cluster 4 holds 12 alpha (GCG,
ARX, IRX2, TTR, LOXL4, DPP4, GLS, RGS4, PEMT, PLCE1, TMEM236, SLC7A2) and 5 beta markers, with SST, PYY and GAST;
cluster 2 holds 7 beta (INS, MAFA, PCSK1, ABCC8, PFKFB2, TSPAN1, ASCL2) and 6 alpha markers (MAFB, IRX1, CRYBA2,
PCSK2, CD36, KCNJ6), with ISL1 itself and PAX6. So cluster 4 leans alpha and cluster 2 leans beta; neither is one
cell type.

Test (step 06): for each library, the signal of two site sets (total, or median per site, which does not depend on
how many sites a set has), log2 ratio, fitted on ~ batch + day with the 9 libraries as replicates. Levels compare the
two sets on one day (between different site sets they also reflect site number and intrinsic strength); changes
compare days. p values uncorrected.

| log2 alpha-leaning / beta-leaning | Sites | D15 | D18 | D22 | Change D22 vs D15 |
|---|---|---|---|---|---|
| Cluster 4 / cluster 2, total | 613 / 629 | -0.01 (p 0.83) | +0.07 (p 0.14) | -0.01 (p 0.87) | +0.00 (p 0.97) |
| Cluster 4 / cluster 2, median per site | 613 / 629 | +0.03 (p 0.62) | +0.08 (p 0.23) | -0.02 (p 0.77) | -0.05 (p 0.58) |
| **Alpha / beta markers, median per site** | 45 / 52 | **-0.01 (p 0.92)** | **-0.13 (p 0.27)** | **-0.26 (p 0.048)** | -0.25 (p 0.14) |
| Alpha markers in cluster 4 / beta markers in cluster 2, median per site | 25 / 10 | +0.10 (p 0.54) | -0.24 (p 0.18) | -0.29 (p 0.11) | -0.39 (p 0.13) |

- **Binding at alpha and beta marker sites is the same at D15; beta sites then pull ahead, and at D22 they bind
  more than alpha sites** (median per site 26.0 vs 32.1; log2 -0.26, p 0.048, uncorrected). Binding is not stronger
  at alpha genes at D15 or D18. The relative shift over time (-0.25 from D15 to D22) is not significant (p 0.14).
- The whole clusters show no crossover: cluster 4 is about 5% above cluster 2 at D18 (p 0.14) and equal at D15 and
  D22.
- Gene by gene, binding near alpha markers of cluster 4 (LOXL4, DPP4, TTR, RGS4) and near beta markers of cluster 2
  (MAFA, PCSK1, ABCC8) peaks at D18 while their control expression keeps rising to D22. GCG binding keeps rising
  (48 / 76 / 84) and ARX is level after D18 (31 / 42 / 41). INS has no linked site.
- In control islets the alpha genes and ISL1 keep rising from D18 to D22 (GCG 35,545 -> 203,556; ARX 485 -> 1,020;
  ISL1 373 / 1,990 / 4,425 at D15 / D18 / D22), so the fall in binding at alpha sites is not from less ISL1 or fewer
  alpha-lineage cells.

### 07 Target genes of the differentially bound sites and their KO effect

The 758 gained sites link to 733 bulk-tested genes and the 12 lost sites to 14 (805 site-gene pairs). Genes are
classed by their sites: gained (>= 1 gained site), lost, only stable sites, no site. KO effect per day: KO vs Ctr
padj < 0.05, lower in the KO (ISL1 activates) or higher (ISL1 represses).

| Genes | n | KO effect any day | Lower / higher in KO, D15 | D18 | D22 | KO effect changes over time |
|---|---|---|---|---|---|---|
| With a gained site | 733 | 46.8% | 7.0 / 7.9% | 21.0 / 13.5% | 22.4 / 11.1% | 22.4% |
| Only stable sites | 9,309 | 36.4% | 5.2 / 6.4% | 14.1 / 12.1% | 11.4 / 9.8% | 14.0% |
| No site | 6,490 | 34.0% | 5.0 / 5.3% | 13.9 / 11.0% | 9.8 / 9.0% | 11.6% |
| With a lost site | 14 | 50.0% | 7.1 / 7.1% | 21.4 / 7.1% | 7.1 / 14.3% | 14.3% |

Gained against only-stable genes, logistic regression with domain length and number of sites (BH):

| Outcome | Gained / stable | Unadjusted OR | Adjusted OR (FDR) |
|---|---|---|---|
| KO effect on D15 | 14.9 / 11.7% | 1.32 | 0.99 (0.96) |
| KO effect on D18 | 34.5 / 26.2% | 1.49 | 1.20 (0.08) |
| **KO effect on D22** | 33.4 / 21.2% | 1.86 | **1.34 (0.007)** |
| KO effect on any day | 46.8 / 36.4% | 1.54 | 1.17 (0.09) |
| KO effect changes over time | 22.4 / 14.0% | 1.77 | 1.24 (0.08) |
| Lower in KO, among genes significant on D22 | 66.9 / 53.7% | 1.74 | 1.35 (0.09) |

- **Genes that gain ISL1 binding from D15 to D18 become ISL1-dependent later:** no excess KO effect at D15, a small
  one at D18 and a clear one at D22 (adjusted OR 1.34, FDR 0.007). Binding is gained before the KO effect is
  strongest.
- Much of the unadjusted difference comes from domain size and number of sites (unadjusted OR 1.86 vs adjusted 1.34
  at D22).
- At D22 the gained-site genes with a KO effect are mostly lower in the KO (67% vs 54% for stable genes; adjusted
  OR 1.35, FDR 0.09): gained binding goes with ISL1 activating its targets.
- 343 genes have a gained site and a KO effect; by KO-effect group: stays 155, cluster 4 54, cluster 2 41, cluster 1
  38, other KO 24, cluster 3 19, clusters 5 and 6 6 each. Strongest at D22: ERO1B, PCP4, LMO2, TMOD1, EDARADD, ETV1,
  TMCC3, CLU, GLS, CHRNA3, DPP4, LOXL4 (all lower in the KO).
- Lost sites are too few (14 genes; e.g. GATA3, NR2E3, TGFBR3) to test.
- Binding at the 313 gained sites of the 343 genes (report boxplot): median 16.5 / 44.7 / 35.3 at D15 / D18 / D22,
  lower at D22 than at D18 for 87% of the sites (all sites: 55%), for genes lower (202 sites) and higher in the KO
  (118 sites) alike. The rise from D15 is how the sites were selected; part of the fall to D22 is regression to the
  mean of sites picked for a high D18 reading.
