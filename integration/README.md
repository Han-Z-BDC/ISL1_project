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
| ISL1 CUT&RUN, control line, D15 / D18 / D22, 3 libraries per day | `cutnrun_ISL1_timecourse` step 06 (sites) and step 07 (IgG greylist, DiffBind), on scratch | 14,431 sites of 401 bp; per-site differential tables (main analysis: DiffBind, background normalization); DiffBind reads per site (`counts_diffbind.tsv`) and size factors (`norm_factors.tsv`, method diffbind, background) |
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
- **Direct target.** A gene the ISL1 KO changes in the bulk RNA-seq (padj < 0.05 in KO vs Ctr at D15, D18 or D22, or
  in LRT_interaction) with at least one ISL1 site (a consensus peak on any day) in its GREAT domain. When the peak
  appears and whether it changes are not considered. **Direct-activated** = lower in the KO, **direct-repressed** =
  higher in the KO, per significant day; overall **mixed** when the days disagree; genes significant only in
  LRT_interaction take the direction of the day with the largest |shrunk KO log2FC| (flagged). ISL1 itself (the
  knocked-out gene) is excluded.
- **Counts and signal.** Every step uses the CUT&RUN step 07 DiffBind counts (paired-end fragments per site, duplicates
  kept) and its background size factors (DESeq2 on 15 kb genome bins). Signal per site and library =
  reads / size factor.
- **Binding classes (CUT&RUN step 07, DiffBind, background normalization).** Gained: FDR < 0.05 up in D18 vs D15 or D22 vs
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
| 01 GREAT domains, site -> gene | `scripts/01_great_site_to_gene.sb` (acpu, 2 CPUs, 7 GB) | `results/01_great/`: `gene_regions.bed`, `site_gene.tsv`, `summary.tsv` |
| 02 Direct targets among the KO-changed genes | `scripts/02_direct_targets.sb` (acpu, 2 CPUs, 7 GB) | `results/02_direct_targets/`: `direct_targets.tsv` (differential and bound genes with labels), `gene_table.tsv` (all genes: KO-effect group, KO log2FC per day, sites, labels), `summary.tsv` |
| 03 Binding change per KO-effect group | `scripts/03_binding_change.sb` (acpu, 2 CPUs, 7 GB) | `results/03_binding_change/`: `camera.tsv`, `site_stats.tsv` (per site: DiffBind and voom log2FC, FDR, gained), `site_groups.tsv` |
| 04 Candidate time-dependent direct targets | `scripts/04_candidates.sb` (acpu, 2 CPUs, 7 GB) | `results/04_candidates/`: `candidate_genes.tsv`, `candidates.tsv` (every matching site), `summary.tsv`, `gained_enrichment.tsv` |
| 05 Binding signal per day at each cluster's sites | `scripts/05_binding_signal.sb` (acpu, 2 CPUs, 7 GB) | `results/05_binding_signal/`: `site_signal.tsv`, `group_signal.tsv`, `wilcoxon.tsv` |
| 06 Alpha-leaning vs beta-leaning binding per day, libraries as replicates | `scripts/06_cluster_contrast.sb` (acpu, 2 CPUs, 7 GB) | `results/06_cluster_contrast/`: `library_values.tsv`, `tests.tsv`, `set_sizes.tsv` |
| 07 Target genes of the differentially bound sites and their KO effect | `scripts/07_diff_site_targets.sb` (acpu, 2 CPUs, 7 GB) | `results/07_diff_site_targets/`: `diff_site_genes.tsv` (every differential site and target gene), `gene_classes.tsv`, `summary.tsv`, `tests.tsv` |
| 08 All ISL1 sites grouped by binding pattern across days | `scripts/08_site_clusters.sb` (acpu, 4 CPUs, 15 GB) | `results/08_site_clusters/`: `site_clusters.tsv` (per site: LRT, VST and z per day, group for k = 2-10), `centers.tsv`, `elbow.tsv`, `heldout.tsv`, `heldout_clusters.tsv`, `links.tsv`, `gene_links.tsv` and `gene_site_groups.tsv` (gene level, GREAT and ±50 kb links), `vst_libraries.tsv` |
| 09 Target genes of the peak-call overlap classes | `scripts/09_peak_overlap.sb` (acpu, 2 CPUs, 7 GB; needs CUT&RUN step 12) | `results/09_peak_overlap/`: `class_groups.tsv` (per class: step 08 VST and groups), `gene_classes.tsv` and `gene_class_links.tsv` (gene level, GREAT and ±50 kb) |
| 10 Pre-specified beta − alpha binding contrast | `scripts/10_beta_alpha_contrast.sb` (acpu, 1 CPU, 3 GB; `MODE=dry_run` now, `primary` for the unified time course) | `results/10_beta_alpha_contrast/`: `sets.tsv`, `set_sites.tsv`, `library_values.tsv`, `contrasts.tsv`, `site_level.tsv` |
| 11 Direct-activated and direct-repressed genes with a stage-specific KO effect | `scripts/11_direct_stage.sb` (acpu, 2 CPUs, 7 GB; needs step 02; runs `11_direct_stage.R` for each label) | `results/11_activated_stage/`, `results/11_repressed_stage/`: `genes.tsv` (per gene of the label: pattern, LRT_interaction padj, scaled profile, clusters k2-k8, group: the k = 4 cluster or constant), `elbow.tsv`, `centers.tsv`, `patterns.tsv`, `groups.tsv` (per group: genes, median log2FC, significant on all days) |
| 12 ISL1 binding of the D15 scRNA-seq trajectory gene clusters | `scripts/12_trajectory_binding.sb` (acpu, 4 CPUs, 15 GB; needs CUT&RUN steps 08 and 13, scRNAseq 7.3h / 7.3l) | `results/12_trajectory_binding/`: `genes.tsv` (per gene: cluster or background, bound per stage, TSS-bound), `binding_by_cluster.tsv` (fraction bound, odds ratio, BH FDR vs background), `tss.bed`, `tss_matrix.gz` (deepTools), `tss_signal.tsv` (per gene and day, ±500 bp), `genes_mode.tsv` and `binding_mode.tsv` (binding mode per stage: fraction, Fisher and adjusted odds ratios), `genes_timing.tsv`, `timing_tests.tsv`, `timing_classes.tsv`, `timing_overall.tsv` (first-bound stage, cluster x stage), `within_clusters.tsv`, `within_chisq.tsv` (each cluster vs the other clusters), `isl1_binding_by_cluster.pdf`, `isl1_binding_mode.pdf`, `isl1_target_timing.pdf`, `isl1_target_within_clusters.pdf`, `isl1_tss_heatmap.pdf`, `isl1_tss_profile.pdf` |

Step 02 reads `timecourse_clusters_LRT_interaction.tsv`, which the bulk report builder writes on the Mac; it is
copied to `bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2/` on Alpine.

Report (Mac, from this folder): `bash scripts/fetch_results.sh`, then `~/anaconda3/bin/python3 build_report.py`
-> `results/report_integration.html` (steps 01-05: tables, direct targets with their labels, cumulative distributions of the
site log2FC per group, candidate lists, boxplots of the binding signal per day and cluster; section 06: alpha-leaning vs
beta-leaning sites per day, and binding and expression of their marker genes, using the bulk normalized counts;
section 07: KO effect of the target genes of the differentially bound sites; section 08: all sites grouped by
binding pattern, with the held-out check and the gene level; section 09: target genes of the peak-call overlap classes
(the classes themselves: CUT&RUN step 12); section 10: the pre-specified beta − alpha contrast).

## Results

### 01 GREAT domains, site -> gene

- 36,559 genes in the GTF; 16,546 of the 16,563 bulk-tested genes have a domain (the other 17 are 13
  mitochondrial genes and 4 on unplaced contigs). Domains average 180 kb; none is over 2 Mb.
- 14,430 of the 14,431 sites link to at least one gene: 3,580 to one gene, 7,639 to two, 3,211 to three or more
  (basal domains overlap). 29,103 site-gene links, 8,790 of them in a basal domain; median distance to the TSS
  23.4 kb.
- **10,056 of the tested genes (60.8%) have at least one site**, a median of 1. Bound-gene fractions of the
  KO-effect groups are therefore compared with this background, not read on their own.

### 02 Direct targets among the KO-changed genes

6,032 genes are differential (8 more have no GREAT domain and are left out; ISL1 itself is excluded); 3,813 of them
(63%) are direct targets.
A gene counts once per day in the per-day rows.

| Bulk test | Direct-activated (of differential) | Direct-repressed (of differential) |
|---|---|---|
| KO vs Ctr, D15 | 539 (862) | 657 (998) |
| KO vs Ctr, D18 | 1,471 (2,376) | 1,223 (1,938) |
| KO vs Ctr, D22 | 1,226 (1,863) | 997 (1,578) |
| KO effect changes over time only | 17 (29) | 20 (35) |
| **Overall** | **1,977 (3,168)** | **1,702 (2,673)**; plus 134 mixed (191) |

- Islet regulators: direct-activated ARX, GCG, MAFA (D22), PDX1 (D18, D22); direct-repressed NEUROG3, NKX6-1, PAX4; MAFB mixed (higher in the KO at D15, lower at D22).
- About 60% of all genes have a site in their GREAT domain, so a site alone separates direct from indirect effects
  weakly; the share is similar for activated and repressed genes and on every day (62-66%).

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

Signal: DiffBind reads per site / background size factor (CUT&RUN step 07), mean of the 3 libraries per day. Paired
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
  fall to D22 (71% of its sites down). Per site, none of these D22 vs D18 changes reaches significance in CUT&RUN step 07.
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

### 08 All ISL1 sites grouped by binding pattern across days

All 14,431 sites, CUT&RUN step 07 DiffBind counts and background size factors. DESeq2 LRT (~ batch + day against ~ batch)
selects the sites whose binding varies across days (p < 0.05, a noise filter: 2,747 sites; 408 at FDR < 0.1); the
other 11,684 are the stable group. Per site: VST, batch 3 effect removed, mean per day, z-score across days; k-means
with k = 4 (share of variance between groups 0.81; k = 5-6 split the same shapes).

| Group | z D15 / D18 / D22 | Sites | Gained (CUT&RUN step 07) | Sites linked to (OR vs stable, FDR) |
|---|---|---|---|---|
| 1 rises by D18, stays | -1.13 / 0.57 / 0.56 | 1,145 | 35.0% | beta markers 3.0 (0.02); KO cluster 4 1.6 (0.003) |
| 2 peaks at D18, falls by D22 | -0.92 / 1.02 / -0.10 | 909 | 36.6% | alpha markers 3.6 (0.01); KO cluster 4 1.9 (3e-4) |
| 3 late, rises at D22 | -0.61 / -0.45 / 1.06 | 419 | 5.5% | none |
| 4 falls at D22 | 0.14 / 0.79 / -0.93 | 274 | 0.4% (10 of the 14 lost sites) | none |

- In the 9 libraries the groups fit early alpha / late beta binding: alpha-marker sites peak at D18 and fall, beta-
  marker sites rise and stay. Marker sites are few (9 and 10).
- **Only the rise after D15 replicates in an independent batch.** Clustering the 6 earlier libraries gives the same
  shapes; the 3 batch 3 libraries rise from D15 in every group but do not follow the D18-vs-D22 shapes (mean batch 3
  profile vs group centre: r 0.89 for "rises and stays", 0.44, 0.11 and -0.28 for the others), and single sites
  match their own group no more often than chance. The earlier libraries come from another differentiation with
  relabelled days (D12 as D15, D15 as D18). The D18-vs-D22 distinction between groups awaits the unified CUT&RUN
  time course; the same caveat applies to the D22 alpha-vs-beta difference of step 06.

Gene level: genes linked to each group's sites against genes whose linked sites are all stable, with two links:
GREAT domains (main) and every gene whose TSS lies within ±50 kb of the site (check). Genes with a varying site have
more sites (median 2 vs 1), so the main figure is the adjusted odds ratio from a logistic regression with
log2(1 + number of sites) (and, for GREAT, log10 domain length) as covariates; BH FDR per link. Tables:
`gene_links.tsv` (per link, site group and KO group) and `gene_site_groups.tsv` (per gene: sites in each group).

| KO gene group | Link | Group 1 | Group 2 | Group 3 | Group 4 |
|---|---|---|---|---|---|
| Cluster 4 (down, emerging and growing) | GREAT | **2.2 (1e-4)** | **2.0 (0.005)** | 1.0 (0.99) | 1.5 (0.83) |
| | ±50 kb | 1.5 (0.19) | **2.9 (2e-5)** | 1.0 (0.96) | **4.0 (6e-4)** |
| Cluster 2 (up early, down by D22) | GREAT | 1.2 (0.83) | 1.2 (0.83) | 1.3 (0.83) | 1.1 (0.93) |
| | ±50 kb | 1.0 (0.96) | **1.8 (0.03)** | 1.0 (0.96) | 1.6 (0.51) |

Adjusted odds ratio (adjusted FDR).

- KO cluster 4 genes are enriched near the varying sites under both links, most consistently near group 2 (peaks at
  D18). KO cluster 2's enrichment (Fisher OR about 2 with GREAT) is mostly explained by the number of sites. Group 3
  (late rise) links to no KO group.
- 279 genes of KO clusters 2 and 4 have a varying site under either link (listed in the report), among them PDX1
  (groups 1 and 2), GCG, PCSK1, ABCC8, DPP4, TTR and, with GREAT only, MAFB (group 2).
- Key genes: within 50 kb, INS, PAX6, PAX4, ARX, MAFA and MAFB have only stable sites, and NKX6-1 one group 1 site.
  With GREAT, INS and PAX6 have no site, PAX4, ARX and MAFA only stable sites, NKX6-1 group 1 and 2 sites.

### 09 Target genes of the peak-call overlap classes

The overlap classes come from CUT&RUN step 12 (`cutnrun_ISL1_timecourse`; class sizes, library call rates and the check
of peak presence against the counts are in that folder's README and report). The classes match the step 08 groups:
D18-only sites are enriched in group 2 (peaks at D18), D22-only sites in group 3 (late rise), and the D15 + D18 and
D15-only sites in group 4 (falls at D22).

Target genes per class: bulk-tested genes linked by GREAT (main) or with the TSS within ±50 kb (check); a gene counts in
every class in which it has a site and is compared with all other linked genes; logistic regression with log2(1 +
number of sites) (and, for GREAT, log10 domain length) as covariates; BH FDR per link.

| Class | Genes (GREAT / ±50 kb) | Adjusted odds ratio (adj. FDR < 0.05) |
|---|---|---|
| all three days | 2,493 / 3,236 | none (GREAT cluster 4 1.46, FDR 0.14) |
| D18 + D22 | 4,160 / 5,521 | none |
| D18 only | 3,386 / 3,164 | no KO effect 0.83 (GREAT), 0.77 (±50 kb); ±50 kb: cluster 1 1.38, cluster 4 1.59, stays 1.21 |
| D22 only | 3,520 / 4,970 | ±50 kb: cluster 4 0.64, no KO effect 1.13 |
| D15 + D18, D15 + D22, D15 only | 160-360 | none (small, unstable odds ratios) |

- No overlap class maps onto one KO cluster. Genes with a D18-only site are more often KO-affected (cluster 4 and
  cluster 1 within 50 kb); genes with a D22-only site are less often in cluster 4. The Fisher enrichments of clusters
  2 and 4 with GREAT are explained by the number of sites.
- Examples (sites in one class only): MAFA (cluster 2) D18 only; NEUROG3 (cluster 6) D22 only; ISL1 (cluster 2)
  D15 + D22; PAX4 (cluster 2), ARX (cluster 4) all three days. All KO cluster 1-6 genes with their classes are listed
  in the report.

### 10 Pre-specified beta − alpha contrast

Fixed on 2026-11-15 (commit 3e89f03), before the unified CUT&RUN time course is analysed; full definition in the header
of `scripts/10_beta_alpha_contrast.R`.

- Hypothesis: ISL1 binding at beta-program genes gains on binding at alpha-program genes late in differentiation.
- Gene sets: `markers.tsv` (28 alpha, 34 beta). Sites: centre within ±50 kb of a marker gene's TSS (29 alpha sites near
  21 genes, 44 beta sites near 23 genes; none shared).
- Per library: y = log2(median signal at beta sites / median at alpha sites), signal = DiffBind reads / background size
  factor. Model: lm(y ~ batch + day), or lm(y ~ day) for one batch.
- **Primary: D22 − D18, one-sided (beta gains), α = 0.05.** Secondary: D22 − D15, D18 − D15, GREAT-linked sites, site-level
  DiffBind log2FC (descriptive).
- The primary test is for the unified time course (`MODE=primary`). The current 9 libraries give a dry run only, since
  step 06 already examined beta vs alpha marker sites in them.

Dry run on the current 9 libraries:

| Sites | Contrast | Estimate (log2) | SE | One-sided p |
|---|---|---|---|---|
| ±50 kb (primary) | D22 vs D18 | -0.03 | 0.19 | 0.55 |
| ±50 kb | D22 vs D15 | 0.08 | 0.19 | 0.34 |
| GREAT | D22 vs D18 | 0.14 | 0.14 | 0.19 |
| GREAT | D22 vs D15 | 0.25 | 0.14 | 0.07 |

- No gain of beta on alpha between D18 and D22 in these libraries. Libraries of the same day differ by up to 0.8 log2
  (D18), more than any shift between days; the GREAT version reproduces step 06.

Replicates needed (planning, `power.tsv`): residual SD of the dry-run primary fit 0.23 log2 (5 df; upper 80% bound 0.40).
In a one-batch design with n libraries per day, the D22 − D18 shift detectable with 80% power (one-sided, α = 0.05) is
0.52 log2 for n = 3, 0.43 for n = 4, 0.38 for n = 5, 0.34 for n = 6 and 0.29 for n = 8 (0.92 / 0.76 / 0.67 / 0.60 / 0.51 at
the upper SD bound). A 0.3 shift needs about 8 libraries per day.

### 11 Direct targets with a stage-specific KO effect

Same method for each label; genes labelled mixed are left out. Direct-activated:

- Genes: the 1,977 direct-activated genes of step 02 (lower in the KO on at least one day, never higher, bound).
- Stage-specific: LRT_interaction padj < 0.05 (the KO effect differs between stages): 857 genes. The other 1,120 have
  a KO effect not shown to change between stages (e.g. G6PC2, SLC30A8). The per-day significance pattern is
  descriptive only.
- Clusters: k-means (seed 1, 25 starts) on each stage-specific gene's shrunk KO log2FC at D15 / D18 / D22 divided by its
  largest |log2FC| (timing, not size); k = 4, chosen from k = 2-8 (78% of the variance between clusters; 71% at k = 3).
- Constant group: all genes that are not stage-specific, kept as a fifth group.

| Group | Genes | Median KO log2FC D15 / D18 / D22 | Significant on all 3 days | Islet genes |
|---|---|---|---|---|
| 1 Sustained, from D18 | 345 | -0.21 / -1.04 / -1.66 | 109 (32%) | ARX, GCG, HHEX, PDX1, PFKFB2, SST, TTR |
| 2 Transient, D18 only | 223 | -0.28 / -1.25 / -0.21 | 42 (19%) | |
| 3 Late, D22 only | 199 | +0.11 / -0.05 / -1.00 | 0 | ABCC8, CD36, MAFA, PCSK2 |
| 4 Early, transient | 90 | -0.29 / -0.44 / +0.20 | 0 | |
| Constant (not stage-specific) | 1,120 | -0.18 / -0.49 / -0.39 | 112 (10%) | G6PC2, SLC30A8 |

- Sustained, from D18: loss from D18 on, still there at D22; Transient, D18 only: loss at D18, gone by D22; Late,
  D22 only: loss at D22 only; Early, transient: small loss at D15-D18, none at D22. A loss mainly at D15 is rare (22 stage-specific genes lose most at D15).
- Constant: on average not flat (largest loss at D18), but only 112 genes are significant on all three days; the
  others are significant on one or two days, with a KO effect not shown to differ between stages.
- INS is not here: it is lower in the KO but has no ISL1 site in its GREAT domain.

Direct-repressed: 1,702 genes (higher in the KO on at least one day, never lower, bound); 478 stage-specific (LRT_interaction
padj < 0.05), 1,224 not shown to change between stages (constant group). k = 4 (77% of the variance between clusters;
71% at k = 3).

| Group | Genes | Median KO log2FC D15 / D18 / D22 | Significant on all 3 days | Endocrine genes |
|---|---|---|---|---|
| 1 Sustained, from D18 | 164 | +0.04 / +0.73 / +0.88 | 31 (19%) | NKX6-1 |
| 2 Early, transient | 147 | +0.82 / +0.50 / -0.03 | 21 (14%) | ASCL2, CHGB, CRYBA2, KCNJ6, MNX1, NEUROD1, NEUROG3, NKX2-2, PAX4 |
| 3 Late, D22 only | 103 | -0.15 / -0.03 / +0.53 | 0 | |
| 4 Transient, D18 only | 64 | +0.23 / +0.70 / -0.05 | 2 (3%) | |
| Constant (not stage-specific) | 1,224 | +0.21 / +0.47 / +0.36 | 212 (17%) | DDC, INSM1, LMX1A |

- Sustained, from D18: gain from D18 on, still there at D22; Early, transient: gain at D15, fading by D22 (the
  endocrine progenitor genes); Late, D22 only: gain at D22 only; Transient, D18 only: gain at D18 only. Unlike the activated genes, an early effect is common
  here (132 stage-specific genes gain most at D15).

### 12 ISL1 binding of the D15 scRNA-seq trajectory gene clusters

User, 2026-11-24: how many genes of each D15 trajectory gene cluster (scRNAseq 7.3h, 3,118 KO vs control genes,
k = 10) have an ISL1 site, per CUT&RUN stage and at any stage, and the ISL1 signal at their TSS.

- Binding = CUT&RUN step 13 links: a merged ISL1 site (pooled macs3_ge2 consensus of D15 / D18 / D22) within ±1 kb of
  any TSS of the gene, or overlapping an ABC element linked to it (ENP2, ENP3, early_SC_alpha; "enhancer only" when no
  TSS link). Per stage: the site is called on that day; any stage: any linked site.
- Background: the 11,374 other genes tested along the D15 lineages (7.3l), expressed but not called. Fisher's exact
  test per cluster and stage on bound / not bound, BH over all cluster x stage tests (separately for TSS-bound).
  20 genes without a GTF gene TSS (mitochondrial, renamed duplicates) dropped.
- TSS signal: one TSS per gene (GTF gene start), ±3 kb in 50 bp bins, deepTools computeMatrix on the per-day mean ISL1
  bigWigs (CUT&RUN step 08); the clusters plus 1,000 random background genes (seed 1234). Heatmap (user, 2026-11-24):
  clusters only, genes sorted within a cluster by mean signal over ±500 bp and the three days; sized to sit next to the
  scRNAseq gene-pattern heatmap (panel area 3.5 in tall, cluster heights in proportion to their genes, each day as
  wide as one beta / alpha curve there, 0.41 in). The mean profiles keep the background (dashed). The bigWigs are copied to node-local disk first (read from
  /scratch, the parallel workers stalled on I/O).

| Bound (TSS or enhancer) | any stage | D15 | D18 | D22 |
|---|---|---|---|---|
| background | 77% | 33% | 65% | 72% |
| C1 | 80% | **42%** | **72%** | 75% |
| C4 | 78% | **44%** | 70% | 72% |
| C10 | **63%** (lower) | 25% | **53%** (lower) | **58%** (lower) |
| C2, C3, C5-C9 | 73-79% | 30-37% | 60-69% | 64-72% (C8 64%, lower) |

- About three quarters of all expressed genes are linked to an ISL1 site, so the bound fraction of every cluster is
  high; the clusters differ only modestly from the background.
- More bound than the background: C1 (beta, ISL1-dependent: ISL1, PDX1, FOXO1, PBX3, LMO2) at D15 (OR 1.50, FDR
  0.002) and D18 (OR 1.38, FDR 0.04); C4 (alpha: ARX, ETV1, PAX6, HHEX) at D15 (OR 1.64, FDR 2e-4). The excess is in
  enhancer-only (ABC) links, not TSS links.
- Less bound: C10 (higher in KO, EC rise in KO) at any stage, D18 and D22 (OR 0.50-0.61); C8 at D22 (OR 0.68).
- TSS-bound alone is lower than the background in most clusters (OR 0.4-0.8; background 53% at any stage), and the
  mean TSS signal is highest in the background block: the background holds many broadly expressed genes with strong
  promoter signal. ISL1 TSS signal rises from D15 / D18 to D22 in every block.

Binding mode, TSS or enhancer (user, 2026-11-24: per stage only; four categories; plus a check for the number of
enhancers a gene has and its expression):

- Per gene and stage, from the sites called that day: TSS = a site within ±1 kb of a TSS of the gene (step 13 source
  TSS or TSS+ABC); enhancer = a site linked only through an ABC element (source ABC). Categories: TSS only, TSS +
  enhancer, enhancer only, all (any of them; the same genes as "bound" above).
- Fisher's exact test per cluster, stage and category against the background (category vs not, all genes); BH within
  a category over cluster x stage.
- Adjusted check: logistic regression per stage and category, category ~ cluster (background = reference) +
  log2(1 + ABC elements of the gene) + mean logCPM along the D15 lineages (7.3l); ABC elements = distinct
  non-self-promoter elements linked to the gene (ABC >= 0.02) in ENP2, ENP3 or early_SC_alpha (median 6 in the
  background, 7-9 in the clusters). 14,449 genes with a logCPM. More ABC elements raise the chance of an
  enhancer-only site (OR 2.3-2.7 per doubling) and lower TSS-only (0.5-0.7).
- Figure `isl1_binding_mode.pdf` (user, 2026-11-24: the bar plot only): stacked bars (TSS only / TSS + enhancer /
  enhancer only) per stage with the total bound and its Fisher stars, background dashed; stage panels 1.32 x 1.58 in,
  page 5.1 x 2.8 in. The odds ratios per category (Fisher and adjusted) are in `binding_mode.tsv`.

Results:
- Enhancer only is the category where clusters exceed the background (background 19% / 27% / 24% at D15 / D18 /
  D22). C4 (alpha) is the clearest: 33% at D15 (OR 2.1, adjusted 1.7, FDR 6e-4), 40% at D18; C1 (beta) 28% at D15
  (OR 1.7, adjusted 1.5, FDR 0.03). Their higher total binding at D15 comes from enhancer-only genes and holds after
  the adjustment (all bound: C1 adjusted OR 1.4, C4 1.5).
- At D18 / D22 most clusters (C1, C4, C6-C9) have more enhancer-only genes than the background in the Fisher test, but
  after the adjustment only C6, C7 and C9 keep it at one or two stages: part of the enhancer excess reflects genes
  with many ABC elements.
- TSS only is lower than the background in most clusters (C6, C8, C10 at every stage); TSS + enhancer is lower at D22
  after the adjustment in C3, C6-C10. C10 (higher in KO) is less bound in every category.
- So the ISL1-dependent beta (C1) and alpha (C4) clusters are bound early (D15) and mainly through distal enhancers;
  no cluster is enriched for TSS binding.

Enrichment among the ISL1 targets and its stage dependence (user, 2026-11-24). Among ISL1 targets (TSS or ABC
enhancer site), a cluster's enrichment is the same 2 x 2 Fisher test as above (the test is symmetric), so the new
question is whether it depends on the stage. Binding accumulates (most D15 sites stay at D18 / D22), so the per-stage
tests are not independent; two tests instead:

- First bound: each gene's first stage with a site (D15 / first at D18 / first at D22 / never). Per cluster: the class
  mix vs the background (chi-square on 2 x 4, BH over clusters) and each class vs the background (Fisher, BH within a
  class).
- Cluster x stage: bound per gene and stage ~ group * stage (background and D15 = reference), logistic GEE with
  independence working model and gene-clustered robust (sandwich) variance, written in base R (lme4 / geepack are not
  installed); Wald tests of all 20 cluster x stage terms jointly and of each cluster's 2 terms (BH over clusters).
- Figure `isl1_target_timing.pdf`: class mix per cluster (stacked) with the two test results, and the per-class odds
  ratios.

| Cluster | Bound from D15 | First at D18 | First at D22 | Never | Class mix FDR | Enrichment D15 (OR) | OR change D22 vs D15 | Interaction FDR |
|---|---|---|---|---|---|---|---|---|
| background | 33% | 33% | 11% | 23% | | | | |
| C1 | 42% | 31% | 7% | 20% | 0.002 | 1.50 | 0.77 | 0.11 |
| C4 | 44% | 27% | 7% | 22% | 2e-4 | 1.64 | 0.61 | 5e-4 |
| C10 | 25% | 29% | 9% | 37% | 5e-4 | 0.69 | 0.77 | 0.27 |

- The cluster x stage terms differ jointly (Wald chi2 = 64.4, 20 df, p = 1.5e-6): the clusters' enrichment among
  targets changes with stage.
- C1 (beta) and C4 (alpha) are enriched among the targets bound from D15 (OR 1.5 / 1.6) and depleted among those first
  bound at D22 (OR 0.6): their genes become targets early. For C4 the enrichment shrinks significantly by D22 (OR
  ratio 0.61, FDR 5e-4); C1 the same direction, not significant (0.77, FDR 0.11). As the background gains binding,
  the early advantage of these clusters fades.
- C2, C5, C6 and C8 also have a significant interaction (FDR 0.009-0.04), each with a lower odds ratio at D22 than at
  D15, without a significant class mix.
- C10 (higher in KO) is enriched among never-bound genes (37% vs 23%, OR 2.0, FDR 2e-4) at every stage (no
  interaction).

Each cluster against the other trajectory clusters (user, 2026-11-24): universe = the 3,118 C1-C10 genes only (no
background); per stage, target = a TSS or ABC-enhancer site called that day. Per cluster, is its share of the targets
different from its share of the genes: 2 x 2 Fisher's exact test, two-sided, cluster vs the other nine x target vs
not, BH over 10 clusters x 3 stages; per stage a 10 x 2 chi-square on the target fraction. Figure
`isl1_target_within_clusters.pdf` (stage panels 1.32 x 1.58 in, as `isl1_binding_mode.pdf`).

| Stage | Chi-square p (10 clusters) | All C1-C10 bound | Enriched (FDR) | Depleted (FDR) |
|---|---|---|---|---|
| D15 | 2.7e-6 | 35% | C1 42% (OR 1.40, 0.02), C4 44% (OR 1.55, 0.003) | C8 30% (OR 0.76, 0.046), C10 25% (OR 0.60, 0.02) |
| D18 | 5.2e-4 | 66% | C1 72% (OR 1.37, 0.046) | C10 53% (OR 0.56, 0.003) |
| D22 | 1.1e-3 | 68% | C1 75% (OR 1.45, 0.02) | C10 58% (OR 0.62, 0.02) |

- Among the trajectory genes, the target fraction differs between clusters at every stage. C1 (beta) is the most bound
  cluster at all three stages; C4 (alpha) only at D15 (at D18 / D22 70-72%, n.s.); C10 is the least bound at all
  stages; C8 is lower at D15 only. The same picture as against the background.

Re-clustering of the D15 targets (user, 2026-11-24: "only take all the D15 direct targets among C1-C10 and redo the
k-means clustering"): `scripts/12_recluster_targets.py` (CellRank env, as 7.3h) -> `results/12_trajectory_binding/
recluster_D15_targets/` (`k.csv`, `clusters_k<k>.csv`, `crosstab_k<k>.csv`, `curves_k<k>.png`, k 4-12). Genes: the
1,090 C1-C10 genes with a site called at D15; features and k-means exactly as 7.3h (smoothed binned curves, z per
gene across the six curves, n_init 50, seed 1234).

- Silhouette is highest at k = 4 (0.26) and falls with k (0.17 at k = 10); k = 4 and 5 are identical across seeds,
  k >= 7 less stable (ARI 0.63-0.90).
- k = 4 (new T1-T4, numbered by where the curves peak): T1 252 early genes falling on every lineage (C6 119, C5 76,
  C2 54); T2 372 late rise on beta and alpha in control, lost in KO (C4 163, C1 133, C3 59); T3 286 rise to the EC
  end, mid-path peak on beta / alpha kept in KO (C8 140, C3 69, C7 34, C10 27); T4 180 transient mid-path peak (C9
  112, C7 41, C10 18).
- k = 10 largely returns the original clusters (ARI 0.65 vs C1-C10): T1 = C1 114 + C4 27 (beta-led late rise, KO
  flat), T5 = C4 124 (alpha rise), T4 = C3 109 + C1 30, T8 = C8 127, T10 = C10 47; the D15 targets do not form new
  patterns. Among the targets, the largest single pattern at k = 4 is T2 (34%): the ISL1-dependent beta and alpha
  rise of C1 and C4 together.
- k = 10 figures (user, 2026-11-24: "use k = 10 and show me heatmap"): scRNAseq `fig_d15_gene_patterns.R` with
  GP_CLUSTERS = `recluster_D15_targets/clusters_k10.csv`, GP_PREFIX = T, GP_TAG = `_D15_targets_k10` ->
  `plot/D15/trajectory/recluster_D15_targets/gene_patterns_{heatmap,curves}_D15_targets_k10.pdf` (same style and size
  as the C1-C10 figures). Of the labelled lineage genes, D15 targets: PDX1, FOXO1, ISL1, LMO2, PBX3 (T1), ARX, PAX6
  (T5), MEIS3 (T6), PAX4 (T9, from C6); not bound at D15: ETV1, MAFB, FEV, MNX1, NEUROD4, NKX6-1, LHX1.

#### TFs and enrichment of the direct targets per cluster

`scripts/12_target_tf_enrich.R` / `.sb` (user, 2026-11-24): direct ISL1 target = a C1-C10 gene bound at any stage
(`genes.tsv` bound_any; no KO expression criterion). TFs = Lambert 2018 list (as scRNAseq 7.3h) -> `target_TFs.csv`
(per TF: cluster, bound / TSS-bound per stage, first-bound stage, sites) and `target_TF_counts.csv`. Over-representation
(fgsea::fora) in MSigDB 2026.1 (the saved scRNAseq table, as 7.1o): GO BP / MF / CC, KEGG Medicus, KEGG Legacy,
Reactome, Hallmark; universe = the 14,490 genes of `genes.tsv` (tested along the D15 lineages), sets of 15-500 universe
genes, BH per cluster x collection -> `target_enrich.csv`, `isl1_target_enrich.pdf` (one page per collection, top 5
terms per cluster with padj < 0.05).

| Cluster | Targets | TFs | Terms padj < 0.05 | Leading terms |
|---|---|---|---|---|
| C1 | 280 | 16 | 46 | calcium ion binding, regulation of hormone levels, Hallmark pancreas beta cells, insulin secretion |
| C2 | 132 | 11 | 57 | G2M checkpoint, E2F targets, chromosome segregation (cell cycle) |
| C3 | 269 | 16 | 142 | ER lumen, Reactome IGF transport / uptake by IGFBPs, hemostasis, synaptic membrane |
| C4 | 288 | 19 | 34 | GPCR ligand binding / signalling, peptide ligand receptors, peptide hormone metabolism |
| C5 | 182 | 14 | 10 | vasculature development, EMT, ECM-receptor interaction, axon guidance |
| C6 | 307 | 30 | 8 | EMT, KRAS signalling down, KEGG MODY, epithelial cell differentiation |
| C7 | 185 | 8 | 6 | motile cilium, microtubule-based movement |
| C8 | 350 | 20 | 17 | apoptosis, interferon gamma response, regulation of hormone levels, glycolysis / pentose phosphate |
| C9 | 246 | 9 | 38 | cilium movement, axoneme assembly, motile cilium |
| C10 | 118 | 11 | 2 | hypoxia, heterophilic cell-cell adhesion |

Paper figure (user, 2026-11-24: top terms across all categories per cluster, as a compact GO dot plot): scRNAseq
`analysis/scripts/fig_d15_target_enrich.R` (local) -> `plot/D15/trajectory/isl1_target_enrich_top.pdf` and
`target_enrich_top.csv`. Per cluster, the padj < 0.05 terms of all collections pooled by padj, a term skipped when its
target genes overlap an already chosen term's by Jaccard >= 0.5, top 6; columns = clusters, rows = the chosen terms
[database], a dot wherever padj < 0.05 (size = target genes, colour = -log10 padj capped at 10); no grid, cells
0.15 x 0.115 in. C7 and C10 have only 3 and 2 non-redundant terms.

Candidate TF summary (user, 2026-11-25: "show candidate TFs per cluster, focus on C1, C4 and C8 ... just a summary of
TF names"): `scripts/12_tf_summary.R` (run at the end of `12_target_tf_enrich.sb`) -> `isl1_target_TF_summary.pdf`
(3.4 x 1.45 in) and `isl1_target_TF_summary.csv`. Layout B (user): three columns of TF names, each with ISL1 binding at
D15 / D18 / D22 (filled = TSS site, open = enhancer only, grey = not bound); header = direct targets, TFs (ISL1 not
counted) and the within-cluster result. TFs chosen by the user from a table ranked by binding (bound at D15 or at >= 2
stages first) then the best-lineage FDR of the 7.3l pseudobulk KO test:

| Column | TFs |
|---|---|
| C1 | SIM1, PDX1, PBX3, FOXO1, TBX3 (top 5) |
| C4 | PAX6, ETV1, ARX, HHEX, CREB3L2 |
| C7 + C8 | FEV, ASCL2 (C7), LMX1A, MNX1, HEY2 (C8) |

ISL1 itself is left out (user: many of its binding sites come from the ISL1 cDNA). C7 and C8 are shown together
(user). FEV, MNX1 and HEY2 are first bound at D18; all five C7 + C8 TFs are enhancer-bound.
Layout D as well (user, 2026-11-25: "a plot that is narrow but tall"; then "8 pt, 80 wide, 150 tall, upright", no
"from D18" labels, no legend): `isl1_target_TF_network.pdf`, 80 x 150 pt, same TFs as a tree ISL1 -> cluster -> TF, C7
and C8 as separate nodes; TF names Arial regular 8 pt upright, node labels Arial bold 5 pt (cairo_pdf, Arial embedded;
the Arial files are in `~/.local/share/fonts` on Alpine). For the figure legend: node colour = within-cluster result
(red enriched at any stage, blue depleted, grey n.s.); TF dot filled = TSS site at any stage, open = enhancer only.

