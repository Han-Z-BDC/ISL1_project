# scRNA-seq: decisions and current state

Decisions made by Han for this analysis, with the reason for each. Steps refer to scripts in `pipeline/scripts/`
(single-cell-pipeline steps) and `analysis/scripts/` (study steps). All paths come from `study_config/config.sh`.

## Data

12 samples from 10x GEM-X On-Chip Multiplexing: D40 (differentiation 34), D15 and D22 (differentiation 37), each day
ctr rep1/rep2 and ISL1 KO rep1/rep2 from one OCM reaction. Sample names and order: `study_config/samples.tsv`,
`study_config/rename_manifest.csv`.

## Per-sample QC (steps 1.0-1.3), same rules for all 12 samples

| Step | Rule | Why |
|---|---|---|
| 1.0 nuclear fraction | keep cells with nuclear fraction >= 0.10 (intronic / (intronic + exonic) reads, Cell Ranger BAM `RE` tag) | removes cytoplasmic fragments without a nucleus (NF ~0.04, little MALAT1 and MT, high ribo); a %ribo cutoff and removing clusters were rejected |
| 1.1 QC | 4 MAD (counts, genes) + at least 2,000 genes + %MT <= 10 | the deeper second sequencing batch adds a low-gene peak (~1,000-1,800 genes) in D15/D22 controls that MAD alone keeps; 3 MAD was rejected (still uneven, cuts clean D40 cells); genes, not counts, because depth differs by day; no upper cutoff |
| 1.2 SoupX | autoEstCont with the pipeline's automatic tfidfMin fallback; never a fixed rho, never one genotype only | D40 KO_1/KO_2 fell back to tfidfMin 0.8 with 2-3 estimates (rho 1-2%, like controls): report as a caveat |
| 1.3 doublets | remove scDblFinder calls only; Scrublet not run | Scrublet's automatic threshold fails on these channels (unimodal simulated scores) |
| no cluster removal | filtering is per cell only | whole-cluster removal was rejected |

Cells kept after 1.3: 46,477 (per-sample table from `analysis/scripts/1.4_qc_cell_accounting.py`).

## Per-day processing (steps 7.1, 7.1c, 7.1b)

- Each day on its own: its 4 samples merged, **no batch correction** (one OCM reaction per day; sample is confounded
  with genotype; replicates mix well).
- 7.1: damaged cells removed at the merge, after doublet removal: **nuclear fraction > 0.6** (high-NF edge streaks,
  mostly controls). Kept: D15 16,741, D22 15,874, D40 12,641. Cell-cycle scores (Seurat `cc.genes`) on log1p.
- **Normalization for the graph: SCT B** (chosen 2026-11-18 after the 7.1c comparison of log1p, SCT, SCT + %MT
  regression and SCT B): SCTransform v2 (glmGamPoi) on the merged day as one sample, S - G2M regressed, 3,000
  variable features with the mitochondrial (MT-) genes removed (they acted as one %MT signal that drove PC1 in D15),
  topped up from the next-ranked genes. 30 PCs, 15 neighbours, UMAP, seed 1234.
- Gene plots, markers and DE use unregressed log1p values.
- 7.1b: Leiden (igraph, 2 iterations, seed 1234) at resolutions 0.2-1.2 on the SCT B graph; one report per day
  (`day_analysis/<day>/clustering/<day>_clusters.html`). Resolution and annotation are decided day by day.

## Compartment subclustering (step 7.1d)

- D15 overview at **resolution 0.2** (7 clusters, 2026-11-18). Endocrine compartment = r0.2 clusters **0, 1, 5**
  (CHGA+/FEV+ endocrine; control-only ISL1+/SCG5+; DCC+/TOX+ endocrine progenitors), minus 290 non-endocrine
  carry-over cells (user, 2026-11-19; first-run r0.6 subcluster 9: FN1, VIM, AFP, FGB, low counts, cycling; barcodes
  in `study_config/annotation/D15_endocrine_excluded_cells.txt`, 7.1d exclusion option): **5,994 cells**.
- The compartment is re-processed on its own: SCT B refit on its cells only (`7.1c_day_sct.R` subset mode, same
  settings), 30 PCs, 15 neighbours, UMAP in 2D and 3D, Leiden 0.2-1.2. Markers vs the rest of the compartment.
  Output: `day_analysis/D15/subcluster_endocrine/` (report, `umap3d.html` interactive 3D UMAP).

## Current state (2026-11-20)

- 7.1b clustering on SCT B done for D15, D22 and D40.
- D15: overview resolution 0.2; endocrine subset (7.1d) rerun without the non-endocrine cells (2026-11-19), 11
  subclusters at r0.6. Names (user, 2026-11-19; `study_config/annotation/D15_endocrine_r0.6.tsv`, no groups):
  0 **pre_beta2**, 1 **ENP2**, 2 **ENP1**, 3 **ENP_PHOX2A+**, 4 **ENP_FEV_high**, 5 **pre_alpha**, 6 **pre_delta**,
  7 **pre_EC**, 8 **pre_beta1**, 9 **pre_alpha_KO**, 10 **KO_only**. The first run's draft labels were removed (the
  rerun renumbered the subclusters); each cell's previous label is kept on scratch
  (`day_analysis/D15/subcluster_endocrine_previous_labels.csv`) and shown in the 7.1e report next to the new
  subclusters. Report: `D15_endocrine_r0.6_report.html` (2D UMAP split by genotype with subcluster numbers and
  names at the centres). Not decided yet: whether to subcluster the other D15 compartments.
- D15 + D22 endocrine integration (step 7.2, user, 2026-11-19): D15 endocrine subset (5,994) + D22 endocrine cells
  (every 7.1h label except non-endocrine, 15,069). SCT B on the merged cells (S - G2M regressed, MT- genes out of the
  3,000 variable genes); batch = stage only, genotype not corrected; Harmony (PCs 1-30 of the merged model), RPCA
  (one SCT B model per stage; anchor genes = the merged model's genes present in both stage models, 2,996 of 3,000),
  scVI (counts of the same 3,000 genes, batch = stage, S - G2M covariate). No clustering. Report
  `day_analysis/integration_D15_D22_endocrine/D15_D22_endocrine_integration.html` (UMAPs by stage / genotype /
  labels, stage mixing, labels and genotype kept, nearest D22 label per D15 state, marker genes). D15 pre_delta
  lands next to D22 SC_delta (Harmony 74%, RPCA 84%, scVI 58% of its cells).
  - Method: **scVI** (user, 2026-11-19: less over-correction, keeps most of the D15 structure). Its numbers in the
    report (30 nearest neighbours): stage mixing D15 32.3% / D22 7.9% (Harmony 36.6 / 6.1, RPCA 22.6 / 11.1);
    same-label neighbours D15 68.3% / D22 90.9% (Harmony 77.5 / 95.1, RPCA 82.3 / 95.9); genotype kept as
    uncorrected. All three agree on pre_alpha -> SC_alpha, pre_beta2 -> SC_beta, pre_delta -> SC_delta,
    pre_EC -> SC_EC2; under scVI 29% of pre_alpha_KO sits next to SC_alpha (Harmony 2%, RPCA 2%).
- D15 trajectory (step 7.3, user, 2026-11-19): question = where the KO stops on the way to alpha / beta and which
  genes switch on right after (candidates for ISL1-controlled steps).
  - 7.3a Slingshot in PCA space (SCT B PCs 1-30 and 1-15, never on a UMAP; curves only drawn on the 2D / 3D UMAPs),
    start ENP1, ends fixed to pre_beta2, pre_alpha, pre_delta, pre_EC, ENP_PHOX2A+ (user: the free-ended tree put
    pre_delta after pre_beta2). Lineages (same for 30 and 15 PCs): ENP1 -> ENP2 -> ENP_FEV_high -> pre_beta1 ->
    pre_beta2 | -> pre_alpha_KO -> pre_alpha | -> pre_beta1 -> pre_EC | -> ENP_PHOX2A+ | -> KO_only;
    ENP1 -> ENP2 -> pre_delta. condiments: KO progression differs on every lineage (strongest alpha, beta, delta;
    weakest PHOX2A, p 0.03), and KO splits between lineages differently (p < 1e-58).
    Fix (2026-11-19): the embedded-curve files (`slingshot_pc{n}_curves_2d/3d.csv`) were empty, because embedCurves
    (slingshot 2.14) returns the curves without lineage names; curves are now named from the lineages and an empty
    table stops the job. 7.3a rerun: 6 lineages x 150 points in each file. Pseudotimes and weights unchanged (two
    reruns identical; the first run differed only at floating-point level, progression statistic 17.204174 vs
    17.204169), so 7.3e tradeSeq (started on the first run's `slingshot_pc30_cells.csv`) was not rerun. The
    differentiation test is classifier-based and varies between runs (PCs 1-30, all lineages: 0.700 / p 3e-66 first
    run, 0.692 / p 2e-61 rerun).
  - **Slingshot stays one shared trajectory for control and KO** (user, 2026-11-20): separate per-genotype
    trajectories would give pseudotimes that cannot be compared and would break the condiments and tradeSeq
    control-vs-KO tests; condiments topologyTest is the formal check if ever needed (not run).
  - 7.3b velocyto per D15 sample (spliced / unspliced, Cell Ranger BAMs) and 7.3c scVelo dynamical model with arrows
    on the 2D and 3D UMAPs: done. **Velocity runs backwards on D15**: latent time is highest in ENP1 (median 0.78) and
    lowest in pre_alpha_KO; Spearman with Slingshot pseudotime -0.43 / -0.82 (beta lineage, control / KO), -0.08 /
    -0.89 (alpha). Not followed up yet.
  - 7.3d CellRank 2 (own env `CELLRANK_ENV`; CellRank 2.0.7 pins Jinja2 < 3.1, single.cell has 3.1.5): kernels
    0.8 velocity + 0.2 connectivity and 0.8 Slingshot pseudotime + 0.2 connectivity; fixed terminal states for fates,
    data-driven macrostates (eigengap 4-12, retried with a lower upper bound when the chosen number splits a
    complex-conjugate eigenvalue pair). Runs (user, 2026-11-20): **all cells, control only and KO only**; the
    per-genotype runs rebuild the 30-neighbour graph within the genotype (the shared graph lets KO cells borrow
    control paths), keep the shared scVelo velocities and Slingshot pseudotime, and fix the ends the genotype has
    (control: pre_beta2, pre_alpha, pre_delta, pre_EC, ENP_PHOX2A+; KO: pre_EC, ENP_PHOX2A+, KO_only, pre_alpha_KO).
    Results:
    - The velocity kernel's only data-driven terminal state is ENP1 in every run (the backwards velocity above), so
      terminal states are read from the pseudotime kernel. KO (pseudotime kernel): 4 macrostates, all terminal:
      ENP_PHOX2A+, KO_only, pre_EC, pre_alpha_KO; pre_beta1 is not an end in KO (KO pre_beta1 -> pre_EC 0.70).
    - Both kernels: KO beta fate is 2-3x lower at ENP_FEV_high and pre_beta1 and goes to pre_EC (and PHOX2A+);
      alpha fate does not differ between genotypes in the early cells. The kernels disagree on where ENP1 / ENP2
      head (velocity: pre_alpha; pseudotime: PHOX2A+ / pre_EC).
    - Drivers (Pearson r of log1p expression with a fate probability in the early cells ENP1, ENP2, ENP_FEV_high;
      Fisher z, BH; the same statistic as CellRank's compute_lineage_drivers). Control program (control-only run,
      r >= 0.2 and FDR < 0.05 in both kernels): pre_alpha 212 genes (FLRT2, ARX, DEPP1, RGS4, ...), pre_beta2 62
      (INS, CRYBA2, HADH, ACVR1C, ...), pre_EC 222, pre_delta 23, ENP_PHOX2A+ 1 (kernels disagree). KO-only run:
      KO_only drivers (TRPM3, CNTNAP2, TOX, DCC, SEMA6A, ...) overlap the control pre_delta drivers; NKX6-1 is a
      KO driver toward PHOX2A+ (velocity kernel); pre_alpha_KO drivers are mostly cilia genes plus ARX.
  - 7.3e tradeSeq on all genes detected in >= 10 cells (control vs KO condition tests per lineage, early DE at
    branches, gene modules by curve shape): running.
  - 7.3f **lost in KO** (user, 2026-11-20: CellRank drivers plus a lost-in-KO comparison; tradeSeq covers the
    curve-based genes). `cellrank_lost_in_KO.csv`, `_profile.csv`, `cellrank_on_after_KO_stop.csv`.
    - Candidate genes = the control program above (both kernels; one-kernel genes flagged).
    - KO vs control by expression, two matchings (user, 2026-11-20: the label match with the proposed thresholds,
      then "add the pseudotime-matched comparison"): by label (beta: KO vs
      control pre_beta1; alpha: KO pre_alpha_KO vs control pre_alpha; both: ENP_FEV_high) and by **Slingshot
      pseudotime** (user, 2026-11-20: Slingshot only, not velocity, which runs backwards): cells with curve weight
      >= 0.5 on the beta / alpha lineage inside the stop window = KO's 75th-95th percentile of pseudotime (beta
      143-151, alpha 140-158). Lost = log2FC <= -0.5, Wilcoxon FDR < 0.05, both KO replicates below both control
      replicates (pseudobulk); gained = the same upward.
    - **Genes control turns on after the KO stop** (user, 2026-11-20): control cells in the first half of the path past
      the stop window vs control cells in the window (log2FC >= 0.5, FDR < 0.05, up in each control replicate);
      missing in KO = KO cells in the window below the control cells after it (same rules).
    - Results. At the stop window control cells are mostly ENP_FEV_high while many KO cells are already pre_EC (beta
      lineage) or pre_alpha_KO (alpha): KO branches off at commitment. Beta lost in both matchings: PDX1, ACVR1C,
      FOXO1, TENM2, GRIN3A, SMC2, SLC7A2, AC109466.1. On after the stop and missing in KO: beta 244 genes (TFs
      TWIST1, HDAC9, LMO2, ISL1; program genes GRIN3A, SMC2, ITGA1, SLC17A6, ABCC8, ACVR1C, ALDH1A1, INS, ...),
      alpha 809 (the after window is nearly all control pre_alpha; TFs ARX, ETV1, PAX6, HHEX, MEIS1, ESRRG, IRX3,
      IRX2, LMO2, GLIS3, FOXO1, ID4, TLE1, ISL1). Gained in KO pre_alpha_KO: NKX6-1 (9% -> 58% of cells), PAX4.
    - Caveats: the alpha pseudotime match compares KO cells that have moved on with uncommitted control cells, so
      for alpha the label match and the after-stop lists are the useful ones; 53 genes are gained in all four label
      comparisons and none lost in all four (e.g. TSPYL5, RBM46, CTSF, CDKN1C; H19 / DLK1 lower early): KO-line-wide
      differences, possibly clonal imprinting, to be excluded or checked against the bulk KO; liver-type genes (ALB,
      ITIH2, APOA2) in the beta after-stop list.
  - D15 endocrine report (7.1e, user, 2026-11-20): the colour menu has the 7.3f gene groups next to the endocrine
    markers (beta on after the stop and missing in KO: TFs + program genes; beta lost at the stop; alpha on after the
    stop: TFs; alpha program genes, top 40 by expression after the stop; gained in KO: NKX6-1, PAX4, GCG, CCK,
    NEUROD4), read from the 7.3f tables, and a gene search box.
  - 7.3g **KO vs control along the beta, alpha and EC lineages** (user, 2026-11-21: use the trajectory version of the
    control vs KO comparison and add the pre_EC lineage; plan accepted as proposed). KO cells do not stop on the EC
    lineage (Slingshot Lineage3, ENP1 -> ENP2 -> ENP_FEV_high -> pre_beta1 -> pre_EC), so each lineage is cut into 8
    windows of equal cell numbers (curve weight >= 0.5), tested when >= 30 cells of each genotype; KO vs control per
    window with the 7.3f rules (log2FC 0.5, Wilcoxon FDR < 0.05, both KO replicates beyond both control replicates);
    window part = shared path / branch (>= half of the cells also on another of the three lineages). Programs: control
    beta / alpha / EC (control-only CellRank runs) and KO EC (KO-only runs). `ko_windows.csv`,
    `ko_vs_ctr_windows.csv.gz`, `ec_programs.csv`, `ko_lineage_profile.csv`.
    - 7.3f gained an EC-program flag (EC_program, EC_r / FDR per kernel) on every row (user, 2026-11-21); its calls were
      checked identical to the earlier run.
    - Results. Beta: KO never reaches the beta branch (window 8, 0 KO cells), so every tested beta window is shared
      path; the beta-program score is lower in KO from ENP2 on (e.g. 0.79 vs 0.95, 0.97 vs 1.11) and the pre_beta2
      fate 2-6x lower while the pre_EC fate is higher (0.92 vs 0.68 in window 7). Alpha: KO reaches window 7
      (pre_alpha_KO) but not 8 (pre_alpha); in window 7 the KO alpha fate is 0.35 vs 0.64 and EC fate 0.44 vs 0.18.
      EC: the control EC-program score is about the same in KO (slightly lower at the end); control EC-program genes
      are more often lower than higher in KO on the EC branch, KO-only EC-program genes more often higher (e.g.
      NEUROD4, MYCN, GAL, VEGFA): KO reaches an EC-like state with its own program rather than an exaggerated control EC
      state. Beta-program genes on the EC lineage are only ever lower in KO (14-19 per window), so KO cells on the EC
      branch do not carry more beta program than control cells there. EC programs (any kernel): 291 genes in both
      genotypes, 202 control only, 164 KO only (both kernels: 112 / 43 / 27).
    - **KO-line-wide genes**: 42 genes are called in >= 15 of the 22 tested windows, from ENP1 on (higher in KO: TSPYL5,
      CTSF, RBM46, LINC00884, BHMT2, ...; lower: GRIN3A, SMC2, SLC7A2, AC109466.1, FGB, FGG, AFP, H19, CXCL11, ...).
      4 of the 8 genes 7.3f called beta-lost in both matchings (GRIN3A, SMC2, SLC7A2, AC109466.1) are among them, so
      they are line-wide differences rather than beta-specific (PDX1, ACVR1C, FOXO1, TENM2 are not). Window
      differences also partly reflect different label mixes of KO and control cells at the same pseudotime.
  - 7.3h **gene patterns** (user, 2026-11-21): unbiased clustering only (no rule-based classes), binned curves now
    (tradeSeq smoothers can replace them later), genes = **only the KO vs control trajectory results**: 7.3f stop calls
    (pseudotime match), 7.3f "on after stop, missing in KO", 7.3g window calls (2,656 genes; CellRank programs are labels,
    not a source; the 42 line-wide genes stay in, flagged). Curves: 10 equal-width pseudotime bins per lineage,
    curve-weighted mean log1p, 3-bin moving average; bins with < 15 cells of a genotype left out (KO never reaches beta
    / alpha bins 8-10). Features z-scored per gene across the 6 curves; k-means, k 4-20 by silhouette.
    - First pass: silhouette is low at every k (0.23 at k = 4 down to 0.13), so the patterns are a continuum; k = 4
      chosen (stable across seeds, ARI 0.995). C1 802 genes: rise late on the control beta and alpha lineages, not on
      EC, KO never reaches them (GCG, SCG5, SCG3, GNAS, TTR, TENM2); C2 602: early progenitor genes falling on every
      lineage (SOX4, DLK1, CDKN1C, TOX, DCC); C3 826: rise toward the EC end and mid beta / alpha, KO the same on EC
      (CHGA, INS, CRYBA2); C4 426: transient mid-path peak that control turns off late on beta / alpha while KO stays
      high (RUNX1T1, NEGR1, SCG2, C9orf16).
    - Finer cuts (user, 2026-11-21: try k 8-10): `gene_patterns_k{8,9,10}.csv` (+ clusters, heatmaps, curves, and
      `gene_pattern_split_k4_to_k*.csv`); k-means agreement across seeds 0.97 / 0.97 / 0.95. At k = 10: C1 370 rise late
      toward beta (less alpha, little EC), KO flat (TTR, APOA2, ISL1, SCGN, CALB2, PTPRN2, TENM2); C2 148 lower in KO on
      all three lineages from the start (13 line-wide; VIM, AFP, ID2, S100A11); C3 343 early progenitor genes falling in
      both genotypes (SOX4, DLK1, CDKN1C, TOX, RFX6); C4 216 dip then late rise on alpha in control, KO no recovery
      (CPE, RGS4, NPAS3, PTPRG); C5 324 alpha-only rise in control, KO starts but stops (GCG, SCG5, SCG3, GNAS, SLC30A8,
      ALDH1A1, GAST); C6 130 high early in KO only (LDHA, SLC2A1, ZEB2, JPT1); C7 427 rise toward beta and the EC end in
      both genotypes, KO follows (INS, CRYBA2, STMN2, APOA1); C8 255 mid-path rise that control turns off late on beta /
      alpha while KO keeps rising, higher in KO on EC (CHGA, LMX1A, MNX1, FEV, NR0B1); C9 241 transient mid-path peak in
      both, turned off late (RUNX1T1, SCG2, cilia genes); C10 202 rise on the EC branch, stronger in KO and up in KO at
      the beta / alpha stop (PTHLH, NETO1, NTM, KCNJ3, BHMT, DUSP4) = "switched to EC in KO".
    - **Working classification: k = 10** (user, 2026-11-21): `gene_patterns_k10.csv`, `gene_pattern_clusters_k10.csv`,
      `gene_patterns_{curves,heatmap}_k10.png` (rerun: `sbatch analysis/scripts/7.3h_gene_patterns.sb 4,8,9,10`; the
      plain-named files are the silhouette-best k = 4). To recheck with tradeSeq (7.3e) curves when they exist.
  - Next: cross the 7.3f lists with ISL1 CUT&RUN and bulk KO; comparison report (Slingshot, velocity, CellRank,
    KO drop-off, tradeSeq).
- D15 + D22 trajectory (step 7.4, user, 2026-11-19): the 7.3 set (Slingshot, velocity, CellRank, tradeSeq) on the
  7.2 integration, scVI latent as the space, scVI UMAP as the primary UMAP; start D15 ENP1, ends fixed to D22 labels;
  tradeSeq on all genes detected in >= 10 cells, conditions = genotype, stage (D15 / D22) as a fixed covariate.
  - First run on all 21,063 cells (ends SC_beta / SC_alpha / SC_delta / SC_EC1 / SC_EC2): 8 lineages, with orderings
    the D15 run did not have (pre_delta after pre_beta2 on the way to SC_delta; the beta and delta lineages through
    pre_EC before pre_beta1; free ends D15 KO_only, D15 ENP_PHOX2A+, D22 GHRL+). tradeSeq, scVelo and CellRank were
    cancelled before they started.
  - **Focus on alpha, beta and EC only** (user, 2026-11-19): drop D15 pre_delta, KO_only, ENP_PHOX2A+ and D22
    SC_delta, ENP, GHRL+ (labels kept: `study_config/annotation/integration_subset_alpha_beta_EC.txt`), and
    **retrain** the integration on the kept cells so the UMAP is made from them alone: the whole 7.2a-d chain
    (SCT B with new variable genes, Harmony, RPCA, scVI, report) with `INT_SUBSET=alpha_beta_EC`, outputs in
    `day_analysis/integration_D15_D22_endocrine_alpha_beta_EC/`; the full integration is left as it was. The new
    scVI is to be checked in its 7.2d report before 7.4 runs on it (`bash analysis/scripts/7.4_submit.sh
    alpha_beta_EC <velocyto jobs> slingshot` first, so the lineages are seen before tradeSeq).
  - Velocyto D22 (7.3b with a day argument) writes the cellsorted BAMs next to the Cell Ranger BAMs; they are kept
    (user, 2026-11-19).
  - alpha_beta_EC subset (2026-11-20): the user checked the 7.2d report of the retrained scVI ("looks good") and
    asked for Slingshot on it. Slingshot: 4 lineages, start D15 ENP1, all through ENP2 -> ENP_FEV_high:
    1 -> pre_EC -> pre_beta1 -> pre_beta2 -> SC_beta; 2 -> pre_alpha_KO -> pre_alpha -> SC_alpha;
    3 -> pre_EC -> SC_EC2; 4 -> SC_EC1. The beta lineage still passes pre_EC before pre_beta1 (in this scVI 82% of
    pre_beta1 cells sit next to SC_EC2); the user accepted the lineages and asked to submit tradeSeq, scVelo and
    CellRank on them (2026-11-20). condiments: progression All 16.2 (p 1e-59; beta 0.76, alpha 0.68, EC2 0.22,
    EC1 0.26); differentiation All 0.78 (p 0), EC2 vs EC1 0.52 (p 0.01).
  - CellRank (7.4d), as 7.3d (user, 2026-11-20): all cells plus control-only and KO-only runs (graph rebuilt within
    the genotype on the scVI latent); per-genotype fixed ends = the Slingshot ends with >= 50 cells of the genotype,
    KO also D15 pre_alpha_KO (alpha_beta_EC: control SC_beta, SC_alpha, SC_EC1, SC_EC2; KO SC_EC1, SC_EC2,
    D15 pre_alpha_KO); same macrostate retry. Queued behind 7.4c scVelo.
  - Job names: every D15 + D22 job starts with `isl1_int_` (D15-only jobs `isl1_`), user 2026-11-19.
  - Environment fixes (CELLRANK_ENV: scvelo 0.3.4, cellrank 2.0.7, numpy 2.4.6), 2026-11-20, in 7.3c/d and 7.4c/d:
    looms read with anndata.read_loom (no scv.read); filter_and_normalize only filters + normalises, so log1p and
    the top 2,000 genes are done with scanpy; CellRank's VelocityKernel passes x= / y= to
    np.testing.assert_array_equal (wrapped in the scripts); fate probabilities with use_petsc=False (the env's
    petsc4py has no Mat.getVecs). The env itself is unchanged.
    2026-11-21: the D15 + D22 scVelo (7.4c) failed at the 30-neighbour graph after 7.5 h (scanpy 1.11.5's default
    pynndescent path calls csr_array.getnnz, gone in scipy 1.17.1; only above a few thousand cells, so the D15
    runs were not affected). 7.4c and the per-genotype graphs of 7.4d now use exact neighbours
    (`transformer="sklearn"`); the D15 scripts are unchanged. 7.4c prints a timestamp per step to find
    where the 7.5 h went (user, 2026-11-21: fix and resubmit scVelo + CellRank).
- **Competency question** (user, 2026-11-20). PLAN ONLY -- the user is still reviewing it; nothing is run.
  - Experimental findings (user): the progenitors competent for alpha and for beta differ, and competency changes
    over time: **alpha competency lasts until ~D18 and is then lost; beta competency lasts until ~D26**. In the
    **ISL1 KO, alpha-competent cells turn on NKX6-1** (not seen in control alpha-competent cells) but keep alpha
    competency until D18.
  - Goal: use the D15 + D22 (alpha_beta_EC) trajectory to define alpha-competent and beta-competent progenitors and
    say why their competency differs. D15 is inside both windows; D22 is after the alpha window and inside the beta
    window, so D22 progenitors should be beta-competent only. D15 and D22 only, so the D18 boundary itself cannot be
    placed. Fate probabilities are model predictions; the gene lists are hypotheses for experiments.
  - **D22 progenitors by expression only** (user, 2026-11-20): the D22 "ENP" cells (135, 128 control) stay out of the
    integration and trajectory; they enter only the competency step (program scores, classifier). The running
    tradeSeq / scVelo / CellRank stay valid.
  - Proposed steps (draft script `analysis/scripts/7.4f_competency.py`, not run):
    1. groups among early D15 progenitors (ENP1, ENP2, ENP_FEV_high), per genotype, from CellRank fates only (both
       kernels must agree): alpha-biased, beta-biased, mixed, EC-directed (proposed: alpha bias = P(alpha) /
       (P(alpha) + P(beta)) >= 2/3, or <= 1/3, and P(alpha) + P(beta) >= 0.2). Never from markers: in KO, NKX6-1
       does not mean beta competence.
    2. checks against the experiments: NKX6-1 on in KO alpha-biased cells, off in control alpha-biased cells; D22 ENP
       should look beta-biased and lack the alpha program.
    3. where / when: D15 state, cell-cycle phase, latent time, Slingshot time per group.
    4. why: per-sample pseudobulk edgeR -- alpha- vs beta-biased (control, paired by sample; and in KO), KO vs
       control within alpha-biased and within beta-biased cells. Alpha / beta programs from the control test; core
       alpha genes = also alpha-up in KO (separates the alpha program from the KO-only NKX6-1 response). Plus
       tradeSeq early-DE (alpha vs beta lineage near the branch) and CellRank drivers.
    5. time: alpha-program genes that drop in D22 ENP (candidates for losing alpha competency); beta-program genes
       kept in D22 ENP (candidates for the longer beta window). D15 vs D22 is also a batch difference: descriptive.
    6. later: TF activity (decoupler + CollecTRI, not installed), cross with ISL1 CUT&RUN and the bulk KO, report.
  - Open: the user to check the plan (thresholds, groups, contrasts) before anything runs.
- D22 **annotated** (user, 2026-11-18). Rules: `study_config/annotation/D22_annotation.tsv`; per-cell labels from
  step 7.1h (`day_analysis/D22/labels/`: labels.csv, composition.csv, umap_labels.png).
  - Base: overview resolution **0.2**: 0 **SC_alpha**, 1 **SC_beta**, 2 **non-endocrine**, 3 **SC_EC1**, 4 **SC_EC2**.
  - Cluster 2 (805 cells, mostly control, about half in S phase) mixes hepatic-like cells (AFP, FGB/FGG, GPC3,
    APOA2; at the tip) and pancreatic progenitor / duct-like cells (SOX9, ONECUT1, PDX1, some NKX6-1; PTF1A off);
    it stays one cluster up to r0.8 apart from 18 liver-like cells, so it is not labelled as liver cells.
  - Override: r1.0 cluster 1 (346 cells, 92% control) holds the delta-like bridge and an endocrine precursor stem
    (7.1f gene check at r1.0). It is split by Leiden restricted to its cells on the SCT B graph (7.1g) at
    sub-resolution **0.1** (same split at 0.2): subcluster 1 = **SC_delta** (143 cells, all control; SST 97%,
    HHEX 36% with log1p > 1; INS also high, so partly SST / INS co-expressing), subcluster 0 = **ENP** (203 cells;
    GHRL, ARX, NEUROG3, NEUROD1, CDKN1C). These cells take these labels over their r0.2 label.
  - GHRL+ (user, 2026-11-19): the ENP stem is split further at sub-resolution **0.3** (nested in the 0.1 split):
    subcluster 1 = **GHRL+** (68 cells, 46 control / 22 KO; GHRL 96%, ARX 72%, NEUROG3 9% with log1p > 1),
    subcluster 0 stays **ENP** (135 cells; NEUROG3 24%, GHRL 30%). Reason: D40 r0.4 cluster 9 (166 cells, 96
    control / 70 KO; GHRL 37%, SST 21%, no NEUROG3) matches the GHRL+ part, so D22 has both ENP and GHRL+ cells and
    D40 only GHRL+ cells. Named GHRL+ (not epsilon), describing the marker only.
  - Not labelled separately: the NEUROG3+ cells on the KO side (mostly r1.0 cluster 4, inside SC_EC1 / SC_EC2);
    the small SST + GCG + INS polyhormonal end of the delta arm (61 cells at sub-resolution 0.5).
  - Caveat: ambient INS is strong on D22 (INS log1p > 1 in 60-100% of cells in the KO clusters). GCG shows the
    same low-level background in KO cells (detected in 64% of KO cells, at low values, vs 28% of control cells).
  - Plots split by genotype (step 7.1i, 2026-11-19): `day_analysis/D22/labels/` umap_labels_split.png (labels: all /
    control / KO) and feature_split_markers.png (2-4 key markers per label, control | KO on one colour scale);
    D22_labels.html. Feature-plot style (user, 2026-11-19, chosen from two versions): all cells drawn, the other
    genotype in very light grey, lighter than the shown genotype's non-expressing cells.
- D40 (in progress, 2026-11-19). Overview resolution **0.2** (7 clusters). Labels (user;
  `study_config/annotation/D40_annotation.tsv`): 0 **proliferating** (MKI67 / TOP2A 75%, CHGA+ with INS / GCG;
  cycling endocrine cells, mostly control), 1 **SC_EC2**, 2 **SC_EC1**, 3 **SC_beta**, 4 **SC_alpha**. The EC names
  follow D22: D40 cluster 1 shares 19 of its top 50 markers with D22 SC_EC2 (TPH1, LMX1A, AFF3, HDAC9, SYT13) and 1
  with SC_EC1; cluster 2 shares 12 with D22 SC_EC1 (PLCL1, RUNX1T1, PRKCH, PDE4B) and 1 with SC_EC2.
  - Cluster overview report: step 7.1j (`day_analysis/D40/clustering/overview/`); gene checks at r0.2 / r0.4 (7.1f).
  - Not labelled yet: r0.2 cluster 5 (516 cells, 92% control; SOX2 / NES / VIM progenitor-like + SOX11 / STMN2 /
    ELAVL4 neuron-like: off-target neural) and cluster 6 (285 cells, 96% control; ONECUT1 / GATA6 / KRT19 / SOX9
    epithelium with hepatic-like features + a COL1A2 / COL3A1 fibroblast-like group). Next: split each with 7.1g.
  - Cluster 6 labelled (user, 2026-11-19): 7.1g split at sub-resolution **0.2** (4 groups; at 0.1 the hepatic and
    progenitor groups are still merged): 0 **acinar** (35 cells; PTF1A detected in 71%, CPA2 86% and CTRB2 57% with
    log1p > 1, PRSS1, CPA1, RBPJL; NR5A2, DLK1; 40% cycling: pre-acinar / tip-like), 1 **hepatic** (68; AFP 49%, FGB
    52%, HNF4A 54%; also CDX2 / CDH17, so hepatic- and gut-like), 2 **pancreatic progenitor** (104; ONECUT1 96%, SOX9
    57%, HNF1B, KRT19, NKX6-1; PDX1 detected in 73% but low), 3 **fibroblast** (78; COL3A1 96%, PDGFRA 85%, COL1A1
    83%). User's choice for group 1: "hepatic, or matching the D22 name"; D22 has no hepatic-only label (its
    non-endocrine cluster mixes hepatic-like and progenitor cells), so "hepatic" is used.
  - Cluster 5 labelled **neural** as one group (user, 2026-11-19; base rule, the 7.1g split is not used for labels):
    at sub-resolution 0.1 it holds cycling early neural progenitors (102), immature neurons (174; DCX, SOX11, STMN2,
    ELAVL4) and radial-glia-like cells (240; SOX2, VIM, FABP7, SLC1A3, NFIA); off-target, 92% control.
  - Per-cell labels (7.1h, 2026-11-19): every r0.2 cluster now has a label, so `day_analysis/D40/labels/` covers all
    12,641 cells (10 labels). The r0.4 cluster 9 groups keep their r0.2 labels until their own labels are decided.
  - No separate delta cluster and no NEUROG3+ progenitors at D40 (NEUROG3 log1p > 1 in at most 0.1% of any cluster).
    SST+ cells are scattered in SC_alpha and concentrated in r0.4 cluster 9 (166 cells, between SC_alpha and the KO
    lobe). 7.1g split of r0.4 cluster 9 (same groups at sub-resolution 0.2, 0.3, 0.5): 0 SST+ delta-like (52 cells,
    50 control; SST 60%, HHEX 15%, PBX3, MEIS1), 1 GHRL+ (36; 18 control / 18 KO; GHRL 89%), 2 GC / PDK3 / CASR /
    PAX6 (35; no SST / GHRL, alpha-lineage-like?), 3 KO EC-like (43, 41 KO; SLC18A1, TPH1, MNX1).
  - r0.4 cluster 9 labelled (user, 2026-11-19) at sub-resolution **0.1** (2 groups): subcluster 0 = **GHRL+** (132
    cells, 70 control / 62 KO; = the sub-resolution 0.2 groups 0, 1 and 3 plus 1 cell; GHRL 46%, SST 27% with
    log1p > 1), subcluster 1 = **SC_delta** (34 cells, 26 control / 8 KO; = sub-resolution 0.2 group 2: GC, PDK3,
    CASR, PAX6, TOX3). Noted when the labels were given: SST and HHEX are not detected above log1p 1 in the SC_delta
    group, and the GHRL+ group also holds the SST+ delta-like cells and 41 KO EC-like cells; the user kept 0.1.
    These cells take these labels over their r0.2 label.
- D22 vs D40 label match (step 7.1k, 2026-11-19; `day_analysis/compare_D22_D40/`): days clustered separately, so
  (1) shared top-50 Wilcoxon markers per label pair and (2) Pearson r of gene-centred label mean profiles on 1,745
  genes highly variable in both days (all cells; also control-only and KO-only).
  - Match (same label is the best partner by both measures): SC_alpha (26 shared markers, r 0.71), SC_beta (24,
    0.72), SC_EC2 (19, 0.68), SC_EC1 (13, 0.56; D22 SC_EC2 also has r 0.53 with D40 SC_EC1), GHRL+ (20, 0.51).
  - SC_delta does not match: D22 SC_delta vs D40 SC_delta 5 shared markers, r 0.13; both are closest to the other
    day's SC_alpha (D22 SC_delta -> D40 SC_alpha 13, r 0.35; D40 SC_delta -> D22 SC_alpha 8, r 0.29).
  - D22 ENP has no clear D40 partner (best: GHRL+, 9 shared, r 0.21). D40 acinar, hepatic, pancreatic progenitor
    and fibroblast all map to D22 non-endocrine (hepatic best, 14 / r 0.60). D40 proliferating is closest to D22
    SC_alpha (r 0.37); D40 neural has no D22 partner.
  - Independent check of 7.1k (2026-11-19, one-off script, not a pipeline step):
    - Centring choice: 7.1k centres each gene on the mean of its own day's labels, and D40 has more
      non-endocrine labels than D22. Re-centring both days on the 6 shared labels gives the same matches (same-label
      r: SC_alpha 0.69, SC_beta 0.73, SC_EC1 0.65, SC_EC2 0.72, GHRL+ 0.63; Spearman 0.57-0.67) and SC_delta still
      does not match (r 0.10).
    - Why SC_delta does not match: the D40 cells that match D22 SC_delta are the SST+ group of r0.4 cluster 9 at
      sub-resolution 0.2 (group 0, 52 cells, 50 control): r 0.51-0.54 and 14 shared top-50 markers with D22
      SC_delta, its best partner by both. Under the sub-resolution 0.1 labels these cells are in D40 **GHRL+**. The
      cells labelled D40 SC_delta are the GC / PAX6 group (sub-resolution 0.2 group 2), closest to D22 SC_alpha
      (r 0.30, 7 markers). Of the other sub-resolution 0.2 groups, group 1 (GHRL+, 36 cells) matches D22 GHRL+
      (r 0.69-0.73, 20 markers) and group 3 (KO EC-like, 43 cells, in D40 GHRL+) has no clear partner (D22 SC_EC2 /
      GHRL+, r 0.29-0.40).
    - Decision (user, 2026-11-19): keep the sub-resolution 0.1 labels as they are, knowing the above.

## Open items

- Step 1.4 does not yet count the NF > 0.6 cells removed at the per-day merge.
- Final report caveats: D40 KO SoupX fallback; scDblFinder run-to-run variation; the 2,000-gene rule; the NF rules;
  the strong ctr / KO separation on D22 and lower cycling in D22 / D40 KO.
- October-only output `HZ_Diff34_D40_scRNA`: decide whether to keep.
