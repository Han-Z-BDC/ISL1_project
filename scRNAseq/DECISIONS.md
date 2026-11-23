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

## Current state (2026-11-23)

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
    branches, gene modules by curve shape): running (job 36009099, 23,279 genes x 6 lineages, 32 cores; past 2 days of
    its 3-day limit on 2026-11-21, progress bar silent in parallel mode). **Cancelled** (user, 2026-11-21) together with
    the D15 + D22 fits 36045487 (~45% after 35 h) and combine 36045488: parallel fitting is the slow part and both would
    have stopped at the rowData bug before the tests; to be rerun serially in chunks (design pending).
    - Timing test (user, 2026-11-21; `7.3e_tradeseq_timing.R`, same inputs): a 1-gene matrix crashes fitGAM ("invalid
      row.names length"), so the smallest test is 2 genes; one core, 6 lineages: 2 genes 1.5 s / gene (ISL1, PDX1),
      10 genes 3.1 s / gene, 50 genes 4.6 s / gene, 3.3 GB memory. All 23,279 genes would be ~30 core-hours, so the
      2-day run is slower than the per-gene cost explains.
  - 7.3j **alpha / beta / EC driver table** (user, 2026-11-21: alpha, beta and EC drivers, then which genes differ by
    genotype along the three trajectories). pre_EC is a terminal state in every 7.3d run; the all-cell runs computed
    drivers for pre_alpha and pre_beta2 only (set on 2026-11-20), the control-only runs for every fate, so the table
    uses the control-only runs (as 7.3f): r >= 0.2, FDR < 0.05 in the control early cells, program both / one kernel.
    No CellRank rerun. Drivers: alpha 212 both + 442 one kernel, beta 62 + 359, EC 222 + 271; 1,299 genes, 269 of them
    drive both beta and EC, none alpha with beta or EC. `drivers_alpha_beta_EC.csv`.
  - 7.3k **subset tradeSeq** (user, 2026-11-21: instead of the full genome, the drivers and the KO vs control genes;
    chunked). Genes = 7.3j drivers (either kernel, any of the three fates) + the 2,656 KO vs control trajectory genes
    of 7.3h + the 7.3e endocrine markers = 3,383 (`tradeseq_subset_genes.csv`); lineages = beta (Lineage1), alpha
    (Lineage2), EC (Lineage3) only, cells with no weight on the three dropped; model and tests as 7.3e; FDR over the
    subset genes. 16 one-core chunks (acpu, 16G) + combine (`7.3k_tradeseq_subset_submit.sh`). Outputs
    `tradeseq_subset_{association,condition,earlyDE,condition_wide,smooth,modules}.csv` (lineage names, driver calls
    and 7.3h cluster on every row).
    - Run 2026-11-21: 5-core MulticoreParam fits took 28.5 min per 212-gene chunk (slower than one core, so the fit
      step is now serial); the first pass stopped after saving the models at `rowData` (SingleCellExperiment not
      loaded), the tests were rerun from the saved models (`reuse`). **The same `rowData` line is in 7.3e and 7.4e
      fit**: those runs will save their models and stop before the tests (not edited while running).
    - Result: 3,382 genes; condition test FDR < 0.05 on any lineage: **58 genes** (alpha 16 only, EC 15 only, beta 2
      only, all three 13, alpha + EC 6, beta + alpha 3, beta + EC 3); 54 of them are in the 7.3h set, 20 are drivers.
      Far fewer than the 7.3g window tests: the 1.5-fold threshold applies over the whole curve, and KO has almost no
      cells on the late beta / alpha segments, where its curve is extrapolated.
  - 7.3l **replicate pseudobulk KO vs control, all genes** (user, 2026-11-21: a simple but comprehensive test; the
    7.3h set misses timing / shape shifts, non-driver genes at the stop and the KO-diverted states, and the both-
    replicates rule is strict). Summed SoupX counts per replicate x position, edgeR quasi-likelihood (filterByExpr,
    TMM, ~0 + group, robust), samples with < 10 cells dropped; call = FDR < 0.05 and |log2FC| >= 0.5.
    - 1 curve: per lineage (beta, alpha, EC), the 8 equal-count windows of 7.3g, windows with all four replicates
      (beta and alpha w8 untested: no KO cells); any (joint test of KO - control in every window), level (mean
      difference), shape (difference changes along the lineage), per-window log2FC / FDR. Called: beta 884, alpha
      1,091, EC 1,269 (1,722 genes) by the joint test; mostly level differences, shape in ~10%. Added (user,
      2026-11-21): a gene is also called when any single window passes (FDR < 0.05, |log2FC| >= 0.5), since the joint
      test dilutes narrow differences: +111 beta, +288 alpha, +155 EC (pattern "window only"); 2,125 genes.
    - 2 past the stop (beta, alpha; 7.3f windows): missing in KO (control switches on right after the stop, KO at the
      stop lacks it) beta 210 / alpha 979; kept on in KO (control switches off, KO keeps it) beta 92 / alpha 655
      (NKX6-1, NEUROG3, FEV, LMX1A, MNX1 kept on in KO on alpha; INS, ISL1, ARX, GCG, SST missing in KO).
    - 3 KO states vs nearest control state: KO_only vs ctr ENP_FEV_high 962 up / 600 down; pre_alpha_KO vs ctr
      pre_alpha 1,395 / 1,629; ENP_PHOX2A+ KO vs ctr 59 / 55. KO-state-only genes (1,876) are a separate category.
    - Re-clustering (7.3h with GENESET=pseudobulk, same curves, k = 10; silhouette-best k = 4): with the single-window
      rule 3,089 genes (curve or past-stop calls), 993 of them new; old genes kept: window-only 1,042 of 1,532,
      after-stop ~all (before the rule: 2,820 genes, ARI old vs new k = 10 on shared genes 0.61). Old genes not
      kept (before the rule, 715): 385 not significant once replicate variance is modelled, 150 single-window only
      (now recovered), 86 too few counts to test, 59 old stop / after-stop only, 35 effect < 0.5; mostly genes
      called in one old window, at about half the expression of the kept genes.
    - (user, 2026-11-21) the 7.3m tradeSeq calls are added to the set: 373 of the 402 were already in it; the 29 called
      only by tradeSeq (all also KO-state calls; TFs MAFB, JAZF1) are tagged "tS only". 3,118 genes; k = 10 clusters
      stable (ARI 0.94 vs before, 97% of genes keep their cluster); the 29 spread over 8 clusters (MAFB C3,
      JAZF1 C8). `gene_patterns_k10_before_tS.csv` keeps the previous assignment.
  - 7.3m **full tradeSeq, D15** (user, 2026-11-21; replaces 7.3e): 22,824 genes (detected in >= 10 of the 5,243 cells on
    the beta / alpha / EC lineages), the three lineages only, 115 chunks of 200 genes, one acpu core (3.7 GB) each,
    serial fit: 2 min and 1.3 GB per chunk (the parallel 7.3e ran > 2 days). Prep writes per-chunk counts so a task
    never loads the whole matrix. Tests per chunk from the saved models: association; conditionTest per lineage,
    l2fc = 0, over the **shared range** (knots up to the KO cells' 95th percentile of pseudotime: beta knots 1-3, alpha
    1-4, EC 1-4) and over the whole curve; earlyDE. Wald eigenvalue threshold 1e-4 (the default 0.01 left 1,000-
    2,000 genes per lineage, e.g. FEV, without a result). Effect = max |log2((KO + 1) / (ctr + 1))| of the smoothed
    curves over the **well-populated KO range**: from the later of the two genotypes' 5th percentile of pseudotime on
    the lineage to the KO 95th percentile (beta 74-151, alpha 73-158, EC 83-184). (user, 2026-11-21: the first
    version measured from the curve start, where few cells sit; 255 of 362 beta calls peaked at pseudotime < 50, e.g.
    INS and CHGA "higher in KO" at 23.) Call = shared FDR < 0.05 and effect >= 0.5; direction and pseudotime of the
    maximum kept; TF flag (Lambert 2018).
    - Result: shared-range FDR < 0.05 beta 863 / alpha 1,018 / EC 1,305 genes; with the effect cut **402 genes**, 342
      of them also pseudobulk curve calls (of 2,125). ISL1 / PDX1 / FEV all three lineages, ARX alpha + EC, ETV1 and
      MNX1 alpha; LMO2 and NKX6-1 not called (their difference lies past the stop: pseudobulk "missing" / "kept on in
      KO").
    - Beta lineage: 184 genes (96 lower in KO, 88 higher; 149 also pseudobulk beta calls). TFs lower in KO: ISL1,
      PBX3, ETS2, PAX6, NPAS3, TSHZ2, PDX1, FOXO1; higher in KO: ESRRG, FEV, ZNF608. INS lower in KO (-1.28).
      `trajectory/tradeseq_full/tradeseq_{association,condition,earlyDE,condition_wide,smooth}.csv`. Outputs `trajectory/pseudobulk/`
      (`pb_curve.csv`, `pb_after_stop.csv`, `pb_ko_states.csv`, `pb_genes.csv`, `gene_patterns/`).
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
    - TFs per cluster (user, 2026-11-21): human TF list = Lambert et al. 2018, Cell 172:650 ("The Human Transcription
      Factors", humantfs.ccbr.utoronto.ca, TF_names_v_1.01.txt, 1,639 genes; downloaded 2026-11-21 with the user's OK,
      CRLF removed) in `study_config/annotation/TF_names_v_1.01.txt`. 169 of the 2,656 genes are TFs (C1 19: ISL1, PDX1,
      FOXO1, TBX3, PBX3, ...; C5 21: ARX, ETV1, PAX6, HHEX, MAF, NR2F1, ESRRG, MEF2C, ...; C8 16: FEV, LMX1A, MNX1, NR0B1,
      ASCL1, CUX2, ...; C10 16: NEUROD4, TWIST1, HEY2, NFIA, TBX5, JUN, ...). LMO2 is not in the list (LIM-only, no
      DNA-binding domain). Lists: scRNAseq/reports/trajectory/gene_patterns_k10/ (local, not in git).
  - Genes along every lineage (helpers `plot_gene_lineages_D15.py`, 6 D15 lineages, and `plot_gene_lineages_D15_D22.py`,
    8 lineages of the full D15 + D22 integration or 4 of the subset; from expression, any gene; user, 2026-11-21):
    LMO2, ETV1, PAX4, PAX6, PDX1, FOXO1, MEIS3, PBX3, NEUROD1, NEUROD4, LHX1, FEV, MNX1, ARX, NKX6-1, MAFB; page
    `results/trajectory/gene_lineages/index.html` (local). The user prefers the **D15-only data** for telling apart
    genes that drive alpha or beta downstream of ISL1. User's reading (hypotheses): ARX rises in KO along the alpha
    lineage and stops only because the KO trajectory stops (ISL1-independent); ETV1 is not switched on in KO (alpha,
    ISL1-dependent); LMO2 beta; PAX6 beta (not alpha); PDX1 beta; FOXO1 both. Checked on the curves: ETV1 also rises on
    the delta lineage and its control rise starts near the KO alpha stop (1-2 overlapping bins); LMO2's beta rise starts
    after the KO beta stop, so not switched on and never reached cannot be told apart yet (LMO2 also rises on alpha and
    delta); PAX6 is ISL1-dependent on beta but follows control on alpha; PDX1 and FOXO1 are lower in KO from the first
    bins on every lineage; NKX6-1 switches off in control pre_alpha but stays on in KO. A systematic D15 classification
    (switched on in control but not KO within the range KO reaches / follows control then truncated / lower from the
    start / undetermined) is proposed, not run.
  - 7.3i **alpha vs beta bias of the progenitors** (user, 2026-11-21). Experimental data (user): ISL1 overexpression,
    even in KO cells, rescues alpha and beta cells and reduces EC, so D15 pre_EC cells are probably still competent
    progenitors; EC is still kept as an end state, because without the manipulation cells do fall into EC (user: an
    alpha / beta-only end set is not correct). **No trajectory rerun**: the existing 7.3d / 7.4d fate probabilities (runs on
    all cells) are read as they are. Per cell bias = P(alpha) / (P(alpha) + P(beta)), S = P(alpha) + P(beta); pseudotime
    kernel = main, velocity = check.
    - D15: alpha = pre_alpha, beta = pre_beta2; cells ENP1, ENP2, ENP_FEV_high and pre_alpha_KO (user: include
      pre_alpha_KO); pre_EC is a D15 end state, so its cells have no bias. D15 + D22: alpha = SC_alpha, beta = SC_beta,
      EC = SC_EC1 + SC_EC2; cells D15 ENP1, ENP2, ENP_FEV_high, pre_alpha_KO and pre_EC (user).
    - Check 1 (pseudotime kernel): every cell has a ratio (S = 0 only for 2 KO D15 + D22 ENP_FEV_high cells). S < 0.05 in
      D15: ENP_FEV_high control 9 / KO 117; in D15 + D22 the D22 EC ends take most of the probability (P(EC) median
      0.78-0.98): ENP_FEV_high control 65 / KO 408, pre_EC control 106 / KO 581 of 712.
    - Median bias control -> KO, D15: ENP1 0.51 -> 0.51, ENP2 0.38 -> 0.41, ENP_FEV_high 0.21 -> 0.32, pre_alpha_KO
      0.77 -> 0.84. D15 + D22: ENP1 0.72 -> 0.73, ENP2 0.51 -> 0.58, ENP_FEV_high 0.23 -> 0.41, pre_alpha_KO 0.96 -> 0.94,
      pre_EC 0.21 -> 0.46 (control pre_EC is beta-leaning; KO pre_EC has mostly S < 0.05). Values are comparable
      within a dataset, not between the two (different end states).
    - No clear alpha / beta separation: one broad peak per subcluster. The bias follows position along the trunk
      (Spearman with pseudotime, control: D15 -0.70, D15 + D22 -0.60; within ENP_FEV_high -0.33 / -0.30); within
      ENP_FEV_high it is unrelated to the EC fate. KO lacks the strongly beta-leaning ENP_FEV_high cells: P(beta) about
      halved, P(alpha) unchanged, so the KO shift toward alpha is lost beta potential. A small alpha-leaning group sits at
      the earliest ENP1 tip in both genotypes. Velocity-kernel bias follows pseudotime almost entirely (velocity runs
      backwards).
    - Not decided: how to group the progenitors (proposed: rank the bias within pseudotime bins in control, check the
      split in expression, place KO cells with an expression classifier).
    - Outputs `day_analysis/D15/trajectory/alpha_bias/{D15,D15_D22}/`; report `results/trajectory/alpha_bias/index.html`
      (`analysis/scripts/build_alpha_bias_report.py`, Mac).
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
    - Done 2026-11-21 (scVelo 36100688, 11 min; CellRank 36100689, 3 min; all 6 runs written; the PETSc "Broken Pipe"
      lines are worker shutdown noise, no traceback; the KO-only pseudotime run found no data-driven macrostates).
      **Velocity runs backwards again**: latent time highest in D15 ENP1 (median 0.84) and the velocity kernel's only
      data-driven terminal state is ENP1 (all, ctr and KO runs); Spearman of latent time with Slingshot pseudotime is
      mixed (control beta 0.35, KO beta -0.70; EC2 control -0.53). The pseudotime kernel finds SC_beta, SC_alpha,
      SC_EC1, SC_EC2 as terminal, so fates and terminal states are read from the pseudotime kernel (as D15).
      Pseudotime kernel, mean fate by label x genotype: beta fate at D15 ENP_FEV_high control 0.10 vs KO 0.03, at
      pre_beta1 0.12 vs 0.03; alpha fate early unchanged (ENP_FEV_high 0.035 vs 0.033); KO early cells head to EC1 / EC2.
  - **7.4g / 7.4h / 7.4i: the D15 trajectory analysis repeated on D15 + D22** (user, 2026-11-21: "Go ahead"). Changes from
    D15: a cell is on a lineage only with curve weight >= 0.5 AND a label on the lineage's path (weight alone put many D22
    EC cells on the beta / alpha curves: beta KO 6,218 vs 2,243 on the path); programs from the pseudotime-kernel drivers
    only; every KO vs control test stays within one stage (D15 or D22 = separate batches); line-wide = called in >= 2/3 of
    the 35 window tests (24 genes, the same ones as D15: TSPYL5, CTSF, GRIN3A, SMC2, SLC7A2, ...).
    - 7.4g: beta stop window (D15 pre_EC / pre_beta1 / ENP_FEV_high) lost only ACVR1C, gained 105 (TBX5, OSR1, NFIA, ZEB2,
      ...), on after the stop and missing in KO 38 (LMO2, SLC17A6, ITGA1, ALDH1A1, APOA2, ...); alpha stop lost 93 (ARX,
      TLE4, RGS4, THRB, NPY, ...), gained 445, on after the stop and missing in KO 539 (GCG, CRH, UCN3, PAX3, PAX5, GAST,
      KIT, LMO2, ...).
    - 7.4h: KO never reaches the D22 beta / alpha branches (0-2 KO cells), so beta / alpha are tested on D15 only; EC2 and
      EC1 are tested at D22 with many KO cells (e.g. EC2 KO 1,000 vs control 180 per window). EC programs control vs KO:
      EC2 236 control only / 9 both / 0 KO only; EC1 45 / 7 / 18.
    - 7.4i: 4,542 genes; silhouette low (0.20 at k = 4); k = 4, 8, 9, 10 written (seed agreement 0.99 / 0.99 / 0.95 / 0.90).
      k = 10: C1 late beta rise at D22, KO flat (TTR, CALB2, SLC30A8, ACVR1C, PDX1); C2 late alpha rise, KO lags (GCG,
      ISL1, ARX, ETV1, PAX6, HHEX, MAF); C3 mid rise, KO higher and still rising at the stop (FEV, LMX1A, MNX1, NKX6-1);
      C4 mid-path peak turned off late on beta / alpha (NKX2-2, INSM1, FOXP2, SSTR2); C5 up in KO at the alpha stop
      (GAST, AMBP, GAL, XBP1, CUX2, NR0B1); C6 early, lower in KO from the start; C7 early progenitor genes falling
      (SOX4, DLK1, NEUROG3, PAX4); C8 high early in KO only (TBX5, SOX11, GATA3); C9 rise toward every end in both
      genotypes (CHGA, CRYBA2, SCG3, ASCL1); C10 rise toward EC2 at D22, KO = control (TPH1, PCSK1, HIF1A, NEUROD4,
      TWIST1). Overlap with the D15 k = 10 clusters (1,972 shared genes) is consistent: D15 C1 -> C1, C5 -> C2, C3 -> C7,
      C9 -> C4, C2 -> C6, C6 -> C8, C8 -> C3 / C5; D15 C10 ("switched to EC in KO") -> C10: at D22 control EC2 cells also
      express these genes, so in KO the EC program comes on early (at D15), not a KO-only program.
      Lists and TFs per cluster: scRNAseq/reports/trajectory/D15_D22/ (local).
      **Working classification for D15 + D22: k = 10** (user, 2026-11-21), as for D15: `gene_patterns_k10.csv`,
      `gene_pattern_clusters_k10.csv`, `gene_patterns_{curves,heatmap}_k10.png` in the subset trajectory folder (rerun:
      `sbatch --export=ALL,INT_SUBSET=alpha_beta_EC analysis/scripts/7.4i_gene_patterns.sb 8,9,10`).
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
  - Cluster 2 split and relabelled (user, 2026-11-22: "subcluster and assign pancreatic progenitor, fibroblast,
    acinar and hepatic like D40"): 7.1g split at sub-resolution **0.2** (as D40 r0.2 cluster 6), 4 groups:
    0 **hepatic** (139 cells, 61 control / 78 KO; AFP and FGB 100%, HNF4A 88%, CDX2 84%; RBP4, APOA2, APOA1),
    1 **hepatic** (427 cells, 302 / 125; AFP 93%, FGB 100%, HNF4A 73%, ONECUT1 / SOX9 ~50%; top markers DST,
    CACNA2D1, PCDH15, CDH2; at sub-resolution 0.5 a hepatic part, a CUBN / TRPM6 part and a small CHGA / SCG3
    endocrine-like part, mostly KO; called hepatic by the user), 2 **pancreatic progenitor** (182 cells, 160 / 22;
    ONECUT1 78%, KRT19 89%, SOX9 50%, PDX1 35%), 3 **fibroblast** (57 cells, 41 / 16; COL3A1 98%, DCN 83%).
    No acinar group: PTF1A, CPA1, CPA2 and CTRB2 in <= 2% of every group. The "non-endocrine" label is gone (7.1h
    now leaves out labels without cells). D22 hepatic 566, pancreatic progenitor 182, fibroblast 57 cells.
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
    - **Merged (user, 2026-11-22):** D40 GHRL+ and SC_delta (both r0.4 cluster 9) are one label, **SC_delta_GHRL**
      (166 cells, 96 control / 70 KO), colour = D22 SC_delta (label_colours.tsv). D22 keeps SC_delta and GHRL+.
      Redrawn: D40 UMAPs, composition test (7.1l) and plots. 7.1l D40 with the merge: SC_delta_GHRL 1.61% -> 1.02%,
      propeller FDR 0.19 (ns); the other D40 FDRs shift slightly (moderated variance re-estimated), all still < 0.01.
      7.1m rerun with the merge (2026-11-22): EC results identical; ec_corr_D40 / ec_dotplot_D40 show SC_delta_GHRL. Its
      KO and control parts differ: control r 0.62 between replicates, SST / ARX / ISL1; KO r 0.63, EC-like (TPH1,
      SLC18A1, TAC1; r 0.31-0.46 with EC pseudobulks, 0.25-0.30 with control SC_delta_GHRL), as noted at labelling.
- D40 marker dot plot: SOX2 (neural) moved to the last column (user, 2026-11-22).

- Label proportions KO vs control (step 7.1l, user, 2026-11-22; `day_analysis/composition_genotype/`): per day
  (D22, D40) and label, the share of each sample's cells; 2 control vs 2 KO samples as replicates (wells of one
  differentiation and one OCM reaction per day, not separate differentiations). Main test **propeller** (arcsine
  square-root proportions, limma ~ genotype, eBayes robust, BH FDR per day; written out with limma, as
  speckle::propeller); check: quasi-binomial GLM per label (NA when a genotype has 0 cells). Effect log2((KO + 0.1%)
  / (control + 0.1%)).
  - D22, propeller FDR < 0.05 for every label: lower in KO SC_alpha (27.5% -> 0.02%), SC_beta (35.1% -> 0.01%),
    SC_delta (2.1% -> 0), ENP (1.9% -> 0.08%), pancreatic progenitor (2.4% -> 0.25%), hepatic (5.3% -> 2.2%),
    fibroblast (0.6% -> 0.17%), GHRL+ (0.68% -> 0.25%, FDR 0.012); higher in KO SC_EC1 (10.0% -> 37.8%), SC_EC2
    (14.4% -> 59.3%).
  - D40: all lower in KO except SC_EC1 (7.6% -> 35.6%) and SC_EC2 (5.3% -> 62.0%) higher; GHRL+ not different
    (1.17% vs 0.91%, FDR 0.44); SC_delta FDR 0.02 (34 cells).
  - propeller FDRs are very small (to 1e-22) because limma borrows variance across labels; the GLM is more
    conservative (D22 FDR 0.002-0.08, e.g. ENP / GHRL+ / hepatic 0.08). The non-endocrine decreases in KO are partly
    compositional (EC labels take ~95% of KO cells).
  - (2026-11-22, later) D40 rerun with SC_delta_GHRL (see the D40 labels): GHRL+ / SC_delta rows above replaced.
  - Figures (local, `scRNAseq/plot/`, `analysis/scripts/fig_composition.R`): composition_{day}.pdf (bar = genotype
    mean, dots = replicates, FDR), composition_log2ratio.pdf.

- KO vs control SC_EC1 / SC_EC2 similarity (step 7.1m, user, 2026-11-22; `day_analysis/ec_similarity/`), per day.
  Label match first (7.1k, unchanged by the D22 relabelling): each day's SC_EC1 / SC_EC2 is the other day's best
  partner (EC2 19 shared markers, r 0.68; EC1 13, r 0.56; D22 SC_EC2 also r 0.53 with D40 SC_EC1 on all cells, 0.18
  on KO cells only).
  - Pseudobulk correlation (sample x label, >= 20 cells, TMM logCPM, day's 3,000 SCT B variable genes, gene-centred):
    KO EC1 vs control EC1 r 0.80-0.82 (D22) / 0.77-0.79 (D40), replicates of one genotype 0.84-0.98, KO EC1 vs control
    EC2 0.58-0.63; KO EC2 vs control EC2 0.90-0.91 / 0.84-0.86. Each KO EC type is closest to the same control EC type.
    Cross-day (KO and control EC): same label across days r 0.55-0.81, other EC label 0.29-0.67.
  - Pseudobulk DE (edgeR QL, 2 vs 2, FDR < 0.05, |log2FC| >= 0.5): KO vs control D22 EC1 1,152 / EC2 858, D40 EC1 982 /
    EC2 1,065; references on the same design: EC2 vs EC1 control 4,227 (D22) / 1,772 (D40), KO 5,343 / 3,766; EC1 vs
    beta (control) 4,605 / 3,858. EC identity genes not or barely changed (TPH1, FEV, LMX1A, SLC18A1, DDC, CHGA within
    |log2FC| 0.5 except D40 EC1 TPH1 -0.7, FEV +0.5); ISL1 RNA higher in KO EC (+0.8 to +2.3; KO allele still
    transcribed); D40 KO EC up: immediate-early genes (FOS, JUNB, EGR1, IER2, DUSP1).
  - kNN genotype mixing within EC1 / EC2 (30 SCT B PCs, K = 30, genotype-balanced subsamples, 20 seeds; 1 = mixed):
    D22 0.68-0.69, D40 0.52-0.53; replicates 0.96-1.01; EC1 vs EC2 0.07-0.16.
  - Reading: KO EC cells keep their EC1 / EC2 identity (same markers, closest to the same control type, far fewer DE
    genes than between cell types) but are not identical to control EC (about 1,000 DE genes; partial mixing).
  - Figures (local, `scRNAseq/plot/`, `analysis/scripts/fig_ec_similarity.R`): ec_corr_{day}.pdf, ec_dotplot_{day}.pdf,
    ec_de_counts.pdf, ec_mixing.pdf. D40 KO GHRL+ cells look EC-like on the dot plot (TPH1, SLC18A1, TAC1), as noted
    at labelling.

- Joint DESeq2 of the EC pseudobulks (step 7.1n, user, 2026-11-22; `day_analysis/ec_deseq2/`): SC_EC1 / SC_EC2 of
  the 8 samples (D22, D40 x control, KO x rep1, rep2) = 16 pseudobulks (SoupX counts, from 7.1m; 153-2,965 cells
  each); 21,022 genes with >= 10 counts in >= 2 pseudobulks. Day = differentiation and 10x run, so day effects
  include batch; the days are read as two independent experiments. Call = padj < 0.05 and |log2FC (MLE)| >= 0.5;
  ashr-shrunk log2FC for ranking / plots.
  - Model A, EC2 vs EC1, paired within sample (`~ group + group:rep + group:celltype`, group = day x genotype):
    D22 ctr 4,257 / D22 KO 3,971 / D40 ctr 1,897 / D40 KO 2,703 genes; average 3,817. Core genes (average call, same
    sign in all 4 groups): 1,535 EC2-high, 1,608 EC1-high. The EC1 / EC2 difference barely depends on genotype (209
    genes, KO - ctr) but strongly on day (3,119, D40 - D22; D22 differences larger). GO: EC1-high = development /
    patterning / epithelial differentiation / ECM (NEUROG3, SOX1, SOX2, ASCL1, PAX6, MEIS1, GLI3, VIM); EC2-high =
    synapse organisation, transmembrane transport, hormone metabolism (TAC1, TRPA1, CALCR, NEUROD4, NR3C2, EBF1).
  - Model B, KO vs control within each cell type (`~ 0 + celltype_genotype_day`, unpaired): EC1 D22 1,127 / D40
    1,068 / average 971; EC2 D22 823 / D40 1,149 / average 979. KO effect EC2 - EC1: 268 genes; D40 - D22: 1,045.
    Classes (average-of-days calls): shared EC1 + EC2 403 (140 up / 263 down; 120 called on both days in both cell
    types), EC1 only 560, EC2 only 568, opposite direction 8. KO effects in EC1 and EC2 correlate (r 0.70, all genes).
    Only 125 genes are called in all 4 day x cell-type tests (UpSet), so most single-test calls are day-specific.
    - Shared down in KO, TFs: ARX, ETV1, TBX3, PBX3, PAX4, SOX2, MECOM, NR2F1, HES6, TWIST1, SMAD9, FOXP2, EPAS1
      (several are ISL1-bound islet / alpha-program TFs). Shared up: ISL1 (KO allele transcribed), STAT1, IRF9, IFI6
      (GO: type I interferon / antiviral), DLX5, DLX6, immediate-early genes (FOS, JUNB, EGR1, EGR4).
    - GO: shared down = axonogenesis, synaptic signalling, muscle contraction; EC1-only down = growth factor / TGF-beta
      signalling, epithelial proliferation, digestive tract development; EC2-only down = regulation of secretion,
      signal release, hormone levels; EC2-only up = amine catabolism, mesenchymal differentiation.
    - GCG EC2 "up in KO" (D40 log2FC 22) rests on near-zero control counts; not interpreted.
  - GO: `7.1n_ec_go.R`, offline (org.Hs.eg.db GO BP incl. descendants, fgsea::fora, universe = tested genes, terms of
    15-500 genes, BH). Environments: DESeq2 / ashr from the DiffBind env, GO from isl1_downstream.
  - Figures (local, `scRNAseq/plot/`, `analysis/scripts/fig_ec_deseq2.R`): ec_deseq2_celltype_volcano.pdf,
    ec_deseq2_KO_scatter.pdf, ec_deseq2_KO_upset.pdf, ec_deseq2_heatmap_celltype.pdf, ec_deseq2_heatmap_KO.pdf,
    ec_deseq2_GO.pdf.

- Hormones in EC cells by genotype (user, 2026-11-22: "is INS higher in KO EC at D22?", then "how about GCG?";
  `7.1m_hormone_check.py`, one-off, `day_analysis/ec_similarity/hormone_check/cells.csv.gz`; plots
  `fig_hormone_check.R INS GCG` -> plot/{INS,GCG}_ec_distribution.pdf). INS:
  - Raw INS is in ~100% of control EC cells (D22 median 21-26 UMIs) and SoupX removed 28-33% of the EC INS in control
    vs 3-12% in KO (KO samples have no beta cells, so little INS background). After SoupX, cells with INS > 0: D22 EC1
    ctr 45% / KO 90%, EC2 31% / 86%; D40 EC1 23% / 78%, EC2 18% / 77%. KO EC cells carry low-level INS of their own
    (D22 EC1 median 25 UMIs); the control EC level is mostly background.
  - Beta-level INS (>= 10th percentile of control beta cells) in EC cells of both genotypes: D22 EC1 ctr 14% / KO 20%,
    EC2 10% / 9%; D40 EC1 11% / 5%, EC2 11% / 1%. These cells are EC cells, not beta doublets: scDblFinder median
    0.008-0.11, beta genes low (IAPP, G6PC2, MAFA, ERO1B, ADCYAP1, HADH mean 0.16-0.28 vs 0.8-1.1 in beta), TPH1 / FEV
    as other EC cells. They hold 65-95% of the summed EC INS, so the pseudobulk mostly reflects them.
  - Reading: D22 KO EC express INS in more cells at a low level and a similar share at beta level, so the per-sample
    total is not different (7.1n: EC1 +0.34, EC2 -0.40, ns); D40 KO EC have fewer beta-level INS cells, hence lower
    pseudobulk INS in KO EC2 (-2.56, padj 1e-7).
  - GCG (reference control SC_alpha, median 777 / 523 UMIs at D22 / D40): control EC GCG is background (raw > 0 in
    72-95% of cells, SoupX removed 78-100%; after SoupX 0-3% of cells). KO EC cells carry trace GCG (SoupX removed
    15-31%; > 0 in 62-65% at D22, 40-43% at D40; median 0-1, 99th percentile 4-8 UMIs). Alpha-level GCG (>= 10th
    percentile of alpha, 24 / 15 UMIs): 0 KO EC cells, 5 control EC cells (the 2 at D40 doublet-like). EC cells do not
    make glucagon; the 7.1n "GCG up in KO EC2" (D40 log2FC 22) is trace KO reads against a control sum of 0.

- Marker figures (user, 2026-11-22; `fig_export_genes.py` on Alpine -> `day_analysis/figure_data/{day}_genes.csv.gz`;
  `fig_feature_dot.R` local -> plot/feature_{day}/{gene}.pdf, plot/dotplot_{day}.pdf). Feature plots: control | KO,
  the other genotype light grey, grey (0) to red, one scale per gene. Dot plots: control and KO side by side per gene,
  no grid / cell numbers, x = group with < 20 cells. Markers chosen on % of control cells expressing:
  beta G6PC2 (17-31% of beta, GCK / SLC2A1 broad), pancreatic progenitor SOX9 (53-84%; PDX1 high in beta / EC; PTF1A
  dropped), fibroblast COL3A1 (97-98%, near absent elsewhere; COL1A1 39% of D22 fibroblasts), D40 neural SOX2 (79% of
  neural, <= 20% elsewhere; DCX / ELAVL4 / TUBB3 higher in endocrine cells).
  - D22: GCG ARX INS G6PC2 NKX6-1 SST GHRL NEUROG3 SLC18A1 TPH1 FEV AFP SOX9 COL3A1.
  - D40: GCG ARX INS G6PC2 NKX6-1 SST GHRL SLC18A1 TPH1 FEV AFP SOX9 PRSS1 COL3A1 SOX2 (no NEUROG3; PRSS1 also high
    in D40 hepatic cells).

- D15 whole-day labels (user, 2026-11-22; before the endocrine subset; `study_config/annotation/D15_annotation.tsv`,
  7.1h -> `day_analysis/D15/labels/`): overview resolution 0.2, base rules only: 0 **ENP** (5,178), 1 **pre_alpha**
  (620, control only), 2 **PP-hepatic** (3,844; AFP / APOB / CDX2 high, PDX1 / NKX6-1 low; named PP-hepatic by the user,
  who does not expect that many hepatic cells in this differentiation), 3 **PP2** (2,374), 4 **PP1** (3,806), 5
  **PP-ENP** (486), 6 **fibroblast** (433). The 290 cells left out of the 7.1d endocrine subset (83 in cluster 0, 4 in
  1, 203 in 5) keep their cluster's label (user). Relation to the 7.1d subtypes: ENP holds ENP2, ENP_FEV_high,
  ENP_PHOX2A+, pre_EC, pre_beta1 / 2, pre_alpha_KO, KO_only; pre_alpha = 526 pre_alpha + 87 pre_delta; PP-ENP = 265 ENP1
  + 203 excluded cells. Colours added to label_colours.tsv (pre_alpha, PP-ENP, PP1, PP2, PP-hepatic; ENP and
  fibroblast reuse the D22 / D40 colours). Note: D22 / D40 call their AFP-high cells "hepatic".
  - Figures (local, plot/D15/): `PLOT_DIR=plot/D15 Rscript analysis/scripts/fig_umap_labels.R D15` -> umap_D15.pdf,
    umap_D15_genotype.pdf, umap_core/umap_D15{,_ctr,_KO}.png; exploratory cluster plots `fig_d15_overview.R`.
- 2026-11-22 (user: "make the D15 genotype composition plots like D22 and D40"): 7.1l rerun as `D15 D22 D40` (BH FDR
  is per day, so D22 / D40 results unchanged). D15 propeller: KO higher ENP (24.6% -> 38.9%, FDR 0.002); KO lower
  pre_alpha (6.7% -> 0%, FDR 1e-4), PP1 (FDR 0.008), PP2 (FDR 0.02); PP-ENP, PP-hepatic, fibroblast ns.
  `fig_composition.R` takes `DAYS` (default D22,D40): `DAYS=D15 PLOT_DIR=plot/D15` -> composition_D15.pdf,
  composition_log2ratio_D15.pdf, composition_stacked.pdf. Panel title wrapped to two lines, stacked caption on three
  lines (no clipping; also applies to the D22D40 versions).
- 2026-11-22 (user): UMAP of the D15 endocrine compartment (7.1d subset, 5,994 cells, its own 2D UMAP, r0.6 names).
  Colours (label_colours.tsv; made more distinct on request, same day): lighter shades of the D22 / D40 colours for
  the precursors: pre_alpha #5fa2d3 and pre_alpha_KO #b4d6ec (SC_alpha), pre_beta2 #fd8d3c and pre_beta1 #fdd49e
  (SC_beta), pre_EC #fb9a99 (SC_EC1), pre_delta #c49c94 (SC_delta); other hues: ENP1 #bcbd22 (olive), ENP2 #2ca02c
  (green), ENP_FEV_high #8e7cc3 (lavender), ENP_PHOX2A+ #4d4d4d (dark grey), KO_only #e7298a (magenta). pre_alpha is
  shared with the whole-day D15 label, so umap_D15*.pdf / png were redrawn too. `PLOT_DIR=plot/D15 Rscript analysis/scripts/fig_umap_labels.R D15_endocrine` ->
  umap_D15_endocrine{,_genotype}.pdf, umap_core/umap_D15_endocrine{,_ctr,_KO}.png.
- 2026-11-22 (user: "compare cell proportion by genotype here", D15 endocrine compartment): 7.1l takes
  `D15_endocrine` (7.1d subset, r0.6 names; sample from the barcode; proportions are of the subset's cells per
  sample); run as `D15 D22 D40 D15_endocrine` (FDR per day, other days unchanged). Propeller, KO vs control: lower
  pre_alpha (18.3% -> 0%), pre_beta2 (12.1% -> 0%), pre_delta (4.0% -> 0%), pre_beta1 (13.0% -> 3.4%); higher pre_EC
  (7.8% -> 23.8%), ENP_FEV_high (19.4% -> 29.1%), ENP_PHOX2A+ (7.0% -> 11.5%), ENP2 (12.0% -> 14.5%, FDR 0.04),
  pre_alpha_KO and KO_only (named for their genotype, so their test is circular); ENP1 ns. All others FDR < 1e-4.
  `DAYS=D15_endocrine PLOT_DIR=plot/D15 Rscript analysis/scripts/fig_composition.R` (lineage order as the UMAP).
  Stacked-bar file per day set: composition_stacked.pdf for D22 + D40, otherwise composition_stacked_{days}.pdf
  (the D15 one is now composition_stacked_D15.pdf); display names without the underscore; one-day star table headed
  "FDR".
- 2026-11-22 (user): marker dot plot of the D15 endocrine compartment (control | KO columns, as D22 / D40), genes
  picked from per-subtype % expressing in control cells (one-off query, not kept) and the 7.1d r0.6 markers:
  general ENP NEUROG3, NEUROD1, SOX4; ENP1 DLK1 (79% vs <= 53%); FEV; ENP_PHOX2A+ PHOX2A (89% vs <= 37%); EC LMX1A,
  ADRA2A (pre_EC 90% vs <= 38%); beta NKX6-1, MNX1, PDX1, INS; pre_beta2 = beta-specific set: ISL1 (95%; pre_EC 43%),
  HADH (89% vs <= 45%, delta 83%), NPTX2 (60% vs <= 33%), with NKX6-1 (86%; alpha 9%, delta 17%); pre_beta1 FXYD2
  (70%; pre_EC 49%, others <= 50%); alpha ARX, GCG; delta SST, HHEX. pre_beta1 has no beta gene that separates it
  well from pre_EC (shares NKX6-1, MNX1, LMX1A, FEV, INS; low ISL1 62%, HADH 44%); its top 7.1d markers are
  liver-like (SERPINA1, APOA1, FTCD). ENP2 marker calling (user, same day; one-off query, table on scratch
  day_analysis/D15/subcluster_endocrine/ENP2_markers.csv, local copy reports/labels/D15/): control cells, ENP2 vs
  every other subtype (mean >= 1.5x the next subtype, >= 10 points more % expressing, >= 30% of ENP2, Wilcoxon vs
  rest padj < 0.01) -> 16 genes, none a classic endocrine gene; the next subtype is mostly ENP_FEV_high. Added to the
  plot (user): TNFRSF11B (50% vs 21%; osteoprotegerin, beta-cell proliferation) and PIEZO2 (53% vs 32%; EC
  mechanosensor). D15 figure data re-exported with 11 more genes (fig_export_genes.py D15, 40 genes).
  `DOT_ONLY=1 PLOT_DIR=plot/D15 Rscript analysis/scripts/fig_feature_dot.R D15_endocrine NEUROG3 NEUROD1 SOX4
  DLK1 TNFRSF11B PIEZO2 FEV PHOX2A LMX1A ADRA2A NKX6-1 MNX1 PDX1 INS ISL1 HADH NPTX2 FXYD2 ARX GCG SST HHEX` -> dotplot_D15_endocrine.pdf.
- 2026-11-22 (user): D22 pseudobulk correlation heatmap (fig_ec_similarity.R, ec_corr_D22.pdf) without GHRL+ (its
  single control pseudobulk); display only, the 7.1m correlations are unchanged (gene centring over all pseudobulks of
  the day). Page 0.3 in wider so the legend title is not clipped (D40 too).
- Figures regenerated on Alpine (user, 2026-11-23: "regenerate them on Alpine, follow the color schemes and all
  settings"; `plot/` was only on the other Mac). The repo's `fig_*.R` scripts, unchanged, with the commands and
  settings above (gene lists, DAYS / PLOT_DIR / DOT_ONLY), in `/scratch/alpine/hzhu10@xsede.org/ISL1_scRNAseq/fig_regen/`:
  `reports/` there links to the matching Alpine outputs (day labels, SCT UMAPs, composition_genotype, ec_similarity,
  ec_deseq2, figure_data), `bin/gzcat` wraps zcat, UpSetR (missing from single.cell) is in its own `Rlib/` (the env is
  unchanged). Job script `/scratch/alpine/hzhu10@xsede.org/fig_regen.sb` (13 commands); 141 files, copied to
  `scRNAseq/plot/` on the Mac.
- **D15 endocrine trajectory figures** (user, 2026-11-23: "make the figures 1-8 from existing outputs"; Slingshot and
  CellRank are the main trajectory analyses; independence of the two is not a concern for the user). Nothing rerun.
  `analysis/scripts/fig_d15_trajectory.R` -> `plot/D15/trajectory/`, style as fig_umap_labels.R (6 cm square UMAP panels,
  label_colours.tsv in lineage order, control | KO with the other genotype grey); every pseudotime is Slingshot's (7.3a,
  PCs 1-30; the CellRank pseudotime kernel uses it too). 1 `slingshot_umap.pdf` lineage curves over the labelled cells;
  2 `lineage_tree.pdf` subclusters at their UMAP centres with the lineage edges; 3 `ko_stop_density.pdf` pseudotime
  density per lineage (curve weight >= 0.5), control vs KO, condiments progression statistic, KO 95th percentile;
  4 `cellrank_fates_umap.pdf` fate probability per end state (pseudotime kernel, run on all cells); 5
  `cellrank_fate_heatmap.pdf` mean fate per subcluster x genotype (>= 20 cells); 6 `cellrank_end_states.pdf`
  data-driven macrostates of the all-cell / control-only / KO-only runs (the control-only run also calls ENP1
  terminal); 7 `alpha_bias_umap.pdf`, `alpha_bias_hist.pdf` (7.3i, D15); 8 `gene_dynamics.pdf` binned curves along the
  beta, alpha and EC lineages as 7.3h for ISL1, PDX1, FOXO1, PBX3, MEIS3, PAX6, LMO2, NKX6-1, MNX1, ARX, ETV1, PAX4,
  MAFB, NEUROD1, NEUROD4, LHX1, FEV (gene list = argument). D15 figure data re-exported with 10 more genes (ETV1,
  LMO2, PAX4, PAX6, FOXO1, MEIS3, PBX3, NEUROD4, LHX1, MAFB; 50 genes, the 40 before unchanged). Inputs:
  `reports/trajectory/D15/` = `day_analysis/D15/trajectory/`. Rerun on Alpine: `sbatch /scratch/alpine/hzhu10@xsede.org/fig_traj.sb`.
- Which gene findings depend on CellRank (2026-11-23): the main ISL1-dependent gene lists (7.3m full tradeSeq, 7.3l
  replicate pseudobulk, the k = 10 patterns of 7.3h on them, 7.3g window calls, the 7.3f after-stop calls) use Slingshot
  and expression only; CellRank gives the lineage labels on genes (driver / program flags, 7.3d / 7.3j), the fates,
  end states and alpha / beta bias (7.3i), the 7.3f / 7.4g "lost in KO" calls (CellRank-gated) and part of the 7.3k
  gene selection.
- **D22 KO EC vs control EC vs control beta** (step 7.1o, user, 2026-11-23). Question (user): KO EC contain INS-high
  cells (experimentally, at D22 the KO INS+ cells are SLC18A1+, i.e. EC-derived, while control INS+ cells are
  SLC18A1- beta cells); how do these KO EC differ from control EC, other KO EC and control beta, and do some KO EC keep
  beta differentiation potential? All differential tests use SoupX counts (also 7.1m / 7.1n; raw counts only in the
  7.1m hormone check).
  - INS distribution (`fig_ins_distribution.R`, `plot/D22D40/ec_ins/INS_distribution.pdf`; INS per 10,000 SoupX
    counts): KO EC carry INS spread over the whole population, no separate INS-high peak (D22 EC1 a flat continuum up
    to beta level; >= control-beta 10th percentile: KO EC1 12%, EC2 4%; control EC1 10%, EC2 7%, mostly 0 after
    SoupX with a sparse positive group); D40 KO EC mostly trace INS (4% / 1%). So no INS cut-off (user: "KO EC in
    general has more INS and it's spread"): a separate D22 model of ctr EC1, ctr EC2, ctr beta, KO EC1, KO EC2, plus a
    single-cell layer (user: "go ahead with the two-layer design").
  - Layer 1 (`7.1o_ec_beta_model.R`): edgeR QL on the 7.1m pseudobulks (>= 20 cells; KO SC_beta 1 cell dropped),
    ~ 0 + group + replicate within genotype; 19,389 genes; call FDR < 0.05, |log2FC| >= 0.5. KO EC1 vs ctr EC1 313 up /
    682 down, EC2 261 / 509; KO EC1 vs ctr beta 2,339 / 2,808, EC2 3,369 / 3,611; reference ctr beta vs ctr EC1
    2,482 / 2,309, EC2 3,649 / 3,239. Shift of KO EC toward beta along the control beta-EC axis (2,834 / 3,906 genes,
    |log2FC| >= 1, INS out; 0 = like control EC, 1 = like beta): EC1 0.19 (Spearman 0.29), EC2 0.09 (0.09). Quadrants of
    the KO calls: EC-high down in KO 263 / 105, beta-high up in KO 59 / 60, beta-high down in KO 142 / 193 (partly beta
    ambient RNA in control samples), EC-high up 56 / 50: the small EC1 shift is mostly loss of EC genes, not gain of beta
    genes. Group-mean correlation with control beta (axis genes, centred): KO EC1 -0.34 vs ctr EC1 -0.13, KO EC2 -0.64
    vs ctr EC2 -0.53.
  - Layer 2 (`7.1o_ec_beta_cells.py`): scanpy score_genes, 100 beta-high and 100 EC-high genes (control beta vs both
    control EC types, INS out); beta index = beta - EC score. Beta-leaning KO EC = beta index above the 99th (95th)
    percentile of the control cells of the EC type: KO EC1 0 (5) of 3,416, EC2 0 (33) of 5,374, fewer than chance at
    the 95th; KO EC beta score maximum -0.054 / -0.056 vs control EC 0.129 / 0.059 and control beta median 0.103; KO EC
    with > 10% control-beta neighbours (SCT B graph) 1 / 0 vs control EC 39 / 18; KO EC cells at or above the control
    beta 1st percentile of the index: 2 / 0. The beta score does not rise with INS (KO EC INS quintiles, median INS
    1 -> 142 / 0 -> 63 per 10k; beta score median -0.173 -> -0.165; Spearman beta index vs INS 0.04 / 0.11). The subset
    DE (beta-leaning vs other KO EC) was skipped: no beta-leaning cells.
  - Reading: at D22 KO EC are EC cells as far from beta as control EC; INS in KO EC comes without the rest of the beta
    program (PDX1, NKX6-1, MAFA, IAPP, G6PC2, HADH ...), and the scRNA data show no KO EC subset with partial beta
    identity. This is a snapshot of state; potential (e.g. the ISL1 overexpression rescue) is not tested by it.
    Control EC carry beta ambient RNA, which raises the control EC range; KO EC also stay below the lowest control
    beta cells.
  - Outputs `day_analysis/ec_beta_D22/` (de_{contrast}.csv, de_wide.csv, beta_shift.csv, quadrants.csv, pb_corr.csv,
    score_genes.csv, cells.csv.gz, cells_all.csv.gz, leaning_counts.csv, cells_summary.csv). Figures (`fig_ec_beta.R`,
    `plot/D22D40/ec_beta/`): ec_beta_quadrant.pdf, ec_beta_scores.pdf, ec_beta_vs_INS.pdf, ec_beta_corr.pdf,
    ec_beta_identity_heatmap.pdf, ec_beta_umap_{beta,EC}_score.pdf.
  - Figure choices (user, 2026-11-23): the quadrant plot keeps every axis gene coloured by quadrant and labels only
    PDX1, G6PC2, IAPP, PCSK1, MAFA, ERO1B, TPH1, SLC18A1, FEV, LMX1A, DDC as large outlined dots (their KO call per EC
    type: ERO1B lower in KO in both, MAFA in EC1 only, the rest not changed; PCSK1 is not on the EC2 axis). Note for
    reading it: the coloured groups are dominated by genes that are not beta / EC identity genes (side hormones
    PYY / NTS / CRH up, likely ambient / contaminant genes such as DCN / GJA1 down, KO-line-wide CTSF / DLX5 / DLX6);
    the message is the near-flat fitted shift.
    Identity-gene heatmap: 18 beta and 8 EC identity genes, scaled group means of the pseudobulks next to log2FC KO EC
    vs control EC (* = called): only ISL1 (KO allele transcribed), SLC30A8 (EC1), ADRA2A (EC1) higher in KO EC; MAFA,
    HADH, ERO1B (and ADCYAP1) lower; PDX1, MAFB, IAPP, G6PC2, GCK, ABCC8 and the EC genes TPH1, SLC18A1, DDC, FEV, LMX1A
    not changed; INS pseudobulk not changed (+0.3 / -0.4). At D22 NKX6-1, NKX2-2 and NEUROD1 are higher in EC than in
    beta (both genotypes). Beta / EC scores on the whole-day UMAP: all D22 cells scored with the same gene sets
    (`cells_all.csv.gz`; the subset scores and thresholds above are unchanged), control | KO, beta score blue and EC
    score green (user: gene expression stays red).

## Open items

- Step 1.4 does not yet count the NF > 0.6 cells removed at the per-day merge.
- Final report caveats: D40 KO SoupX fallback; scDblFinder run-to-run variation; the 2,000-gene rule; the NF rules;
  the strong ctr / KO separation on D22 and lower cycling in D22 / D40 KO.
- October-only output `HZ_Diff34_D40_scRNA`: decide whether to keep.
