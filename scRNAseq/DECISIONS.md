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

## Current state (2026-11-19)

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
  labels, stage mixing, labels and genotype kept, nearest D22 label per D15 state, marker genes). Method not chosen
  yet. D15 pre_delta lands next to D22 SC_delta (Harmony 74%, RPCA 84%, scVI 58% of its cells).
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
  - 7.3b velocyto per D15 sample (spliced / unspliced, Cell Ranger BAMs); 7.3c scVelo dynamical model with arrows on
    the 2D and 3D UMAPs; 7.3d CellRank 2 with a velocity kernel and a Slingshot-pseudotime kernel (fates, terminal
    states, driver genes for pre_alpha / pre_beta2); 7.3e tradeSeq on all genes detected in >= 10 cells (control
    vs KO condition tests per lineage, early DE at branches, gene modules by curve shape). CellRank runs in its own
    env (`CELLRANK_ENV`; CellRank 2.0.7 pins Jinja2 < 3.1, single.cell has 3.1.5). Running on 2026-11-19.
  - Next: comparison report (Slingshot curves, velocity, CellRank fates, KO drop-off per lineage, tradeSeq modules);
    then cross the candidate genes with ISL1 CUT&RUN and bulk KO.
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
