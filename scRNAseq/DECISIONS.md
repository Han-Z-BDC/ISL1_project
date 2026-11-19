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
  (CHGA+/FEV+ endocrine; control-only ISL1+/SCG5+; DCC+/TOX+ endocrine progenitors): 6,284 cells.
- The compartment is re-processed on its own: SCT B refit on its cells only (`7.1c_day_sct.R` subset mode, same
  settings), 30 PCs, 15 neighbours, UMAP in 2D and 3D, Leiden 0.2-1.2. Markers vs the rest of the compartment.
  Output: `day_analysis/D15/subcluster_endocrine/` (report, `umap3d.html` interactive 3D UMAP).

## Current state (2026-11-18)

- 7.1b clustering on SCT B done for D15, D22 and D40.
- D15: overview resolution 0.2; endocrine subclustering (7.1d) done. Report with interactive 3D UMAP (7.1e) at
  subcluster resolution 0.6 with **draft** labels (`study_config/annotation/D15_endocrine_r0.6.tsv`; status column
  = draft until confirmed). The report has interactive 2D (all / control / KO) and 3D UMAPs.
- D15 is **paused** here (2026-11-18, user: complicated, come back later). Not decided yet: the subcluster
  resolution (0.6 is a working choice), the endocrine labels (all draft), and whether to subcluster the other D15
  compartments. Summary so far: hormone+ cells (beta-, alpha-, delta-like) are almost control-only; KO cells sit in
  FEV+, EC-like (LMX1A+) and two KO-only states.
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

## Open items

- Step 1.4 does not yet count the NF > 0.6 cells removed at the per-day merge.
- Final report caveats: D40 KO SoupX fallback; scDblFinder run-to-run variation; the 2,000-gene rule; the NF rules;
  the strong ctr / KO separation on D22 and lower cycling in D22 / D40 KO.
- October-only output `HZ_Diff34_D40_scRNA`: decide whether to keep.
