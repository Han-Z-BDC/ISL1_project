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
  = draft until confirmed). Next: confirm the resolution and labels.

## Open items

- Step 1.4 does not yet count the NF > 0.6 cells removed at the per-day merge.
- Final report caveats: D40 KO SoupX fallback; scDblFinder run-to-run variation; the 2,000-gene rule; the NF rules;
  the strong ctr / KO separation on D22 and lower cycling in D22 / D40 KO.
- October-only output `HZ_Diff34_D40_scRNA`: decide whether to keep.
