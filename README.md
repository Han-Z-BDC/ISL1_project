# single-cell-pipeline

The lab's standardized single-cell (RNA / multiome) analysis pipeline.
Config-driven, one module per QC/processing step, built by comparing the
lab's own validated CellRanger + notebook methodology against IGVF-UCSD's
`multiome-sample-qc` / `CellCommander` and the sc-best-practices.org
reference, then picking one evidence-based default per step.

This repo is self-contained -- clone it into your own folder and run it
there. Nothing in it points at any specific person's storage; every path
you need is set once, in your own `config.sh`, which you create from the
template and never commit.

**This is the standardized pipeline for future analyses.** The lab's
existing per-project notebooks and job templates
(`/projects/hzhu10@xsede.org/{jobs/templates/single_cell,notebooks/new_notebooks}`)
are kept exactly as they are, unchanged, as the reference record for
publication -- they are not superseded in place, and this repo does not
modify them.

## Status

| Step | Module | Status |
|---|---|---|
| 1.1 | RNA QC: count/gene pre-threshold + MT%/Ribo% cutoff | **Implemented** -- `scripts/1.1_rna_qc.py` |
| 1.2 | Ambient RNA removal (SoupX) | **Implemented** -- `scripts/1.2_remove_background.py` |
| 1.3 | Doublet detection (scDblFinder + Scrublet) | **Implemented** -- `scripts/1.3_detect_doublets.py` |
| 2.1 | ATAC QC, AMULET | **Implemented** -- `scripts/2.1_atac_amulet.py`, run separately (slow) |
| 2.2 | ATAC QC (feature space, general QC, doublets) | **Implemented** -- `scripts/2.2_atac_qc.py` |
| 3.1 | Reconcile RNA + ATAC barcodes per sample | **Implemented** -- `scripts/3.1_merge_modalities.py` |
| 3.2 | Clean up intermediate h5ad, once 3.1 succeeds | **Implemented** -- `scripts/3.2_cleanup_intermediate_h5ad.sh` |

Numbering: `1.x` = RNA QC track, `2.x` = ATAC QC track (independent,
can run in parallel with 1.x), `3.x` = reconciling the two tracks once
both are done. Steps are being designed and built one at a time, in
order. Within the RNA track, 1.1's two criteria (count/gene threshold +
MT%/Ribo%) could be combined freely because both are just boolean cell
masks (order-independent); 1.2→1.3 are NOT interchangeable, so each
depends on the previous one actually being decided first, not guessed at
ahead of time. Normalization (what would have been 1.4) isn't scaffolded
-- not needed yet; add it when it's actually being designed.

Step 1.2's method (SoupX with `autoEstCont`, no marker-gene file) was
chosen by comparing the lab's own notebook, CellCommander/IGVF-UCSD, and
sc-best-practices.org -- see `configs/remove_background.yaml` for the full
reasoning. Note that IGVF-UCSD's SoupX step *is* CellCommander (their
pipeline calls the CellCommander CLI directly), so that comparison was
really two implementations plus one reference, not three independent ones.

Step 1.3 uses the lab's own validated method as-is: scDblFinder + Scrublet,
where Scrublet's `expected_doublet_rate` is calibrated per sample off
scDblFinder's own call rate, and a cell is removed if EITHER tool flags it
(union, not consensus) -- see `configs/detect_doublets.yaml`.

**ATAC QC** was compared the same way (lab notebook vs. IGVF-UCSD/
CellCommander vs. sc-best-practices.org) and turned out architecturally
different across all three, not just parameter-different -- see
`configs/atac_qc.yaml` for the full three-way comparison and reasoning.
Headline decisions: feature space is SnapATAC2's genome-wide tile matrix
built from fragments (not CellRanger's own called peaks, which the lab's
notebook and sc-best-practices both use -- the one place this pipeline
matches only 1 of 3 references); doublets are a union of three detectors
(SnapATAC2's own scrublet + AMULET + scDblFinder with
`aggregateFeatures=TRUE`); and scDblFinder runs with `clusters=FALSE`
rather than sc-best-practices' own tutorial default of `TRUE`, since their
own guidance says cluster-based scoring isn't appropriate for
continuous-trajectory data like this lab's differentiation time-course.

**AMULET runs as its own separate step** (`scripts/2.1_atac_amulet.py`),
not inline inside `2.2_atac_qc.py` -- it's slow and its runtime doesn't
depend on QC thresholds (it scans the whole raw fragment file regardless
of any filter). The lab's own production job for this
(`jobs/templates/single_cell/batch_amulet.sb`, calling
`software/python_scripts/run_Amulet.py`) requests 3 days for a sequential
loop over every sample -- confirming this needed decoupling, not folding
into one combined job. `2.2_atac_qc.py` reads AMULET's cached CSV output
rather than calling it itself.

## 1. Set up the environment

```bash
acompile --ntasks=4 --time=01:00:00      # build environments on a compute node, not the login node
module load miniforge
conda env create -f environment.yml      # builds an env named "sc-pipeline"
conda activate sc-pipeline
```

`config.sh` (next step) points `TOOLS` at this environment's `bin/`
directory so every pipeline script calls tools by full path rather than
relying on `conda activate` inside a batch job.

## 2. Configure

```bash
cp config.example.sh config.sh
nano config.sh
```

Fill in `IN_DIR` (your CellRanger outputs root), `OUT_DIR`, and `TOOLS`
(the env `conda env create` just built -- find it with `conda env list`).
`config.sh` is gitignored on purpose.

## 3. Describe your samples

```bash
cp samples.example.tsv samples.tsv
```

One row per sample: `sample_id` and the path to its
`raw_feature_bc_matrix.h5` (relative to `IN_DIR`) -- this is the anchor
identifying the sample's CellRanger `outs/` directory, and every other
file this pipeline needs (`filtered_feature_bc_matrix.h5`,
`atac_fragments.tsv.gz`) is derived from it as a sibling in that same
directory, assuming CellRanger's standard layout.

## 4. Set QC parameters

Edit `configs/qc_RNA.yaml` -- it's versioned, not gitignored, because it's
the record of what thresholds a given run actually used, not
account-specific config. Defaults (5-MAD count/gene threshold, MT% cutoff
at 3%, no Ribo% cutoff) are the lab's agreed starting point; see the
comments in that file for the reasoning and how to override per sample.

## 5. Run it, in order

```bash
bash scripts/submit_1.1_rna_qc.sh            # array job -- don't sbatch 1.1_rna_qc.sb directly
bash scripts/submit_1.2_remove_background.sh # array job -- needs step 1.1 finished first
bash scripts/submit_1.3_detect_doublets.sh   # array job -- needs step 1.2 finished first

# ATAC QC -- independent of steps 1.1-1.4, can run in parallel with them.
bash scripts/submit_2.1_atac_amulet.sh       # array job -- start this first, it's slow (see 2.1_atac_amulet.py)
bash scripts/qc_only_atac.sh <sample_id> # interactive (acompile), not sbatch -- repeat until threshold_review.png looks right
bash scripts/submit_2.2_atac_qc.sh           # array job, full pipeline -- needs the AMULET job done; reuses qc_only_atac.sh's cached checkpoint

# Once BOTH RNA (through 1.3) and ATAC (through 2.2) are done for a sample:
bash scripts/submit_3.1_merge_modalities.sh  # array job -- reconciles RNA + ATAC barcodes, one final h5ad per modality

# Once 3.1 is done for a sample -- single-sample only, NOT cross-sample/replicate integration (a future step):
bash scripts/submit_4.1_rna_normalize.sh     # array job -- RNA normalize + UMAP, two parallel methods (log, SCTransform)
bash scripts/submit_4.2_atac_normalize.sh    # array job -- ATAC normalize + UMAP (SnapATAC2 spectral embedding)
```

Step 1.1, for every sample, writes an annotated h5ad (all cells + QC
metrics + pass/fail flag), a filtered h5ad (QC-passing cells only), a
before/after distribution plot, and a summary CSV -- check the plot and
summary before trusting the filter for a new dataset; the MAD default
self-calibrates per sample, but a very unusual sample (e.g. very low
sequencing depth) is worth eyeballing rather than assuming.

Step 1.2 (SoupX) consumes step 1.1's filtered h5ad plus the same raw h5,
and writes a SoupX-corrected h5ad, a pre/post UMAP diagnostic plot, and a
summary JSON (total counts removed, top genes affected) per sample --
matching the lab's own notebook, which doesn't report a contamination
percentage either.

Step 1.3 (doublet detection) consumes step 1.2's SoupX-corrected h5ad,
and writes an annotated h5ad (both tools' scores + calls, for every
cell), a doublets-removed h5ad (the union-filtered result actually meant
for downstream use), a diagnostic plot (scDblFinder score histogram + a
per-tool doublet-call breakdown), and a summary CSV.

ATAC QC's AMULET step consumes `atac_fragments.tsv.gz` directly and
writes `{sample}.amulet_result.csv` (barcode, p.value, q.value,
amulet_predicted_doublet) plus a score histogram -- run this first, for
every sample, since it's slow and independent of everything else. The
main ATAC QC step then also consumes `atac_fragments.tsv.gz` (assumes
CellRanger ARC's output layout -- a sibling of the RNA h5 in the same
`outs/` directory) plus that AMULET CSV, and
writes an annotated h5ad (all cells, all three detectors' scores + a
diagnostic UMAP/clustering), a doublet-union-filtered final h5ad, plots
(fragment size distribution, threshold review, doublet-call breakdown,
UMAP), and a summary CSV. It does not depend on RNA QC steps 1.1-1.4 and can run in
parallel with them.

Fragment-file parsing (the slow part) is cached to
`{sample_id}.import_tsse_checkpoint.h5ad` the first time either
`qc_only_atac.sh` or the full pipeline runs for that sample -- a later run
of either one reuses it instead of reparsing, so tuning thresholds with
`qc_only_atac.sh` costs nothing extra when you move on to the full run.
Delete that checkpoint file to force a clean re-import (e.g. if the
fragment file itself changed).

Once both tracks are done for a sample, `3.1_merge_modalities.py`
intersects RNA's and ATAC's final QC-passing barcodes (stripping ATAC
QC's `{sample_id}#` prefix first) and subsets both objects to the shared
cells -- matching the lab's own notebook's method
(`3.single.timepoint.joint.analysis.ipynb`) exactly. Treated as just
another filter step like everything before it: writes
`{sample}.rna_final.h5ad` and `{sample}.atac_final.h5ad` -- plain h5ad,
same convention as every step before it, no MuData bundling (not needed
yet; trivial to add later if it becomes useful --
`mu.MuData({"rna": rna, "atac": atac})`). This is where the two
independent QC tracks actually come together; everything past it
(replicate merging, normalization, integration-method comparison) is a
future step's scope, not built yet.

`4.1_rna_normalize.py` consumes 3.1's `{sample}.rna_final.h5ad` and runs
two parallel single-sample normalize+UMAP methods, matching the lab's own
notebook's two methods ("I did both sctransform and log method") --
NOT cross-sample/replicate integration, which stays a future step's scope:
  - **log**: scanpy-native `normalize_total -> log1p -> HVG -> PCA -> UMAP
    -> Leiden`. No scaling and no mt%/ribo%/cell-cycle regression --
    matches current scanpy's own tutorial exactly (it dropped both; see
    `configs/normalize_rna.yaml` for the full comparison against
    sc-best-practices and a live scanpy/Seurat disagreement on this).
  - **sct**: `Seurat::SCTransform` via rpy2 (the one piece of math with no
    scanpy equivalent to the notebook's validated method), then
    `RunPCA`/`RunUMAP` in R matching the notebook's own pipe chain exactly,
    transferred back into the same AnnData object. Clustering on this path
    is still done in scanpy (on the transferred `sct.pca`), so both paths
    share one clustering method rather than mixing Leiden and
    Seurat's FindClusters.

Cell-cycle scores are computed once (`sc.tl.score_genes_cell_cycle`, using
`cc.genes` pulled fresh from the R Seurat package so both paths agree on
the same gene set) and stored regardless of whether a path regresses on
them. Writes `{sample}.rna_normalized.h5ad` (both paths' PCA/UMAP/Leiden +
cell-cycle scores in one object), a 4-panel UMAP comparison plot, and a
summary JSON. Confirmed on this lab's real D12_rep1 sample: SCTransform is
too heavy for the login node directly (killed by its watchdog after ~2min)
-- must run through sbatch/srun, matching every other compute step.

`4.2_atac_normalize.py` consumes 3.1's `{sample}.atac_final.h5ad` and runs
SnapATAC2's spectral embedding (`spectral -> umap -> knn -> leiden`),
matching the lab's own actual normalize notebook
(`prepare.peakcall.snapatac2.ipynb` -- not `9.4.D12.atac.snapatac2.ipynb`,
which turned out to be peak-calling-by-cluster only, no spectral/UMAP at
all). NO fragment reloading, no `add_tile_matrix`/`select_features` --
3.1's output already carries the tile matrix and selected features
forward from step 2.2 untouched, so this step only recomputes the
embedding on the clean (doublet-removed, barcode-intersected) population,
overwriting the stale one that rode along from 2.2's own diagnostic pass
(computed there on the pre-doublet-removal population). Two deviations
from the notebook's bare defaults: `umap`/`knn` exclude spectral dim 0
(matches 2.2's own convention), and `leiden`'s `objective_function` is set
explicitly to `"modularity"` (Galaxy/scverse's SnapATAC2 docs warn CPM
gives poorly-separated UMAP clusters -- already this installed
snapatac2's default, so this documents the choice rather than changing
behavior). Writes `{sample}.atac_normalized.h5ad`, a 3-panel UMAP plot
(Leiden, TSS enrichment, fragment count), and a summary JSON. Confirmed
too heavy for the login node directly (killed by its watchdog after ~35s)
-- must run through sbatch/srun, same as 4.1.

## Output layout

```
$OUT_DIR/
├── rna_qc/
│   ├── h5ad/     {sample}.qc_annotated.h5ad, {sample}.qc_filtered.h5ad
│   ├── plots/    {sample}.qc_before_after.png
│   └── summary/  {sample}.qc_summary.csv
├── remove_background/
│   ├── h5ad/     {sample}.soupx_corrected.h5ad
│   ├── plots/    {sample}.soupx_pre_post_umap.png
│   └── summary/  {sample}.soupx_summary.json
├── detect_doublets/
│   ├── h5ad/     {sample}.doublets_annotated.h5ad, {sample}.doublets_removed.h5ad
│   ├── plots/    {sample}.doublet_diagnostics.png
│   └── summary/  {sample}.doublets_summary.csv
├── atac_qc/
│   ├── amulet/   {sample}.amulet_result.csv, {sample}.amulet_score.png
│   ├── h5ad/     {sample}.import_tsse_checkpoint.h5ad, {sample}.atac_qc_annotated.h5ad, {sample}.atac_qc_final.h5ad
│   ├── plots/    {sample}.frag_size_distr.png, {sample}.threshold_review.png,
│   │             {sample}.doublet_breakdown.png, {sample}.umap_doublets_clusters.png
│   └── summary/  {sample}.atac_qc_summary.csv
├── merged/
│   ├── h5ad/     {sample}.rna_final.h5ad, {sample}.atac_final.h5ad
│   └── summary/  {sample}.merge_summary.csv
├── rna_normalize/
│   ├── h5ad/     {sample}.rna_normalized.h5ad
│   ├── plots/    {sample}.rna_normalize_umap.png
│   └── summary/  {sample}.rna_normalize_summary.json
└── atac_normalize/
    ├── h5ad/     {sample}.atac_normalized.h5ad
    ├── plots/    {sample}.atac_normalize_umap.png
    └── summary/  {sample}.atac_normalize_summary.json
```

## Cleaning up intermediate files

Every step above keeps its own h5ad output -- deliberately, so a
decision (the MAD threshold, the doublet union call) can be revisited
without recomputing something expensive. That adds up fast, especially
for ATAC's genome-wide tile matrix. Once you've reviewed the plots and
summaries and are satisfied with a sample's result:

```bash
bash scripts/3.2_cleanup_intermediate_h5ad.sh          # dry run -- lists what would be deleted + total size
bash scripts/3.2_cleanup_intermediate_h5ad.sh --yes    # actually deletes
```

Only touches a sample once its `merged/h5ad/{sample}.rna_final.h5ad` and
`{sample}.atac_final.h5ad` already exist -- skips (does not delete
anything for) a sample that hasn't reached that point yet. Deletes every
step's intermediate h5ad (including `detect_doublets`'s
`doublets_removed.h5ad` and `atac_qc`'s `atac_qc_final.h5ad`, both
superseded by `merged/h5ad`'s versions), but never touches any plot,
summary CSV, AMULET's result CSV, or anything under `merged/`. This is
destructive and not easily reversible -- rerunning from scratch means
redoing everything each deleted file represents, including AMULET, which
is slow -- so it defaults to a dry run and only deletes with `--yes`.

## Before trusting your QC calls

- Look at every sample's before/after plot, not just the summary numbers --
  the MAD default is meant to be universal, but "universal" can still be
  wrong for one weird sample, and that's what `filtering_strategy: manual`
  in `configs/qc_RNA.yaml` is for.
- MT%/Ribo% are computed and logged for every sample regardless of whether
  a cutoff is enforced -- check them even when `ribo_pct_max` is null.
- RNA and ATAC QC run as independent modules and each produce their own
  QC-passing cell list -- `scripts/3.1_merge_modalities.py` reconciles
  them afterward (intersects barcodes, subsets both objects, writes
  final RNA + ATAC h5ad per sample), matching the lab's own notebook's method
  exactly. Check `{sample}.merge_summary.csv`'s retention percentages --
  a low overlap almost always means a barcode format problem, not
  genuinely poor multiome capture; the script warns loudly and prints
  example barcodes from both sides when that happens.
- For ATAC QC, general QC uses fixed thresholds on purpose, not MAD --
  unlike RNA QC, all three references (this lab's notebook, IGVF-UCSD,
  sc-best-practices) use fixed thresholds for ATAC, and MAD's assumption
  (one well-behaved population with rare tail outliers) doesn't hold well
  against a raw ATAC fragment file, where the empty-droplet/background
  population is typically a whole second mode, not a small tail. Every
  run writes `{sample}.threshold_review.png` -- the current config's
  cutoff lines drawn directly on the fragment-count and TSSe distributions
  plus their joint scatter, matching the lab's own notebook's
  plot-then-decide workflow. Look at it per sample before trusting the
  filter: `configs/atac_qc.yaml`'s numbers are IGVF-UCSD's own validated
  defaults, not yet checked against this lab's own data, and are meant to
  be edited once the plot gives a reason to. Use `scripts/qc_only_atac.sh
  <sample_id>` to do that tuning fast, per sample, without paying for the
  full pipeline each iteration -- see "Run it, in order" above.
