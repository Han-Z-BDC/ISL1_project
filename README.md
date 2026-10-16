# ISL1 CUT&RUN

ISL1 CUT&RUN (with matched IgG) across the ISL1 stem-cell islet differentiation, control vs
ISL1-KO, at D12, D15 and D22. Analysis for publication, built on the general lab pipeline
[cutnrun-pipeline](https://github.com/ZhuLab-CUAnschutz/cutnrun-pipeline), vendored unchanged
in `pipeline/` (copied from commit `4845283`, 2026-10-16).

Repo (Alpine): `/pl/active/ZhuLab_BDC/GitHub/ISL1_cutNrun` -- local git, not pushed to GitHub yet.

## Data

Admera, two deliveries under `/pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run/`:

| Batch | Folder | Libraries | md5 |
|---|---|---|---|
| 1 | `7535-0001-08202026_144755` (2026-08-20) | rep1: D12 Ctr, D12 KO, D15 Ctr, D15 KO, D22 Ctr (ISL1 + IgG each, 10) | not yet checked |
| 2 | `7535-0002-10082026_215649` (2026-10-08) | D12 Ctr rep2; D12 KO rep2-4; D15 KO rep3-4; D22 Ctr rep2-3 (16) | 32/32 OK (job 34347701) |

Batch 2 arrived in the same folder as the ISL1 scRNA-seq delivery and was moved here 2026-10-16
(see the ISL1_scRNAseq repo).

Samples per condition (ISL1 + matched IgG each):

| | Ctr | KO |
|---|---|---|
| D12 | rep1, rep2 | rep1, rep2, rep3, rep4 |
| D15 | rep1 | rep1, rep3, rep4 |
| D22 | rep1, rep2, rep3 | -- |

Confirmed 2026-10-16: D15 Ctr has only the one replicate, and there is no D15 KO rep2
(numbering kept as delivered, rep3/rep4).

ISL1 CUT&RUN in the KO doubles as the antibody-specificity control: true ISL1 peaks should be
lost or strongly reduced there.

Open questions on the sample sheet (to confirm with the bench):
- Batch 1 libraries are named `...-3d-1` (D12/D15) or `D22-Ctr-1`; treated here as rep1.
- "Ctr" is `Ctr-GFP` at D12/D15 and plain `Ctr` at D22 -- confirm the same control line.

## Layout

```
pipeline/       cutnrun-pipeline, unchanged (its config.sh and logs/ are symlinks, gitignored)
study_config/   config.sh (all paths), samples.tsv (pipeline sample sheet),
                sample_metadata.tsv (day / genotype / batch per sample), setup_links.sh
scripts/        study-specific steps (md5 check now; downstream analysis later)
```

`samples.tsv` keeps the pipeline's 5-column format; anything extra (day, genotype, batch,
delivery folder) is in `sample_metadata.tsv`, keyed by `sample_name`.

## Setup (once, on Alpine)

1. Build the pipeline env (on a compute node):
   ```
   acompile --ntasks=4 --time=01:00:00
   module load miniforge
   conda env create -f pipeline/environment.yml     # -> envs/cutnrun, = TOOLS in config.sh
   ```
2. Link config, logs and fastqs (both batches into one `RAW_DIR` on scratch), then validate:
   ```
   bash study_config/setup_links.sh
   ```

Reference: hg38 bowtie2 index `/pl/active/ZhuLab_BDC/genome/hg38_bt2/hg38`, blacklist
`/pl/active/ZhuLab_BDC/genome/GRCh38_unified_blacklist.bed`. IgG background: `matched`.

## Pipeline

Step 0 (this repo): `sbatch scripts/00_md5_check.sb 7535-0001-08202026_144755`

Then the general pipeline, submitted from `pipeline/` (its Slurm logs go to `logs/` relative to
the submit directory, which `setup_links.sh` points at scratch):

```
cd pipeline
sbatch scripts/00_fastqc.sb
bash   scripts/submit_align.sh            # array, one task per library (26)
sbatch scripts/02_build_igg_background.sb
sbatch scripts/03_call_macs3.sb
sbatch scripts/04_tracks.sb
```

All outputs on scratch: `/scratch/alpine/hzhu10@xsede.org/ISL1_cutNrun/results/`
(purged after 90 days -- copy keepers to `/pl/active`).

## Local review copy (Mac)

`~/Documents/Claude/ISL1_cutNrun` is a clone (remote `alpine`) holding code plus HTML reports
only -- no BAMs, bigwigs or peak files unless asked for. Refresh with:

```
bash scripts/fetch_reports.sh
```
