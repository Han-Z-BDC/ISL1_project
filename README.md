# ISL1 CUT&RUN

ISL1 CUT&RUN (with matched IgG) across the ISL1 stem-cell islet differentiation, control vs
ISL1-KO, at D12, D15 and D22. Analysis for publication, built on the general lab pipeline
[cutnrun-pipeline](https://github.com/ZhuLab-CUAnschutz/cutnrun-pipeline), vendored unchanged
in `pipeline/` (copied from commit `1afe307`, 2026-10-16).

Repo (Alpine): `/pl/active/ZhuLab_BDC/GitHub/ISL1_cutNrun` -- local git, not pushed to GitHub yet.

## Data

Admera, two deliveries under `/pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run/`:

| Batch | Folder | Libraries | md5 |
|---|---|---|---|
| 1 | `7535-0001-08202026_144755` (2026-08-20) | rep1: D12 Ctr, D12 KO, D15 Ctr, D15 KO, D22 Ctr (ISL1 + IgG each, 10) | 20/20 OK (job 34348706) |
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
`/pl/active/ZhuLab_BDC/genome/GRCh38_unified_blacklist.bed`.

## Pipeline

Step 0 (this repo): `sbatch scripts/00_md5_check.sb <batch folder>` (both batches OK, 52/52),
`sbatch scripts/00b_build_env.sb`.

Then the general pipeline, submitted from `pipeline/` (the scripts find `config.sh` and `logs/`
through `SLURM_SUBMIT_DIR`; `setup_links.sh` points `logs/` at scratch):

```
cd pipeline
sbatch scripts/00_fastqc.sb
sbatch --array=0-25 --cpus-per-task=8 --mem=30G --time=06:00:00 scripts/01_align.sb   # see Resources
sbatch --mem=8G --time=04:00:00 scripts/02_build_igg_background.sb                    # own IgG
sbatch --mem=8G --time=04:00:00 --export=ALL,IGG_METHOD=pooled scripts/02_build_igg_background.sb  # pool of 13
sbatch --mem=12G scripts/03_call_macs3.sb
sbatch --time=12:00:00 --export=ALL,IGG_METHOD=pooled scripts/03_call_macs3.sb
sbatch scripts/04_tracks.sb
cd .. && sbatch scripts/03b_compare_igg_methods.sb
```

### Resources

Requests are overridden on the `sbatch` command line rather than edited in `pipeline/`, so the
vendored copy stays identical. Alignment: 8 CPU, 16 GB, 6 h per task instead of the pipeline's
32 GB / 24 h. On `acpu` each core carries 3.84 GB, so 32 GB is charged as 9 cores; the real peak is
about 7 GB (bowtie2 with hg38, then `samtools sort -@ 8`). The March 2026 lab CUT&RUN run
(job 24475155) took about 46 min per 31 M-pair library on 12 cores, so 6 h leaves a wide margin.
First submission (job 34353064, pipeline defaults) was cancelled while still pending;
the alignment that ran is job 34354886: all 26 tasks completed in 1.6-3.0 h. Peak RSS was
10-16 GB, higher than the estimate above; the largest task (index 18) was within 1.5 MB of the
16 GB limit. For a rerun, request `--mem=30G`: 8 cores already carry 30.7 GB on `acpu`,
so the charge stays at 8 cores.

Later steps (job IDs in the run log below): 02 at 4 CPU / 8 GB (samtools view only),
03 at 2 CPU / 12 GB matched and 2 CPU / 16 GB pooled (MACS3 is single-threaded; the pooled
control holds 117 M fragments), 04 at the pipeline defaults.

### Alignment (job 34354886)

Bowtie2 overall alignment 90-98% for every library; 23-49 M read pairs. After the MAPQ 20 /
proper-pair filter, ISL1 libraries keep 18.7-36.6 M fragments (duplicate rate 9.5-41.6%).
Batch-1 IgGs are the most duplicated and the shallowest once deduplicated:

| IgG | Batch | Dup rate | Unique fragments |
|---|---|---|---|
| D22_Ctr_rep1 | 1 | 58.0% | 10.1 M |
| D12_KO_rep1 | 1 | 54.7% | 9.4 M |
| D15_Ctr_rep1 | 1 | 42.4% | 13.9 M |
| D12_Ctr_rep1 | 1 | 36.5% | 13.5 M |
| D15_KO_rep1 | 1 | 33.1% | 16.3 M |
| batch-2 IgGs (8) | 2 | 9.0-16.4% | 20.6-34.8 M |

`POOL_TARGET_FRAGMENTS` = 9,000,000: just below the smallest IgG, D12_KO_rep1, which has
9,350,997 unique fragments.

### IgG background: matched vs pooled

Every ISL1 sample is called twice: against its own IgG (`matched`) and against one pool of all
13 IgGs, each deduplicated and downsampled to `POOL_TARGET_FRAGMENTS` with a fixed seed
(`pooled`, the pipeline's own method). `IGG_METHOD` (default `matched`) picks the output tree:
`results/` or `results_pooled/`, the latter sharing `results/bam` through a symlink.

`scripts/03b_compare_igg_methods.sb` then compares the two peak sets per sample (overlap both
ways, top-1000 recovery, bp Jaccard) and flags samples that agree much less than the cohort.
A flagged sample points at its own IgG. Report: `results/igg_method_comparison/`.

All outputs on scratch: `/scratch/alpine/hzhu10@xsede.org/ISL1_cutNrun/results/` and
`results_pooled/` (purged after 90 days -- copy keepers to `/pl/active`).

## Local review copy (Mac)

`~/Documents/Claude/ISL1_cutNrun` is a clone (remote `alpine`) holding code plus HTML reports
only -- no BAMs, bigwigs or peak files unless asked for. Refresh with:

```
bash scripts/fetch_reports.sh
```
