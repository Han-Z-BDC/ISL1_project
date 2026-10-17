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
sbatch scripts/00c_build_idr_env.sb                      # IDR env (python 3.9 / numpy 1.23)
sbatch scripts/02c_build_group_igg.sb                    # final calls: IgG pooled per group
sbatch --array=0-12 scripts/03c_call_macs3_group.sb      # each ISL1 replicate vs its group IgG
sbatch scripts/04c_idr.sb                                # reproducible peaks per group
```

### Final peak calls: group IgG

Each ISL1 replicate is called against its biological group's IgG background: the day x
genotype group's deduplicated IgGs merged, with no downsampling. D12_KO_rep1_IgG is excluded by
the pre-specified matched-vs-pooled criterion (below); it is also the shallowest, least complex
IgG. D12_KO_rep1_ISL1 is therefore called against the other three D12 KO IgGs. D15 Ctr has one
IgG, so its group background is that IgG alone (13.9 M). Settings in `study_config/config.sh`
(`GROUP_IGG_EXCLUDE`); output in `results_group/` with `igg_background/group_igg_manifest.tsv`.

| Group | IgGs used | Approx. unique fragments |
|---|---|---|
| D12_Ctr | rep1, rep2 | 42.7 M |
| D12_KO | rep2, rep3, rep4 (rep1 excluded) | 67.1 M |
| D15_Ctr | rep1 | 13.9 M |
| D15_KO | rep1, rep3, rep4 | 83.2 M |
| D22_Ctr | rep1, rep2, rep3 | 71.3 M |

Each replicate is called twice against its group IgG:
- **`macs3/`:** q 0.01, the pipeline's settings, for browsing and per-sample counts.
- **`macs3_relaxed/`:** p 0.01, one row per peak, top 300k by p value, blacklist filtered.
  This is the input for IDR.

Reproducible peaks (`scripts/04c_idr.sb`, IDR 2.0.4.2, ranked by p value, global IDR < 0.05):
every replicate pair in a group is run, and the group set is the pair with the most passing
peaks (ENCODE's choice for more than two replicates). D15 Ctr has a single replicate, so IDR
runs on two pseudoreplicates (fragments split by read-name parity, pairs kept together).
Its set shows self-consistency only and should be reported as unreplicated.
Output: `results_group/idr/<group>.reproducible.narrowPeak` and `idr_summary.tsv`.

Run 2026-10-17: IDR env 34433836, group IgG 34433838 (13 min), MACS3 array 34433840
(25-54 min per task, peak RSS 12.0 GB, at the 12 GB request then; the script now asks 16G), IDR 34433842.

| Group | q 0.01 peaks per replicate | IDR pairs (passing peaks) | Reproducible set |
|---|---|---|---|
| D12_Ctr | rep1 12,726; rep2 9,996 | rep1-rep2 2,867 | 2,867 |
| D12_KO | rep1 9,815; rep2-4 599-912 | with rep1 5-11; rep2-4 pairs 58-82 | 82 |
| D15_Ctr | rep1 32,748 | pseudoreplicates 9,344 | 9,344 (unreplicated) |
| D15_KO | rep1 47,154; rep3/4 531, 435 | 48-77 | 77 |
| D22_Ctr | rep1 8,602; rep2 35,726; rep3 18,680 | rep1-rep2 1,282; rep1-rep3 1,550; rep2-rep3 8,557 | 8,557 |

**Batch 1 (all rep1) reproduces poorly against batch 2.** Every batch-1 ISL1 library fills the
300k cap of the relaxed calls; batch-2 libraries reach only 120-184k. That means much more weak
signal at p 0.01. In the q 0.01 calls, 79-87% of a rep1's top 1000 peaks are in the batch-2
replicate, but only 59-61% of D22 rep2/rep3's top 1000 are in D22 rep1. The two batch-2 D22
replicates recover 99-100% of each other's top 1000. So rep1 shares the strongest sites but
is noisier and misses many strong ones. D12 Ctr has no batch-2 pair, so its set (2,867) is
limited by rep1 (but see the ENCODE ratios below: D12 rep1 itself is fine). The noisy KO libraries (D12_KO_rep1, D15_KO_rep1; about 10 k and 47 k peaks
against the group IgG) are removed by IDR: 5-11 and 72-77 peaks pass.

### Replicate QC report and ENCODE IDR ratios

`scripts/03d_pseudoreps.sb` (self: job 34471362, pooled: 34471363), `scripts/04d_idr_qc.sb`
(34471364), `scripts/05_replicate_peak_report.sb` (34471365). Report:
`results/replicate_peak_report/replicate_peak_report.html`, with per-replicate stats, FRiP,
Venn/UpSet overlaps, pairwise recovery, IDR and a cross-replicate Jaccard map. Pooled
pseudoreplicate tasks peaked at 5.3 GB, so the 32G override was unnecessary; 16G is enough.

| Group | Self-pseudorep IDR peaks per replicate | Self-consistency | Nt | Np | Rescue | ENCODE |
|---|---|---|---|---|---|---|
| D12_Ctr | rep1 2,641; rep2 2,416 | 1.09 | 2,867 | 4,498 | 1.57 | pass |
| D22_Ctr | rep1 904; rep2 7,919; rep3 4,512 | 8.76 | 8,557 | 10,397 | 1.22 | borderline |
| D15_Ctr | rep1 9,344 | -- | -- | -- | -- | unreplicated |
| D12_KO | 0-17 | inf | 82 | 180 | 2.20 | fail (not meaningful: no binding) |
| D15_KO | 0-17 | inf | 77 | 252 | 3.27 | fail (not meaningful: no binding) |

- **Batch 1 is not uniformly weak.** D12 Ctr rep1 is as self-consistent as rep2, and D12 Ctr
  passes. D15 Ctr rep1 (batch 1) has the highest self-pseudorep count of any replicate (9,344),
  and FRiP 5.0%. Its peaks overlap D22 rep2/rep3 more than any other pair across days
  (Jaccard 0.28-0.29).
- **The weak control replicate is D22 Ctr rep1.** Self-pseudorep 904, against 4,512-7,919 for
  rep2/rep3. FRiP 0.9%, the lowest of the control replicates. It misses up to 31% of its
  partners' top 1000 peaks. It alone makes D22 borderline; the rescue ratio (1.22) is fine.
- **D12 Ctr rep1 has low FRiP (1.5%) but holds up on IDR,** so it is not a failing replicate.
- **KO groups:** near-zero self-pseudorep counts are what a specificity control should show.

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

Later steps, all completed:

| Step | Job | Request | Elapsed | Peak RSS |
|---|---|---|---|---|
| 02 matched | 34358598 | 4 CPU / 8 GB | 1:30 | <0.1 GB |
| 02 pooled | 34358599 | 4 CPU / 8 GB | 2:06 | <0.1 GB |
| 03 matched | 34358600 | 2 CPU / 12 GB | 2:05 | 8.3 GB |
| 03 pooled | 34358601 | 2 CPU / 16 GB, 12 h | 6:58 | 6.7 GB |
| 04 tracks | 34358602 | 8 CPU / 24 GB (defaults) | 2:50 | 3.6 GB |
| 03b compare | 34388152 | 1 CPU / 4 GB | <1 min | -- |

02 could run at 1 CPU / 2 GB.

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
ways, top-1000 recovery, bp Jaccard). It flags a sample when under 10% of one peak set is in the
other, or when its pooled/matched peak-count ratio is a cohort outlier (>3 scaled MADs on log2).
Report: `results/igg_method_comparison/`.

Result (2026-10-17):

- **12 of 13 samples agree.** In Ctr samples, 93-100% of the top 1000 peaks are recovered
  either way. The pool calls 1.0-3.4x more peaks for every sample (median log2 ratio 0.60).
- **One flag: D12_KO_rep1.** 28 peaks vs own IgG, 10,803 vs the pool; 0.15% of pooled peaks
  are in the matched set.
- **The flag is a depth effect, not the IgG's content** (`scripts/diag/depth_matched_d12ko1.sb`,
  job 34388231). Against other controls downsampled to its own IgG's 9.35 M fragments,
  D12_KO_rep1 gives 16 peaks (D12_KO_rep2 IgG) and 37 peaks (pool). MACS3 scales the larger
  file down (`--scale-to small`), so a control shallower than the ISL1 library (9.4 M vs
  27.0 M) lowers sensitivity. The matched calls for the five shallow batch-1 IgGs are
  therefore conservative.
- **The 10.8 k pooled peaks are not ISL1 binding.** Only 2% are at Ctr ISL1 peaks (D12 Ctr
  union); other D12 KOs are at 20-27%. They are narrow (median 164 bp), genome-wide in
  proportion to chromosome size (not a CNV), and weak: about 16 fragments each, 0.6% of
  fragments. That points to a noisier ISL1 library, not ISL1 signal.
- **D15_KO_rep1 is a separate ISL1-library problem.** 35,123 peaks even against its own IgG,
  vs about 450 for D15 KO rep3/4, and only 7% at D15 Ctr ISL1 sites. Same narrow,
  genome-wide pattern; 41.6% duplicates, the lowest-depth ISL1 library (10.9 M unique). Not
  flagged, because both IgG methods agree on it.
- **Duplicates do not explain either** (`scripts/diag/dup_in_peaks.sb`, job 34387009). In-peak
  duplicate share is 1.2-1.6x the genome-wide rate in every sample, controls included.
- **KO residual peaks (D12 KO rep2-4).** 600-1,200 peaks hold 3-10% of fragments at 58-73%
  duplicates: high-signal artifact regions the blacklist misses. Handle before differential
  analysis.
- **Contigs.** The hg38 index includes unplaced and alt contigs (`GL000195.1`, `KI270438.1`).
  Publication peak sets should probably keep chr1-22, X and Y only.

All outputs on scratch: `/scratch/alpine/hzhu10@xsede.org/ISL1_cutNrun/results/` and
`results_pooled/` (purged after 90 days -- copy keepers to `/pl/active`).

## Local review copy (Mac)

`~/Documents/Claude/ISL1_cutNrun` is a clone (remote `alpine`) holding code plus HTML reports
only -- no BAMs, bigwigs or peak files unless asked for. Refresh with:

```
bash scripts/fetch_reports.sh
```
