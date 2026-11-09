# ISL1 KO / H3 control CUT&RUN time course

CUT&RUN for ISL1, Flag, H3K27ac and H3K27me3 (each with an IgG control) across islet
differentiation (D9, D15, D18, D22) in the H3 control line and the ISL1 KO line. Analyzed
on its own, independent of the earlier ISL1 CUT&RUN batches (repo `ISL1_cutNrun`). Built on
the general lab pipeline
[cutnrun-pipeline](https://github.com/ZhuLab-CUAnschutz/cutnrun-pipeline), vendored unchanged
in `pipeline/` (copied from commit `1afe307`).

Repo (Alpine): `/pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun` -- local git, not pushed
to GitHub yet.

## Data

One Admera delivery: `/pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run/7535-0003-11022026_200951/`.
It was downloaded 2026-11-05 with a copy of `jobs/templates/transfer/admera.download.sb` on
`acompile`. All 130 delivered files match Admera's md5.txt (`md5_check.txt`). The 124 CUT&RUN
fastqs were moved there from the scRNA folder; `md5.txt` lists them and
`md5.admera_full_delivery.txt` is the full delivery list.

62 libraries in 18 blocks (day x line x replicate); each block has its own IgG:

| Line | Days | Replicates | Antibodies |
|---|---|---|---|
| H3CTR (H3 control) | D9, D15, D18, D22 | 1, 2 | ISL1, H3K27ac, H3K27me3, IgG |
| ISL1KO | D9, D15, D18, D22 | 1, 2 | H3K27ac, H3K27me3, IgG |
| CtrGFP | D15 | 2 | ISL1, Flag, IgG |
| KOISL1 | D15 | 2 | ISL1, Flag, IgG |

- **D15 CtrGFP and D15 KOISL1 are set aside (2026-11-06).** Their 6 libraries (ISL1, Flag,
  IgG) are a separate question and stay out of the analysis for now. Their peaks were called
  with the rest and are left in `results/`; downstream steps use the H3CTR/ISL1KO time course
  only.
- **No ISL1 rescue in this batch** (confirmed 2026-11-05). The Flag target is still to be
  confirmed.
- **Antibody names normalized** from two delivery typos: `H327Kac` -> H3K27ac (D15 ISL1-KO),
  `H3K27Me` -> H3K27me3 (D22 H3-CTR).
- **Sample names** are `<day>_<line>_rep<n>_<antibody>`, built by `study_config/make_samples.py`.

## Layout

```
pipeline/       cutnrun-pipeline, unchanged (its config.sh and logs/ are symlinks, gitignored)
study_config/   config.sh (all paths), make_samples.py -> samples.tsv (pipeline sample sheet)
                and sample_metadata.tsv (day / line / antibody), setup_links.sh
scripts/        study-specific steps; fetch_reports.sh (Mac)
```

## Setup (once, on Alpine)

```
python3 study_config/make_samples.py /pl/active/ZhuLab_BDC/raw_data/ISL1_cut_N_run/7535-0003-11022026_200951 \
  7535-0003-11022026_200951 study_config/samples.tsv study_config/sample_metadata.tsv
bash study_config/setup_links.sh
```

## Pipeline

Submitted from `pipeline/` (the scripts find `config.sh` and `logs/` through `SLURM_SUBMIT_DIR`).
All analysis jobs run on `acpu`. `acompile` is kept for data downloads only (8 cores per user).
The first FastQC run (job 35274348) went to `acompile` before that rule was set.

```
cd pipeline
sbatch scripts/00_fastqc.sb
sbatch --array=0-61 --cpus-per-task=8 --mem=30G --time=06:00:00 scripts/01_align.sb
```

### Alignment QC (2026-11-06)

FastQC/MultiQC (job 35274348) and alignment (array 35274349, 62 tasks) all completed. Peak
memory was 17.4 GB and the longest task took 3 h 29 min, against 30 GB / 6 h requested. Stats
come from `results/qc/alignment/` (cutadapt, bowtie2, flagstat):

| Antibody | Libraries | Align % | Dup % | Unique fragments (M), range (median) |
|---|---|---|---|---|
| ISL1 | 10 | 88.0-95.0 | 16.1-31.2 | 6.8-31.3 (25.5) |
| Flag | 2 | 91.8-93.5 | 18.8-23.5 | 29.0-30.5 |
| H3K27ac | 16 | 92.7-97.5 | 12.4-18.8 | 19.2-40.1 (25.9) |
| H3K27me3 | 16 | 89.0-97.1 | 13.0-29.1 | 17.1-39.1 (29.7) |
| IgG | 18 | 85.4-96.0 | 15.9-22.4 | 18.0-37.8 (25.8) |

- **One shallow library:** D22_H3CTR_rep1_ISL1 has 9.8 M read pairs and 6.8 M unique fragments.
  Every other library has at least 17 M unique fragments.

### Peak calling (ENCODE-style, matched IgG)

Each TF row is called against its own block IgG, deduplicated by the pipeline's step 02
(`IGG_METHOD=matched`). Target BAMs keep their duplicates, as in the pipeline.

```
cd pipeline
sbatch scripts/02_build_igg_background.sb
cd ..
sbatch --array=0-43 --dependency=afterok:<02 job> scripts/03_call_peaks.sb
```

- **ISL1, Flag, H3K27ac:** `results/macs3/`. MACS3 narrow, q 0.01, `--call-summits`, the
  pipeline's command, blacklist filtered.
- **H3K27me3:** `results/macs3_broad/`. MACS3 `--broad`, q 0.01, broad cutoff 0.1, blacklist
  filtered.
- **ISL1 and Flag, for IDR:** `results/macs3_relaxed/`. p 0.01, top 300k peaks. The
  single-replicate D15 CtrGFP and D15 KOISL1 groups also get two pseudoreplicates.
- **Per-sample stats:** `results/qc/peaks/<sample>.peak_stats.tsv` (peaks and FRiP).

Results (2026-11-06; jobs 35295284 and 35295285, all completed; at most 9.1 GB and 1 h 15 min
per task). Regions are distinct blacklist-filtered peaks, and FRiP counts all filtered fragments:

| Antibody | Mode | Regions per library | FRiP |
|---|---|---|---|
| H3K27me3 | broad | 54.8k-119.2k | 33-49% |
| H3K27ac | narrow | D15/D22: 32.7k-80.9k; D9/D18: 7.2k-42.2k | D15/D22: 2.1-8.3%; D9/D18: 0.3-3.0% |
| ISL1 | narrow | 18.4k-37.8k | 1.2-4.6% |
| Flag | narrow | 12.1k, 12.4k | 1.3% |

ISL1 checks (top 1000 peaks by p value):
- **Replicate pairs agree poorly.** 20-40% of rep1's top 500 peaks fall in any rep2 peak (D15,
  D22).
- **Match to the earlier batches' D15 Ctr conservative ISL1 set (9,329 peaks):**
  - Strong: D15 CtrGFP 86%, D18 H3CTR rep1 74%, D22 H3CTR rep2 60%.
  - Weaker: D22 rep1 26%, D18 rep2 25%, D15 KOISL1 21%, D15 H3CTR rep1/rep2 16% / 8%.
  - D9 (both replicates): 3%.

  The comparison is a QC reference only; the analysis stays standalone.
- **Flag peaks do not match ISL1** (0.2-0.7% vs the earlier set, 1-2% vs ISL1 in the same
  block). They look like background.

### Reproducible peaks (ENCODE design)

H3CTR and ISL1KO groups only (day x line x antibody: 4 ISL1, 8 H3K27ac, 8 H3K27me3); D15 CtrGFP
and KOISL1 are set aside. Target duplicates are kept throughout, as is usual for CUT&RUN (ENCODE
removes them for ChIP-seq).

```
sbatch --array=0-39 scripts/04_relaxed_pseudoreps.sb
sbatch --array=0-19 --dependency=afterok:<04 job> scripts/05_reproducible_peaks.sb
```

- **04, per library:** relaxed MACS3 call vs its own IgG (p 0.01; top 300k for ISL1, top 500k
  for histones; blacklist filtered). The library is then split in half by read name, and each
  half is called the same way.
- **05, per group:**
  - **Pooled calls:** both replicates vs both IgGs, plus the two pooled halves.
  - **Comparisons (each keeps the reference peaks supported by two lists):**
    - Nt: pooled peaks vs rep1 and rep2 (replicated).
    - Np: pooled peaks vs the pooled halves (pseudoreplicated).
    - N1 / N2: each replicate vs its own halves (QC only).
  - **Methods:**
    - Overlap (all marks, ENCODE histone rule): at least 50% of either peak.
    - IDR (ISL1, ENCODE TF): IDR < 0.05, pooled peaks as `--peak-list`. The earlier ISL1
      batches ran IDR without a peak list.
  - **Not done here:** no automatic conservative/optimal pick; both the Nt and Np sets are kept.
  - **Ratios:** self-consistency N1/N2 and rescue Np/Nt; both at most 2 = pass.
- **Output:** `results/reproducibility/{overlap,idr}/<group>.{Nt,Np}.narrowPeak` and
  `results/reproducibility/qc/<group>.repro_qc.tsv`.

Results (2026-11-06; jobs 35309994 and 35309995, all 60 tasks completed). Step 04 used at most
7.6 GB and 1 h 44 min per task, and step 05 at most 14.3 GB and 1 h 33 min.

- **Relaxed lists:** ISL1 215k-281k peaks per library, H3K27ac 183k-269k, H3K27me3 90k-182k. No
  list reached its cap (300k / 500k), so the caps never applied.

Reproducible sets. Self = N1 / N2 ratio and rescue = Np / Nt ratio. FRiP is for rep1 / rep2 in
Nt:

| Group | Method | N1 / N2 | Nt | Np | Self | Rescue | Flag | Nt Mb | Nt FRiP |
|---|---|---|---|---|---|---|---|---|---|
| D9_H3CTR_ISL1 | overlap | 23,676 / 21,282 | 18,792 | 28,684 | 1.11 | 1.53 | pass | 12.4 | 1.7 / 1.4% |
| D9_H3CTR_ISL1 | IDR | 0 / 7 | 11 | 33 | inf | 3.00 | fail | 0.0 | 0.7 / 0.4% |
| D15_H3CTR_ISL1 | overlap | 20,540 / 34,334 | 24,273 | 40,649 | 1.67 | 1.67 | pass | 17.0 | 2.8 / 2.7% |
| D15_H3CTR_ISL1 | IDR | 44 / 17 | 85 | 159 | 2.59 | 1.87 | borderline | 0.1 | 1.3 / 1.2% |
| D18_H3CTR_ISL1 | overlap | 27,900 / 22,759 | 22,144 | 22,094 | 1.23 | 1.00 | pass | 15.5 | 3.2 / 2.4% |
| D18_H3CTR_ISL1 | IDR | 855 / 92 | 485 | 965 | 9.29 | 1.99 | borderline | 0.6 | 1.7 / 1.1% |
| D22_H3CTR_ISL1 | overlap | 25,826 / 31,838 | 23,509 | 33,081 | 1.23 | 1.41 | pass | 17.4 | 3.2 / 3.2% |
| D22_H3CTR_ISL1 | IDR | 123 / 474 | 597 | 938 | 3.85 | 1.57 | borderline | 0.7 | 1.0 / 1.6% |
| D9_H3CTR_H3K27ac | overlap | 21,517 / 17,760 | 20,427 | 25,806 | 1.21 | 1.26 | pass | 14.4 | 1.3 / 1.3% |
| D9_ISL1KO_H3K27ac | overlap | 19,537 / 17,865 | 13,773 | 22,796 | 1.09 | 1.66 | pass | 8.7 | 0.7 / 0.7% |
| D15_H3CTR_H3K27ac | overlap | 53,230 / 56,334 | 71,628 | 78,798 | 1.06 | 1.10 | pass | 97.2 | 10.8 / 11.4% |
| D15_ISL1KO_H3K27ac | overlap | 50,577 / 29,545 | 50,164 | 81,025 | 1.71 | 1.62 | pass | 59.2 | 6.5 / 5.4% |
| D18_H3CTR_H3K27ac | overlap | 42,960 / 22,450 | 38,764 | 63,940 | 1.91 | 1.65 | pass | 37.4 | 4.0 / 3.1% |
| D18_ISL1KO_H3K27ac | overlap | 16,587 / 22,111 | 27,103 | 34,376 | 1.33 | 1.27 | pass | 20.7 | 2.2 / 2.0% |
| D22_H3CTR_H3K27ac | overlap | 55,046 / 42,029 | 58,227 | 67,493 | 1.31 | 1.16 | pass | 74.3 | 10.5 / 7.6% |
| D22_ISL1KO_H3K27ac | overlap | 45,092 / 58,363 | 69,599 | 84,060 | 1.29 | 1.21 | pass | 94.6 | 11.8 / 10.8% |
| D9_H3CTR_H3K27me3 | overlap | 32,089 / 27,258 | 27,476 | 34,004 | 1.18 | 1.24 | pass | 99.6 | 29.5 / 42.6% |
| D9_ISL1KO_H3K27me3 | overlap | 27,313 / 25,412 | 26,450 | 27,169 | 1.07 | 1.03 | pass | 88.9 | 31.1 / 42.2% |
| D15_H3CTR_H3K27me3 | overlap | 65,551 / 47,095 | 56,610 | 56,812 | 1.39 | 1.00 | pass | 169.3 | 34.0 / 45.0% |
| D15_ISL1KO_H3K27me3 | overlap | 51,879 / 51,785 | 56,069 | 62,487 | 1.00 | 1.11 | pass | 166.2 | 39.2 / 42.5% |
| D18_H3CTR_H3K27me3 | overlap | 55,698 / 49,899 | 53,900 | 55,914 | 1.12 | 1.04 | pass | 167.8 | 38.5 / 45.3% |
| D18_ISL1KO_H3K27me3 | overlap | 55,128 / 59,776 | 60,960 | 62,742 | 1.08 | 1.03 | pass | 170.6 | 35.4 / 36.1% |
| D22_H3CTR_H3K27me3 | overlap | 67,488 / 53,890 | 63,028 | 65,196 | 1.25 | 1.03 | pass | 181.8 | 39.2 / 34.2% |
| D22_ISL1KO_H3K27me3 | overlap | 60,239 / 65,139 | 64,708 | 79,556 | 1.08 | 1.23 | pass | 208.8 | 39.5 / 39.3% |

- **Histone marks:** all 16 groups pass by overlap, and Nt and Np are close. H3K27ac is strong at
  D15 and D22 and weak at D9 and D18, as in the per-sample calls.
- **ISL1, overlap:** all four groups pass, with 18.8k-24.3k Nt peaks.
- **ISL1, IDR:** few peaks pass, from 11 Nt peaks (D9) to 597 (D22). For D18 Nt, IDR matched
  27,835 pooled peaks in both replicates; the fit puts about 6% of them in the reproducible
  component, and 485 pass IDR < 0.05. The replicates rank peaks differently, which fits the poor
  top-peak agreement in the per-sample checks.
- **chrM:** 0-6 chrM peaks are in each set; the blacklist has no chrM entries. They are not
  removed yet.
- **Histone marks: Nt chosen (2026-11-07).**

### ISL1 consensus peaks (nf-core/cutandrun style, 2026-11-07)

Standard CUT&RUN practice is to call each replicate, merge, and keep peaks found in both. That is
nf-core/cutandrun with `--replicate_threshold 2`.

```
sbatch --array=0-7 scripts/06_isl1_seacr.sb
sbatch --dependency=afterok:<06 job> scripts/07_isl1_consensus.sb
```

- **A_macs3:** step 03 MACS3 calls (q 0.01). Summit rows are collapsed to regions; rep1 and rep2
  are merged (overlapping or book-ended), and regions with a peak from both replicates are kept.
- **A_seacr:** SEACR 1.3 `norm stringent` against the deduplicated block IgG. Fragment bedGraphs
  use fragments up to 1000 bp, with target duplicates kept. Same consensus rule as A_macs3.
- **B / C:** the step 05 overlap and IDR sets, copied with chrM removed.
- **chrM dropped** from every set. Output: `results/isl1_consensus/`, with `compare.tsv`.
- **Resources:** jobs 35344508 and 35344509 completed; at most 4.5 GB and 39 min per task.

Peaks per set. "Ref" is the share of the top 1000 peaks that overlap the earlier batches' D15 Ctr
conservative ISL1 set (QC reference only; shuffled peaks match 0-0.8%):

| Group | A_macs3 (rep1 / rep2 / both) | Ref | B Nt | Ref | C IDR Nt / Np | Ref (Np) |
|---|---|---|---|---|---|---|
| D9 | 18,937 / 18,734 / 157 | 13% | 18,788 | 6% | 7 / 29 | 14% |
| D15 | 21,797 / 37,831 / 530 | 34% | 24,268 | 33% | 80 / 154 | 63% |
| D18 | 25,400 / 18,355 / 976 | 70% | 22,139 | 81% | 480 / 960 | 81% |
| D22 | 36,142 / 35,817 / 1,500 | 74% | 23,505 | 81% | 593 / 933 | 80% |

- **SEACR fails on half the libraries.** SEACR's IgG-based threshold gave 0 peaks for D9 rep1,
  D15 rep2 and D22 rep2 (empirical FDR 0.21-0.30). D22 rep1 (the shallow library) gave 23k peaks
  at an FDR of 0.46. Only D18 has SEACR peaks in both replicates: 3,913 and 2,659, with 482 in
  both.
- **Few peaks are shared by the replicates.** At any cutoff, the both-replicate count is only
  1.2-2.3 times what shuffled rep2 peaks give. At p 0.01 about 70% of shared peaks would be
  chance. That suggests about 1-3k truly shared sites at D18 and D22.
- **B's deeper peaks.** For D18 Nt, the share matching the reference by rank is: top 1000 81%,
  1k-5k 33%, 5k-10k 11%, beyond 10k 5% (chance 0.4%).
- **D9 and D15 H3CTR ISL1 are weak** against the reference at every set.
- **Unplaced contigs:** a few peaks sit on unplaced contigs (e.g. `KI270442.1`); not removed yet.

SEACR relaxed (2026-11-08; jobs 35433987 and 35433988, completed; at most 1.7 GB and 3.5 min
per task). The bedGraphs were reused and the same consensus rule applied (scripts/06 with
`SEACR_MODE=relaxed`, then scripts/08). "Chance" is the both-replicate count after shuffling
rep2's peaks. Consensus Mb, median width and the match to the reference are given in the rows
below:

| Group | Method | Rep1 | Rep2 | Union | Both | Chance | Both / chance | Median width | Ref top 1000 / all |
|---|---|---|---|---|---|---|---|---|---|
| D9 | MACS3 | 18,937 | 18,734 | 37,513 | 157 | 99 | 1.6 | 578 | 13% / 13% |
| D9 | SEACR stringent | 0 | 4,641 | 4,641 | 0 | - | - | - | - |
| D9 | SEACR relaxed | 0 | 38,207 | 38,207 | 0 | - | - | - | - |
| D15 | MACS3 | 21,797 | 37,831 | 59,093 | 530 | 215 | 2.5 | 616 | 34% / 34% |
| D15 | SEACR stringent | 1,440 | 0 | 1,440 | 0 | - | - | - | - |
| D15 | SEACR relaxed | 9,146 | 0 | 9,146 | 0 | - | - | - | - |
| D18 | MACS3 | 25,400 | 18,355 | 42,768 | 976 | 134 | 7.3 | 717 | 70% / 70% |
| D18 | SEACR stringent | 3,913 | 2,659 | 6,086 | 482 | - | - | 6,063 | 71% / - |
| D18 | SEACR relaxed | 18,880 | 22,463 | 37,551 | 3,644 | 1,140 | 3.2 | 5,336 | 66% / 46% |
| D22 | MACS3 | 36,142 | 35,817 | 70,448 | 1,500 | 392 | 3.8 | 701 | 74% / 60% |
| D22 | SEACR stringent | 22,953 | 0 | 22,953 | 0 | - | - | - | - |
| D22 | SEACR relaxed | 22,953 | 0 | 22,953 | 0 | - | - | - | - |

- **Relaxed mode keeps the same three libraries at 0 peaks** (D9 rep1, D15 rep2, D22 rep2;
  SEACR FDR 0.26-0.44). Their IgG-overlap step removes every target region. D22 rep1 is
  unchanged (FDR 0.46).
- **Only D18 gains.** It has 3,644 SEACR relaxed peaks in both replicates against 1,140 by
  chance. The regions are broad (median 5.3 kb), and 46% of them overlap the reference set.
- **The MACS3 consensus is mostly above chance:** 7.3x at D18 and 3.8x at D22. Here the chance
  baseline uses the step 03 per-replicate calls, not the looser step 04 lists tested earlier.

### ISL1 replicate correlation and PCA (2026-11-08)

`sbatch scripts/09_isl1_replicate_correlation.sb` (job 35434729, completed in 29 min, 0.7 GB).

- **Peaks:** step 03 MACS3 calls, chr1-22/X/Y. Two kinds of union: rep1 + rep2 per day
  (37.5k-70.4k), and all 8 ISL1 replicates (200,232).
- **Counts:** fragments per peak with featureCounts (`-p --countReadPairs`). Duplicates are kept,
  and the 8 block IgGs are counted the same way.
- **Analysis:** log2(CPM + 1); Pearson and Spearman correlation; PCA (centered).
- **Report:** `results/qc/replicate_correlation/isl1/report.html`.

| Day | Union peaks | rep1 vs rep2 Pearson | Spearman | rep vs own IgG (Pearson) | Fragments in peaks |
|---|---|---|---|---|---|
| D9 | 37,500 | -0.45 | -0.53 | 0.16 / 0.05 | 1.6 / 1.5% |
| D15 | 59,065 | -0.35 | -0.44 | 0.08 / 0.07 | 2.3 / 2.7% |
| D18 | 42,756 | 0.13 | -0.27 | 0.25 / 0.17 | 2.4 / 1.9% |
| D22 | 70,412 | -0.18 | -0.45 | 0.11 / 0.19 | 5.1 / 3.5% |

- **Negative within-day correlation is a selection effect.** Most union peaks are called in only
  one replicate, so they are high there and low in the other. Only D18 and D22 have a tail of
  strong peaks on the diagonal.
- **All-sample union:** Pearson is 0.11-0.44 between ISL1 libraries, about the same as against
  the IgGs. The pairs are D9 0.24, D15 0.21, D18 0.41, D22 0.25. D18 rep1 and D22 rep2 correlate
  best (0.44).
- **PCA follows signal strength, not day.** On the top 5,000 variable peaks with IgG included:
  - The IgGs cluster together, and D9 rep1/rep2 sit beside them.
  - D15 and D18 rep2 are in between.
  - D18 rep1 and D22 rep2 are furthest from IgG.
  - D22 rep1, the shallow library at 6.8 M fragments, is an outlier on its own.
- **Replicates pair by day only at D9,** where both look like background.

### Group IgG test (2026-11-08)

Question: does a deeper, shared IgG per day x line raise the number of reproducible peaks?

**Method**
- **Group IgG** (step 10, job 35436454): the deeper of the two block IgGs is subsampled to the
  shallower one's depth (seed 101, fraction 0.73-0.98), then merged. Each group IgG holds
  35.9-60.3 M fragments, about 2x a single block IgG. Details: `results/group_igg/group_igg_manifest.tsv`.
- **Calls against the group IgG:**
  - Step 04 with `CONTROL=group`: relaxed calls plus the ISL1 q 0.01 summit call (job 35437434).
  - Step 12: per-sample histone calls (job 35437435).
  - Step 05 with `CONTROL=group`: reproducible sets (job 35437436).
  - Step 11: matched vs group comparison (job 35437437).
- **Settings:** all as with the matched IgG; only the control changes. Steps 04, 05 and 12 ran on
  amem with 240 G and 12 CPUs, peaking at 4-9 GB.
- **Outputs:** `results/reproducibility_groupIgG/compare/{isl1_A.tsv,sets.tsv}`.

**ISL1, matched -> group (chr1-22/X/Y)**

| Day | A both-replicate | A chance | A ratio | B overlap Nt | C IDR Nt | C IDR Np | Ref match (A, top 1000) |
|---|---|---|---|---|---|---|---|
| D9 | 157 -> 256 | 99 -> 150 | 1.6 -> 1.7 | 18,781 -> 23,973 | 7 -> 6 | 29 -> 29 | 13 -> 10% |
| D15 | 530 -> 698 | 215 -> 274 | 2.5 -> 2.5 | 24,261 -> 28,095 | 80 -> 85 | 153 -> 152 | 34 -> 30% |
| D18 | 976 -> 1,151 | 134 -> 144 | 7.3 -> 8.0 | 22,132 -> 26,855 | 480 -> 536 | 960 -> 963 | 70 -> 72% |
| D22 | 1,499 -> 1,688 | 392 -> 474 | 3.8 -> 3.6 | 23,494 -> 26,745 | 593 -> 659 | 932 -> 926 | 74 -> 75% |

Per-replicate q 0.01 calls rise 5-30% (for example, D9 rep2 goes from 18.7k to 24.4k).

**Histone overlap Nt, matched -> group.** All groups pass both ways.

| Group | H3K27ac Nt | H3K27me3 Nt |
|---|---|---|
| D9 H3CTR | 20,426 -> 25,668 (+26%) | 27,448 -> 29,322 (+7%) |
| D9 ISL1KO | 13,772 -> 16,757 (+22%) | 26,421 -> 28,192 (+7%) |
| D15 H3CTR | 71,628 -> 77,629 (+8%) | 56,584 -> 60,492 (+7%) |
| D15 ISL1KO | 50,163 -> 55,459 (+11%) | 56,040 -> 59,723 (+7%) |
| D18 H3CTR | 38,763 -> 43,942 (+13%) | 53,877 -> 58,086 (+8%) |
| D18 ISL1KO | 27,101 -> 31,773 (+17%) | 60,932 -> 66,176 (+9%) |
| D22 H3CTR | 58,221 -> 65,376 (+12%) | 63,001 -> 67,915 (+8%) |
| D22 ISL1KO | 69,592 -> 75,456 (+8%) | 64,674 -> 70,223 (+9%) |

Per-sample histone calls (step 03 vs step 12), matched -> group:
- **H3K27ac:** regions +5-16%, FRiP up by 4-15% of its value.
- **H3K27me3:** regions +2-6%, FRiP up by 1-3% of its value.
- Stats are in `results/qc/peaks_groupIgG/`.

**Reading**
- **Every call set grows with the group IgG,** by 7-28% for Nt and 13-63% for ISL1 option A.
  - A deeper control gives a less noisy local lambda, so fewer regions are lost to IgG noise.
  - The largest relative gains are in the weakest groups: D9, and H3K27ac at D9 and D18.
- **Reproducibility itself does not improve.**
  - For ISL1 A, chance overlap grows as much as the real overlap, so the ratio is unchanged
    (1.6-8.0).
  - The share matching the earlier batches' reference is unchanged.
  - ENCODE self-consistency and rescue ratios move by under 0.2, and the flags are the same.
- **Np is unchanged (±3%),** as expected: the matched pooled call already used both block IgGs
  as its control.
- **ISL1 remains limited by its signal, not by the control.** D9 is still near chance under
  every method, and C stays at tens to hundreds of peaks.

### Earlier control ISL1 libraries, own IgG, by collection day (2026-11-09)

The earlier batches' control ISL1 libraries are treated as further Ctr replicates. Batch 3 H3CTR
rep1/rep2 count as Ctr rep1/rep2, and the earlier libraries become rep3 onward, named by
collection day. The mapping is in `study_config/ctr_matched_samples.tsv`, which also keeps the GFP
electroporation day and the prep/sequencing batch. KO libraries are not included.

`sbatch --array=0-6 scripts/13_ctr_matched_calls.sb` (job 35478797, amem; 6-37 min, 8.9 GB):
- **Inputs:** BAMs are read in place from ISL1_cutNrun and batch 3.
- **Command:** the same MACS3 call as step 03: q 0.01, summits, duplicates kept, own block IgG,
  blacklist filtered.
- **Output:** `results/ctr_matched/`.
- **Check:** all 7 are identical to the October calls of the same libraries under their old names.

| Sample | Was | Electroporated | Batch (prep / seq) | Fragments | IgG | Peak rows | Regions | FRiP |
|---|---|---|---|---|---|---|---|---|
| D15_Ctr_rep3 | D12 Ctr+GFP rep1 | D12 | 1 / 1 | 22.5 M | 13.5 M | 14,019 | 12,117 | 1.7% |
| D15_Ctr_rep4 | D12 Ctr+GFP rep2 | D12 | 2 / 2 | 34.6 M | 29.2 M | 7,809 | 6,308 | 4.4% |
| D18_Ctr_rep3 | D15 Ctr+GFP rep1 | D15 | 1 / 1 | 26.3 M | 13.9 M | 32,748 | 25,869 | 5.0% |
| D18_Ctr_rep4 | D15_CtrGFP_rep2 (batch 3 run) | D15 | 2 / 3 | 38.3 M | 37.8 M | 41,834 | 36,669 | 4.6% |
| D22_Ctr_rep3 | D22 Ctr rep1 | - | 1 / 1 | 19.7 M | 10.1 M | 8,878 | 8,133 | 1.8% |
| D22_Ctr_rep4 | D22 Ctr rep2 | - | 2 / 2 | 36.6 M | 28.5 M | 31,146 | 21,226 | 6.0% |
| D22_Ctr_rep5 | D22 Ctr rep3 | - | 2 / 2 | 33.6 M | 32.7 M | 18,286 | 13,699 | 3.2% |

Control ISL1 replicates per collection day: D9 2, D15 4, D18 4 and D22 5, for 15 in total.

### Control ISL1 sets from all replicates per day (2026-11-09)

**Inputs**
- `sbatch --array=0-6 scripts/14_ctr_relaxed.sb` (job 35482355): relaxed calls for Ctr rep3
  onward, each against its own IgG (p 0.01, top 300k). This gives 141k-289k peaks per library.
- `sbatch --array=0-2 scripts/15_ctr_multirep_sets.sb` (job 35482357; 33-111 min, max 23 GB):
  D15, D18 and D22. The replicate list is `study_config/ctr_replicates.tsv`.
- **Output:** `results/ctr_multi/<day>/`, chr1-22/X/Y only.

**Methods**
- **A:** merged q 0.01 calls, keeping regions with peaks from >= k replicates.
- **B:** pooled relaxed peaks that overlap the relaxed peaks of >= k replicates by >= 50%.
- **IDR:** every replicate pair against the pooled peak list.
  - k = 2: passes in at least one pair.
  - k = 3: passes in pairs that together involve >= 3 replicates.
- **Chance (A only):** every replicate shuffled independently.
- **Reference:** the earlier batches' D15 Ctr conservative set. It was built from what is now D18
  rep3, so it is not independent of the earlier libraries.

**Results**

| Day (n) | A >=2 (chance) | A >=3 (chance) | B >=2 | B >=3 | IDR >=2 | IDR >=3 | IDR best pair |
|---|---|---|---|---|---|---|---|
| D15 (4) | 3,884 (506) | 903 (2) | 51,127 | 10,366 | 3,191 | 407 | rep3-rep4, 3,166 |
| D18 (4) | 11,455 (1,170) | 4,191 (13) | 69,377 | 21,374 | 7,861 | 3,008 | rep3-rep4, 7,480 |
| D22 (5) | 11,215 (1,550) | 4,953 (18) | 81,340 | 27,040 | 9,535 | 4,180 | rep4-rep5, 9,121 |
| D22, 4-replicate mean | 8,664 (936) | 3,262 (7) | 59,620 | 18,182 | 7,174 | 2,649 | |

IDR pairs, peaks passing:

| Day | batch 3 rep1-rep2 | batch 3 x earlier | earlier x earlier |
|---|---|---|---|
| D15 | 73 | 82-346 | 3,166 |
| D18 | 431 | 409-2,708 | 7,480 |
| D22 | 522 | 547-1,783 | 3,200-9,121 |

**Reading**
- **The earlier libraries agree with each other far better** than batch 3 rep1/rep2 do, even
  though each earlier pair spans batches 1 and 2.
  - Among the batch 3 libraries, D15 rep2 and D18 rep2 agree worst with everything.
  - D18 rep1 is the best, with 2.7k IDR peaks against each earlier library.
- **A >=3 is highly specific.** Chance is 2-18 regions, or 0.2-0.4%.
  - For the 4-replicate days it needs support across batches by construction.
  - Median width is 0.85-1.0 kb, and 79-89% of the top 1000 fall in the reference set.
- **A >=2 is about 10-14% chance overlap.** Most of its regions are supported by the earlier
  libraries only.
- **A >=3 and IDR agree.**
  - 87-93% of A >=3 regions are in IDR >=2, and 40-70% are in IDR >=3.
  - 87-98% of IDR >=3 peaks are in A >=2.
- **B is too permissive** (10k-81k). Only 11-36% of all its regions fall in the reference set,
  although its top 1000 are as good as the other methods'.
- **D15 is the weakest day** (A >=3 903). With D22 cut to 4 replicates, D18 and D22 are similar
  (A >=3 4,191 vs 3,262).

### Control ISL1 from three chosen replicates per day (2026-11-09)

**Which three:** the better batch 3 library plus the earlier libraries.
- D15 and D18: rep1 + rep3 + rep4.
- D22: rep2 + rep4 + rep5.
- Basis: agreement with the earlier libraries. For D22, rep2 + rep4 + rep5 shares 91-94% of its
  >=2 regions with rep3 + rep4 + rep5.

**Run:** `SUBSET="D15:1,3,4 D18:1,3,4 D22:2,4,5" OUT_NAME=ctr_trio sbatch --array=0-2 --export=ALL
scripts/15_ctr_multirep_sets.sb` (job 35494855, 23-29 min).
- The pooled peak list for IDR and B is built from the three chosen libraries only.
- Output: `results/ctr_trio/<day>/`.

| Day | MACS3 >=2 of 3 (chance) | MACS3 3 of 3 | IDR, any pair | IDR, all 3 | IDR pairs: batch 3 x earlier, earlier x earlier |
|---|---|---|---|---|---|
| D15 | 3,299 (127) | 642 | 3,365 | 405 | 279 / 386, 3,349 |
| D18 | 10,945 (748) | 3,812 | 7,873 | 2,988 | 2,651 / 2,811, 7,501 |
| D22 | 9,804 (539) | 2,779 | 9,577 | 1,927 | 1,850 / 1,595, 9,483 |

- **IDR "any pair" is mostly the earlier x earlier pair.** The batch 3 library adds 16-372 peaks.
- **B is too permissive:** 33k-54k regions at >=2, of which only 15-19% fall in the reference set.

**Kept sets (2026-11-09): MACS3 >=2 of 3 and IDR any pair.** `bash scripts/18_isl1_ctr_final_peaks.sh`
writes them to `results/isl1_ctr_final/{macs3_ge2,idr_any_pair}/<day>_Ctr_ISL1.bed`, on scratch
and copied under the repo on /pl/active (not tracked by git).
- **Columns:** chr, start, end, name, n_reps, reps, max -log10 p, n_batch3, n_earlier.
- **Counts:** macs3_ge2 has 3,299 / 10,945 / 9,804 peaks (D15 / D18 / D22; median 612-811 bp).
  idr_any_pair has 3,365 / 7,873 / 9,577 (median 975-1,380 bp).

### H3K27ac differential analysis, D15/D18/D22 (2026-11-09)

**Inputs** (steps 16-17)
- **Regions:** union of the 6 groups' matched-IgG overlap Nt sets, giving 179,438 regions
  (238 Mb). Of these, 178,007 pass the pre-filter (>= 10 fragments in >= 2 libraries).
- **Counts:** featureCounts on 12 libraries.
- **QC PCA:**
  - PC1 (55%) is day.
  - PC2 (35%) tracks each library's fraction of fragments in the regions (r = -0.90), so it is
    likely technical.
  - rep2 sits below rep1 on PC1 in all 6 groups, which supports the replicate term.

**Model** (`scripts/histone_deseq2.R`, job 35496947, 7 min, 19 GB)
- DESeq2 `~ rep + genotype + day + genotype:day`, median-of-ratios normalization.
- Wald contrasts use `~ rep + group`, with ashr-shrunk log2FC.
- Outputs: `results/diff/H3K27ac/deseq2/` (tables, summary, ko_patterns, time_clusters,
  report.html).
- Tables write NA as empty fields.

**Results.** Counts are regions at padj < 0.05; up = higher in the first-named group.

| Test | Tested | padj < 0.05 (up / down) | and shrunk abs(log2FC) >= 0.585 | and >= 1 |
|---|---|---|---|---|
| KO vs Ctr, D15 | 8,901 | 27 (8 / 19) | 3 | 1 |
| KO vs Ctr, D18 | 178,007 | 2 (0 / 2) | 0 | 0 |
| KO vs Ctr, D22 | 53,766 | 1,034 (602 / 432) | 659 | 40 |
| Interaction D22 vs D15 | 19,254 | 50 | 38 | 1 |
| Interaction D18 vs D15 / D22 vs D18 | 178,007 | 1 / 0 | 0 | 0 |
| Ctr, D18 vs D15 | 43,412 | 2,605 | 1,593 | 145 |
| Ctr, D22 vs D18 | 33,059 | 4,391 | 2,106 | 14 |
| Ctr, D22 vs D15 | 91,728 | 15,573 | 15,290 | 5,149 |
| KO, D22 vs D15 | 84,826 | 12,985 | 12,666 | 3,150 |

Likelihood-ratio tests at padj < 0.05: time 21,314; genotype (any day) 907; interaction 41.

**Reading**
- **KO vs Ctr differences appear only at D22,** with about 1,000 regions.
  - D15 has 27.
  - D18 has 2. D18 has the weakest H3K27ac libraries, so this is low power, not proof of no
    effect. Dropping D18 H3CTR rep2 still gives 0.
- **The D22 effect differs significantly from D15** for 50 regions. It cannot be separated from
  D18, which has no power.
- **ashr shrinkage with 2 replicates is conservative.** The shrunk |log2FC| >= 1 cut keeps only
  40 of the 1,034 D22 regions, while 734 have a raw |log2FC| >= 1.
- **A cluster of KO-lost regions** lies at chr9:101-104 Mb (GRIN3A, LINC03094, LINC00587; also
  the top genotype-LRT hits). IgG coverage there matches Ctr (ratio 0.89-1.02 in both), so it is
  not a deletion in the KO line.
- **Library-size normalization** gives median KO vs Ctr log2FC of -0.21 at D15, -0.16 at D18 and
  +0.05 at D22. This follows each library's fraction of fragments in the regions, so it is not
  clear evidence of a global change.

**Re-run with a 1.5-fold cut and two-tier labels (job 35497518).** Significant = padj < 0.05 and
shrunk abs(log2FC) >= 0.585.

KO vs Ctr, up / down:

| Day | Main | Library-size normalization |
|---|---|---|
| D15 | 0 / 3 | 2 / 25 |
| D18 | 0 / 0 | 0 / 0 |
| D22 | 374 / 285 | 438 / 202 |

KO patterns and tiers:

| Pattern | Tier | Regions |
|---|---|---|
| D22 only, up | partly_confirmed (differs from D15) | 24 |
| D22 only, up | restricted_only | 350 |
| D22 only, down | partly_confirmed | 18 |
| D22 only, down | restricted_only | 265 |
| D15 only, down | confirmed | 1 |
| D15 + D22 only, down | restricted_only | 2 |

- **No D22-only region can be fully confirmed:** the D22 vs D18 interaction has no power.
- **Ctr D18 vs D15 (1,593 regions) depends on D18 H3CTR rep2.** Without that library, only 227
  of the 1,593 remain.

**More sensitive KO vs Ctr tests (step 19, job 35498604, 4 min).**
`MARK=H3K27ac sbatch --export=ALL scripts/19_histone_diff_sensitive.sb`. Output:
`results/diff/H3K27ac/sensitive/`.

Regions at padj < 0.05, up / down in KO:

| Test | No fold-change cut | 1.5-fold lfcThreshold test |
|---|---|---|
| D15 | 8 / 19 | 0 / 1 |
| D18 | 0 / 2 | 0 / 0 |
| D22 | 602 / 432 | 0 / 5 |
| Shared across days (~ rep + day + genotype) | **972 / 1,385** | 0 / 0 |

- **500 bp windows find essentially nothing.** Of 527,217 windows, 444,458 pass the pre-filter.
  No window is significant per day; the shared test gives 3 up / 11 down. The Simes combination
  per region finds 0. At about 10 fragments per window there is too little signal, and the
  multiple-testing burden is large.
- **Pooling days in the shared test is the real gain:** 2,357 regions.
- **The lfcThreshold test is very strict with 2 replicates.**

## Local review copy (Mac)

`~/Documents/Claude/ISL1_KO_histone_cutNrun` is a clone (remote `alpine`) holding code plus
HTML reports only. Refresh with `bash scripts/fetch_reports.sh`.
