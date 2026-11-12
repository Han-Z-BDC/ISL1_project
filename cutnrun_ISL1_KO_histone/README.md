# ISL1 KO / H3 control CUT&RUN time course

CUT&RUN for ISL1, Flag, H3K27ac and H3K27me3 (each with an IgG control) across islet
differentiation (D9, D15, D18, D22) in the H3 control line and the ISL1 KO line. Analyzed
on its own, independent of the earlier ISL1 CUT&RUN batches (repo `ISL1_cutNrun`). Built on
the general lab pipeline
[cutnrun-pipeline](https://github.com/ZhuLab-CUAnschutz/cutnrun-pipeline), vendored unchanged
in `pipeline/` (copied from commit `1afe307`).

Repo (Alpine): `/pl/active/ZhuLab_BDC/GitHub/ISL1_project/cutnrun_ISL1_KO_histone` -- part of
`Han-Z-BDC/ISL1_project` (GitHub, private).

**Current focus (2026-11-12): ISL1 CUT&RUN only.** The histone marks are set aside: the H3K27ac
tracks are poor, and H3K27me3 is of limited use without H3K27ac. Histone steps ran in the
10x-filtered re-run, but their outputs stay on scratch unreviewed. The histone sections below
(H3K27ac / H3K27me3 differential analysis, DiffBind / csaw) describe the run before the 10x filter.

**Primary ISL1 peak set (2026-11-12): `macs3_ge2`**, MACS3 consensus of three control replicates
per day (one batch 3 library plus two earlier libraries), regions with peaks in >= 2 of 3. See
"Control ISL1 from three chosen replicates per day".

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
Analysis jobs run on `amem` (`--qos=mem-normal`, 12 CPUs, 240 G). `acompile` is kept for data
downloads only. The first FastQC run (job 35274348) went to `acompile` before that rule was set.

```
cd pipeline
sbatch scripts/00_fastqc.sb
sbatch --array=0-61 --cpus-per-task=8 --mem=30G --time=06:00:00 scripts/01_align.sb
```

### 10x read filter and re-run (2026-11-10 to 2026-11-12)

Some read pairs carry 10x single-cell library structure (a poly-T run in R1 after the barcode
and UMI). Step 00f (`scripts/00f_filter_10x.sb`) drops a pair when R1 has >= 10 T at bases 29-38
or >= 24 T in bases 29-58, and writes the kept pairs to `clean_fastq/`.

- **All 62 libraries:** 2,542,696,439 pairs in, 3,380,087 removed (0.133%), 794,594 of them with
  a 10x cell barcode. Per library 0.088-0.379%; the ~0.09% floor is the expected loss of genuine
  pairs.
- **Most affected:** the batch 3 rep1 H3CTR libraries (D9 rep1 ISL1 0.38%, D18 rep1 H3K27ac
  0.37%, D9 rep1 H3K27me3 0.33%).
- Per-library table: `${SCRATCH}/tenx_filter/<sample>.tsv` (pairs_in, removed,
  removed_with_10x_barcode, pairs_out, pct_removed).

The settled analysis was re-run on the filtered reads (`scripts/rerun_after_10x_filter.sh`);
results are in `${SCRATCH}/results/`. Unless a section says otherwise, results below are from
this re-run. Exploratory steps were not re-run; their outputs from before the filter are in
`${SCRATCH}/results_with10x/`: 06-12 (SEACR, step 07 pair consensus, replicate correlation,
group IgG), 20, 22 and `ctr_multi`.

| Step | Job(s) | Note |
|---|---|---|
| 01 align, 02 IgG background | 35568628, 35568629 | |
| 03 peaks, 04 relaxed + pseudoreps | 35568630, 35568631 | |
| 05 reproducible sets | 35568632 | 11 of 20 tasks finished while Slurm was down (2026-11-11 ~13:10-14:30) and were recorded as TIMEOUT; their outputs are complete |
| 13, 14 earlier control ISL1 | 35568633, 35568634 | |
| 15 `ctr_trio` | 35586302 | after 35586301, which rebuilt `isl1_consensus/{reference.bed,genome.txt}` (step 07's helper files; now `scripts/07a_isl1_reference_files.sb` in the re-run script) |
| 18 kept ISL1 sets | 35609970 | |
| 16, 17, 19, 21 histone | 35609971-74, 35610439-42 | 19 must wait for 17 |
| 23a-c tracks | 35609979, 35609980, 35610443 | |
| 24 before/after comparison | 35610444 | |
| 25 ISL1 peak-methods report | 35611761 | |

Steps 18, 23c and 24 copy their keepers under this folder on /pl/active (`results/`, not tracked).

**Effect on ISL1 (step 24, `results/tenx_check/`).** "Old" is the analysis on the raw fastqs.
- **Libraries:** the batch 3 H3CTR rep1 ISL1 libraries carry the 10x reads: 0.25-0.38% of pairs
  removed (D9 0.38%, D22 0.31%, D15 0.29%, D18 0.25%), 18k-108k of them with a 10x barcode; their
  IgGs lost 0.17-0.26%. Every other ISL1 library is at the ~0.1% floor. 10x fragments made up
  0.10-0.27% of the old filtered BAMs, and step 03 peaks per library changed by at most 0.4%.
  The earlier-batch libraries (rep3 onward, ISL1_cutNrun) were not filtered.
- **Peak sets** (lost = old region with no new overlap; gained = the reverse):

| Set | Day | Old | New | Lost | Gained | Lost on a 10x fragment | Old regions on a 10x fragment |
|---|---|---|---|---|---|---|---|
| macs3_ge2 | D15 | 3,299 | 3,299 | 1 | 1 | 1 | 728 (22%) |
| macs3_ge2 | D18 | 10,945 | 10,947 | 19 | 21 | 19 | 2,631 (24%) |
| macs3_ge2 | D22 | 9,804 | 9,803 | 1 | 0 | 1 | 2,761 (28%) |
| idr_any_pair | D15 | 3,365 | 3,368 | 1 | 4 | 1 | 1,087 (32%) |
| idr_any_pair | D18 | 7,873 | 7,870 | 25 | 23 | 17 | 2,661 (34%) |
| idr_any_pair | D22 | 9,577 | 9,577 | 1 | 2 | 0 | 3,965 (41%) |
| D9 rep1 + rep2 (MACS3 2 of 2) | D9 | 157 | 158 | 1 | 2 | 1 | 25 (16%) |

- **All 21 lost macs3_ge2 regions overlap a 10x fragment** (against 22-28% of all old regions),
  so the filter removed the regions the 10x reads had created. D18 changes most because its
  trio includes D18 rep1, one of the contaminated libraries; still 0.2% of the set. The
  replicate choice and the primary set are unaffected. Lost regions: `results/tenx_check/lost/`.

### Alignment QC (2026-11-06, before the 10x filter)

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

Results (10x-filtered re-run, job 35568630, all 44 tasks completed). Regions are distinct
blacklist-filtered peaks, and FRiP counts all filtered fragments:

| Antibody | Mode | Regions per library | FRiP |
|---|---|---|---|
| H3K27me3 | broad | 54.8k-119.1k | 33-49% |
| H3K27ac | narrow | D15/D22: 32.6k-80.9k; D9/D18: 7.2k-41.9k | D15/D22: 2.1-8.3%; D9/D18: 0.3-3.0% |
| ISL1 | narrow | 18.4k-38.0k | 1.2-4.6% |
| Flag | narrow | 12.1k, 12.4k | 1.3% |

ISL1 checks (top 1000 peaks by p value; before the 10x filter):
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

Results (10x-filtered re-run; jobs 35568631 and 35568632, all 60 tasks produced their outputs;
see the re-run table for the 11 step 05 tasks recorded as TIMEOUT).

- **Relaxed lists:** ISL1 214k-281k peaks per library, H3K27ac 183k-269k, H3K27me3 91k-182k. No
  list reached its cap (300k / 500k), so the caps never applied.

Reproducible sets. Self = N1 / N2 ratio and rescue = Np / Nt ratio. FRiP is for rep1 / rep2 in
Nt:

| Group | Method | N1 / N2 | Nt | Np | Self | Rescue | Flag | Nt Mb | Nt FRiP |
|---|---|---|---|---|---|---|---|---|---|
| D9_H3CTR_ISL1 | overlap | 23,756 / 21,251 | 18,694 | 28,581 | 1.12 | 1.53 | pass | 12.4 | 1.7 / 1.4% |
| D9_H3CTR_ISL1 | IDR | 0 / 7 | 11 | 18 | inf | 1.64 | borderline | 0.0 | 0.7 / 0.4% |
| D15_H3CTR_ISL1 | overlap | 23,594 / 34,325 | 24,223 | 40,409 | 1.45 | 1.67 | pass | 16.9 | 2.7 / 2.6% |
| D15_H3CTR_ISL1 | IDR | 36 / 16 | 86 | 146 | 2.25 | 1.70 | borderline | 0.1 | 1.3 / 1.2% |
| D18_H3CTR_ISL1 | overlap | 27,954 / 22,732 | 22,215 | 22,073 | 1.23 | 1.01 | pass | 15.5 | 3.2 / 2.4% |
| D18_H3CTR_ISL1 | IDR | 853 / 92 | 488 | 956 | 9.27 | 1.96 | borderline | 0.6 | 1.7 / 1.1% |
| D22_H3CTR_ISL1 | overlap | 25,717 / 32,027 | 23,460 | 33,129 | 1.25 | 1.41 | pass | 17.4 | 3.2 / 3.2% |
| D22_H3CTR_ISL1 | IDR | 121 / 477 | 598 | 940 | 3.94 | 1.57 | borderline | 0.8 | 1.0 / 1.6% |
| D9_H3CTR_H3K27ac | overlap | 21,554 / 17,757 | 20,408 | 25,838 | 1.21 | 1.27 | pass | 14.4 | 1.3 / 1.2% |
| D9_ISL1KO_H3K27ac | overlap | 19,531 / 17,875 | 13,756 | 22,747 | 1.09 | 1.65 | pass | 8.7 | 0.7 / 0.7% |
| D15_H3CTR_H3K27ac | overlap | 53,340 / 56,235 | 71,628 | 78,811 | 1.05 | 1.10 | pass | 97.2 | 10.8 / 11.4% |
| D15_ISL1KO_H3K27ac | overlap | 50,571 / 29,474 | 50,082 | 81,091 | 1.72 | 1.62 | pass | 59.0 | 6.5 / 5.4% |
| D18_H3CTR_H3K27ac | overlap | 42,785 / 22,312 | 37,795 | 63,799 | 1.92 | 1.69 | pass | 36.6 | 3.9 / 3.1% |
| D18_ISL1KO_H3K27ac | overlap | 16,561 / 22,097 | 27,109 | 34,302 | 1.33 | 1.27 | pass | 20.7 | 2.2 / 2.0% |
| D22_H3CTR_H3K27ac | overlap | 55,328 / 42,010 | 58,295 | 67,445 | 1.32 | 1.16 | pass | 74.4 | 10.5 / 7.6% |
| D22_ISL1KO_H3K27ac | overlap | 45,030 / 58,462 | 69,596 | 84,064 | 1.30 | 1.21 | pass | 94.5 | 11.8 / 10.8% |
| D9_H3CTR_H3K27me3 | overlap | 32,317 / 27,289 | 27,455 | 33,999 | 1.18 | 1.24 | pass | 99.5 | 29.6 / 42.6% |
| D9_ISL1KO_H3K27me3 | overlap | 27,321 / 25,412 | 26,446 | 27,185 | 1.08 | 1.03 | pass | 88.9 | 31.1 / 42.2% |
| D15_H3CTR_H3K27me3 | overlap | 65,635 / 47,147 | 56,649 | 56,854 | 1.39 | 1.00 | pass | 169.2 | 34.1 / 45.0% |
| D15_ISL1KO_H3K27me3 | overlap | 51,892 / 51,814 | 56,116 | 62,520 | 1.00 | 1.11 | pass | 166.1 | 39.2 / 42.5% |
| D18_H3CTR_H3K27me3 | overlap | 55,690 / 49,918 | 53,955 | 55,942 | 1.12 | 1.04 | pass | 167.7 | 38.5 / 45.3% |
| D18_ISL1KO_H3K27me3 | overlap | 55,117 / 59,779 | 60,980 | 62,745 | 1.08 | 1.03 | pass | 170.5 | 35.4 / 36.1% |
| D22_H3CTR_H3K27me3 | overlap | 67,526 / 53,994 | 63,098 | 65,262 | 1.25 | 1.03 | pass | 181.6 | 39.2 / 34.2% |
| D22_ISL1KO_H3K27me3 | overlap | 60,266 / 65,191 | 64,735 | 79,610 | 1.08 | 1.23 | pass | 208.6 | 39.5 / 39.3% |

- **Histone marks:** all 16 groups pass by overlap, and Nt and Np are close. H3K27ac is strong at
  D15 and D22 and weak at D9 and D18, as in the per-sample calls.
- **ISL1, overlap:** all four groups pass, with 18.7k-24.2k Nt peaks.
- **ISL1, IDR:** few peaks pass, from 11 Nt peaks (D9) to 598 (D22). The two batch 3 replicates
  rank peaks differently, which fits the poor top-peak agreement in the per-sample checks.
- **chrM:** a few chrM peaks are in each set; the blacklist has no chrM entries. Later ISL1 sets
  (steps 15, 18) keep chr1-22/X/Y only.
- **Histone marks: Nt chosen (2026-11-07).**

### ISL1 consensus peaks (nf-core/cutandrun style, 2026-11-07)

*Before the 10x filter; not re-run (outputs in `${SCRATCH}/results_with10x/`).*

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

*Before the 10x filter; not re-run (outputs in `${SCRATCH}/results_with10x/`).*

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

*Before the 10x filter; not re-run (outputs in `${SCRATCH}/results_with10x/`).*

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

`sbatch --array=0-6 scripts/13_ctr_matched_calls.sb` (first run job 35478797; 10x-filtered re-run
job 35568633; 6-37 min, 8.9 GB):
- **Inputs:** BAMs are read in place from ISL1_cutNrun and batch 3.
- **Command:** the same MACS3 call as step 03: q 0.01, summits, duplicates kept, own block IgG,
  blacklist filtered.
- **Output:** `results/ctr_matched/`.
- **Check:** the six batch 1/2 libraries are identical to the October calls under their old names.
  D18 rep4 is a batch 3 run, so its re-run uses the 10x-filtered reads.

| Sample | Was | Electroporated | Batch (prep / seq) | Fragments | IgG | Peak rows | Regions | FRiP |
|---|---|---|---|---|---|---|---|---|
| D15_Ctr_rep3 | D12 Ctr+GFP rep1 | D12 | 1 / 1 | 22.5 M | 13.5 M | 14,019 | 12,117 | 1.7% |
| D15_Ctr_rep4 | D12 Ctr+GFP rep2 | D12 | 2 / 2 | 34.6 M | 29.2 M | 7,809 | 6,308 | 4.4% |
| D18_Ctr_rep3 | D15 Ctr+GFP rep1 | D15 | 1 / 1 | 26.3 M | 13.9 M | 32,748 | 25,869 | 5.0% |
| D18_Ctr_rep4 | D15_CtrGFP_rep2 (batch 3 run) | D15 | 2 / 3 | 38.3 M | 37.8 M | 41,779 | 36,624 | 4.6% |
| D22_Ctr_rep3 | D22 Ctr rep1 | - | 1 / 1 | 19.7 M | 10.1 M | 8,878 | 8,133 | 1.8% |
| D22_Ctr_rep4 | D22 Ctr rep2 | - | 2 / 2 | 36.6 M | 28.5 M | 31,146 | 21,226 | 6.0% |
| D22_Ctr_rep5 | D22 Ctr rep3 | - | 2 / 2 | 33.6 M | 32.7 M | 18,286 | 13,699 | 3.2% |

Control ISL1 replicates per collection day: D9 2, D15 4, D18 4 and D22 5, for 15 in total.

### Control ISL1 sets from all replicates per day (2026-11-09)

*Before the 10x filter; not re-run (outputs in `${SCRATCH}/results_with10x/`).*

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
- Basis: agreement with the earlier libraries (chosen on the data before the 10x filter). For
  D22, rep2 + rep4 + rep5 shared 91-94% of its >=2 regions with rep3 + rep4 + rep5.

**Run:** `SUBSET="D15:1,3,4 D18:1,3,4 D22:2,4,5" OUT_NAME=ctr_trio sbatch --array=0-2 --export=ALL
scripts/15_ctr_multirep_sets.sb` (10x-filtered re-run: job 35586302).
- The pooled peak list for IDR and B is built from the three chosen libraries only.
- Output: `results/ctr_trio/<day>/`, chr1-22/X/Y only.

Regions per set. Chance: every replicate's regions shuffled independently and counted alike (A
only). Reference: share of the top 1000 / all regions in the earlier batches' D15 Ctr conservative
set (QC only; built from what is now D18 rep3).

| Day | Method | k | Regions | Chance | Mb | Median bp | Batch 3 + earlier | Earlier only | Reference |
|---|---|---|---|---|---|---|---|---|---|
| D15 | A: MACS3 | >= 2 | 3,299 | 126 | 2.22 | 611 | 904 | 2,395 | 85 / 74% |
| D15 | A: MACS3 | 3 | 639 | 0 | 0.59 | 857 | 639 | 0 | 82 / 82% |
| D15 | B: overlap | >= 2 | 33,546 | | 21.08 | 542 | 22,163 | 11,383 | 85 / 19% |
| D15 | B: overlap | 3 | 5,846 | | 5.92 | 920 | 5,846 | 0 | 85 / 49% |
| D15 | IDR | >= 2 (any pair) | 3,368 | | 3.59 | 974 | 424 | 2,944 | 85 / 77% |
| D15 | IDR | 3 | 408 | | 0.58 | 1,349 | 408 | 0 | 84 / 84% |
| D18 | A: MACS3 | >= 2 | 10,947 | 723 | 8.12 | 675 | 5,136 | 5,811 | 89 / 60% |
| D18 | A: MACS3 | 3 | 3,813 | 3 | 3.66 | 899 | 3,813 | 0 | 89 / 76% |
| D18 | B: overlap | >= 2 | 50,058 | | 33.68 | 565 | 35,273 | 14,785 | 89 / 17% |
| D18 | B: overlap | 3 | 15,496 | | 15.77 | 924 | 15,496 | 0 | 89 / 43% |
| D18 | IDR | >= 2 (any pair) | 7,870 | | 9.16 | 1,076 | 3,358 | 4,512 | 89 / 76% |
| D18 | IDR | 3 | 2,995 | | 4.13 | 1,301 | 2,995 | 0 | 89 / 84% |
| D22 | A: MACS3 | >= 2 | 9,803 | 526 | 8.69 | 811 | 3,713 | 6,090 | 86 / 58% |
| D22 | A: MACS3 | 3 | 2,779 | 4 | 3.26 | 1,083 | 2,779 | 0 | 86 / 76% |
| D22 | B: overlap | >= 2 | 54,302 | | 48.65 | 751 | 35,972 | 18,330 | 87 / 15% |
| D22 | B: overlap | 3 | 13,768 | | 19.51 | 1,296 | 13,768 | 0 | 87 / 40% |
| D22 | IDR | >= 2 (any pair) | 9,577 | | 14.32 | 1,380 | 2,016 | 7,561 | 87 / 59% |
| D22 | IDR | 3 | 1,925 | | 3.60 | 1,721 | 1,925 | 0 | 86 / 81% |

IDR peaks passing per pair: D15 rep1-rep3 283, rep1-rep4 385, rep3-rep4 3,352; D18 rep1-rep3
2,661, rep1-rep4 2,822, rep3-rep4 7,500; D22 rep2-rep4 1,848, rep2-rep5 1,598, rep4-rep5 9,482.

- **No region is supported by the batch 3 library alone** in any set.
- **IDR "any pair" rests mostly on the earlier x earlier pair:** 87% (D15), 57% (D18) and 79%
  (D22) of its peaks come from that pair only.
- **B is too permissive:** 33k-54k regions at >= 2, of which only 15-19% fall in the reference set.
- **A 3 of 3 is near chance-free** (0-4 shuffled regions).

**Kept sets:** `bash scripts/18_isl1_ctr_final_peaks.sh` (job 35609970) writes MACS3 >= 2 of 3 and
IDR any pair to `results/isl1_ctr_final/{macs3_ge2,idr_any_pair}/<day>_Ctr_ISL1.bed`, on scratch
and copied under this folder on /pl/active (not tracked by git).
- **Columns:** chr, start, end, name, n_reps, reps, max -log10 p, n_batch3, n_earlier.

| Day | macs3_ge2 | Median bp | idr_any_pair | Median bp | macs3_ge2 inside IDR | IDR inside macs3_ge2 |
|---|---|---|---|---|---|---|
| D15 | 3,299 | 611 | 3,368 | 974 | 76% | 73% |
| D18 | 10,947 | 675 | 7,870 | 1,076 | 67% | 91% |
| D22 | 9,803 | 811 | 9,577 | 1,380 | 82% | 82% |

**Primary ISL1 set (2026-11-12): `macs3_ge2`.** The nf-core/cutandrun replicate-threshold rule;
it uses the batch 3 library more (27-47% of regions are batch 3 + earlier, against 13-43% for
IDR), has a measured chance rate (3.8-6.6%) and narrower regions. `idr_any_pair` stays on disk.
For a high-confidence subset, MACS3 3 of 3 (`results/ctr_trio/<day>/A/k3.bed`).

**Report:** `scripts/25_isl1_peak_methods_report.sb` (job 35611761) writes all of the above, with
per-library calls and the batch 3 rep1/rep2 sets, to `results/qc/isl1_peak_methods/report.html`.

### H3K27ac differential analysis, D15/D18/D22 (2026-11-09)

*Histone marks set aside (2026-11-12). Results below are from before the 10x filter; the re-run outputs (steps 16, 17, 19, 21) are in `${SCRATCH}/results/diff/`, not reviewed.*

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

### H3K27me3 differential analysis, D15/D18/D22 (2026-11-09)

*Histone marks set aside (2026-11-12). Results below are from before the 10x filter; the re-run outputs (steps 16, 17, 19, 21) are in `${SCRATCH}/results/diff/`, not reviewed.*

**Inputs**
- **Domains:** union of the 6 groups' matched-IgG Nt sets, merged with a 3 kb gap (step 16,
  `MERGE_GAP=3000`). This gives 51,412 domains (369 Mb); 51,352 pass the pre-filter.
- **QC PCA:** PC1 (81%) is day. **PC2 (11%) is genotype:** all KO libraries sit at +12 to +14,
  Ctr at D18/D22 at -13 to -24.

**Original method**
- Step 17 (whole domains, job 35499125): 1.5-fold shrunk cut. Drop-one check uses D22 H3CTR rep2,
  the lowest-FRiP library.
- Step 19 (`WIN=5000`, job 35499126): padj only, shared test, and 5 kb windows.

**csaw-style method** (step 20, `scripts/histone_csaw_style.R`, job 35499127, 11 min)
- 2 kb and 5 kb windows over the domains, kept at >= 2x the 10 kb-bin background.
- Three normalizations: background bins, library size, median-of-ratios.
- Simes per domain; sub-domains are joined 2 kb windows, capped at 30 kb.
- Size factors from the three normalizations agree within about 10%, so there is no sign of a
  global shift.

**Results.** Domains significant at padj < 0.05, total (up / down; up = later day or higher in KO).

| Test | Original, whole domain | Original, 1.5x shrunk | csaw-style, background | csaw-style, libsize | csaw-style, median-of-ratios |
|---|---|---|---|---|---|
| Ctr D18 vs D15 | 2,235 (1,303 / 932) | 1,022 | 1,300 | 2,027 | 1,782 |
| Ctr D22 vs D18 | 9,739 (3,595 / 6,144) | 6,336 | 6,855 | 9,591 | 8,970 |
| Ctr D22 vs D15 | 14,434 (6,122 / 8,312) | 11,436 | 13,375 | 15,947 | 14,788 |
| KO vs Ctr D15 | 357 (230 / 127) | 179 | 263 (211 / 52) | 420 | 278 |
| KO vs Ctr D18 | 3,991 (1,236 / 2,755) | 1,591 | 2,499 (427 / 2,072) | 3,135 | 3,110 |
| KO vs Ctr D22 | 4,810 (2,114 / 2,696) | 3,218 | 3,459 (1,327 / 2,132) | 4,522 | 4,391 |
| KO vs Ctr shared | 6,753 (2,793 / 3,960) | - | 4,780 | 6,366 | 6,697 |

**Day-specific KO patterns** (1.5x cut, with tiers):
- D18+D22 lower in KO: 117 confirmed, 158 partly confirmed.
- D22 only, lower / higher in KO: 64 / 46 confirmed.
- D18 only, lower in KO: 21 confirmed.
- All days, consistent: 45 lower, 69 higher.
- Interaction LRT: 1,455.

**Reading**
- **The KO effect on H3K27me3 is much larger than on H3K27ac.** It is small at D15, appears at D18
  (mostly loss in KO) and persists at D22.
- **The methods agree.**
  - 73-89% of the csaw-style background calls are also original whole-domain calls.
  - Windows add few domains: 74-1,113 per test, versus 125-1,331 found only at whole-domain level.
  - Most of the difference between methods comes from normalization (background gives the fewest
    calls), not from windows vs domains.
- **D22 KO vs Ctr depends partly on D22 H3CTR rep2:** without it, 2,010 of 3,218 calls (1.5x cut)
  remain.

### Cross-checks with DiffBind and csaw (2026-11-09)

*Histone marks set aside (2026-11-12). Results below are from before the 10x filter; the re-run outputs (steps 16, 17, 19, 21) are in `${SCRATCH}/results/diff/`, not reviewed.*

**Setup**
- Env: `scripts/00e_build_diffbind_env.sb` (job 35499514) builds DiffBind 3.12.0, csaw 1.36.0,
  edgeR 4.0.16.
- DiffBind: `scripts/21_histone_diffbind.sb`. For H3K27me3, use broad per-library calls and no
  summit re-centring.
- csaw: `scripts/22_histone_csaw.sb`.
- Everything uses the same 12 libraries. Tests are at FDR / padj < 0.05, with no fold-change cut.
- DiffBind counting with 12 cores ran out of 240 GB; it now uses 4.

**H3K27ac, KO vs Ctr**

| Method | D15 | D18 | D22 | Shared |
|---|---|---|---|---|
| Ours (DESeq2, full-width union regions) | 27 | 2 | 1,034 | 2,357 |
| DiffBind on our union regions, full width, median-of-ratios | 42 | 20 | 1,381 | 3,234 |
| DiffBind own consensus, 401 bp summit sites (rip / libsize / background) | 1 / 1 / 1 | 0 / 4 / 9 | 0 / 0 / 0 | 60 / 240 / 293 |
| csaw, 150 bp + 1 kb genome windows (composition / efficiency) | 4 / 22 | 2 / 7 | 27 / 302 | 380 / 402 |

**H3K27ac, Ctr D22 vs D15**

| Method | Regions |
|---|---|
| Ours | 15,573 |
| DiffBind on our union | 17,232 |
| DiffBind consensus (rip) | 3,400 |
| csaw (efficiency) | 3,209 |

**H3K27ac reading**
- **DiffBind on our union regions reproduces our calls.** It contains 1,013 of our 1,034 D22 regions
  and all 27 at D15.
- **DiffBind's summit sites cover the same places but miss the effect.** They overlap 938 of our 1,034
  D22 regions, all in the same direction. The 401 bp summit windows measure much smaller fold changes
  (e.g. log2FC 1.74 in the region vs 0.64 at the summit). The D22 KO change lies across the regions,
  not at the summits.
- **csaw tests only 9,533 regions** (windows at >= 3x background), the strongest sites. With efficiency
  normalization it finds 302 D22 regions, 153 of them inside our calls.
- **Library-size, composition and background normalization make D18 time changes one-directional.**
  For example, Ctr D18 vs D15 under DiffBind libsize gives 6 up / 2,227 down. This follows the weak
  D18 libraries.

**H3K27me3, KO vs Ctr**

| Method | D15 | D18 | D22 | Shared |
|---|---|---|---|---|
| Ours (whole domains) | 357 | 3,991 | 4,810 | 6,753 |
| DiffBind on our domains, median-of-ratios | 802 | 4,690 | 4,936 | 6,609 |
| DiffBind own broad consensus, median-of-ratios | 842 | 5,213 | 5,488 | 7,045 |
| DiffBind own broad consensus, background | 487 | 3,061 | 3,462 | 4,340 |

**H3K27me3 reading:** every method and region set agrees. About 75-85% of DiffBind's calls fall in
our calls; csaw-style background calls are close to DiffBind background.

**H3K27me3 with csaw** (job 35501130, 85 min)
- Windows: genome-wide 2 kb and 500 bp, kept at >= 2x background (849k and 2.9M windows), merged
  into 85,356 regions of up to 30 kb.
- Normalization factors are 0.89-1.11 under both methods, so again there is no global shift.

| Normalization | D15 | D18 | D22 | Shared |
|---|---|---|---|---|
| Composition | 325 (256 up / 69 down) | 3,542 (649 / 2,892) | 4,797 (1,980 / 2,817) | 6,095 |
| Efficiency | 659 | 6,248 | 7,811 | 6,734 |

Ctr D22 vs D15: 18,876 (composition) / 23,900 (efficiency).

About 70-80% of csaw's regions with composition normalization fall in our calls. As with every
other method, the KO effect on H3K27me3 starts at D18 and is mostly loss in KO.

## Local review copy (Mac)

`~/Documents/git_repository/ISL1_project/cutnrun_ISL1_KO_histone` (clone of `Han-Z-BDC/ISL1_project`)
holds code plus HTML reports in `reports/` and downloaded tracks in `results/tracks/` (both
gitignored).
