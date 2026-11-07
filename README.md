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

## Local review copy (Mac)

`~/Documents/Claude/ISL1_KO_histone_cutNrun` is a clone (remote `alpine`) holding code plus
HTML reports only. Refresh with `bash scripts/fetch_reports.sh`.
