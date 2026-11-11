# Handoff: batch 3 CUT&RUN re-run after the 10x read filter (2026-11-11)

Written by the Claude session that submitted and watched the re-run, for the session that takes
over in `ISL1_project`. Read all of it before running anything. The older handoff (2026-11-10)
covers the background; this note replaces it for everything about the re-run.

## Standing rules from the user (follow exactly)
- **Do not push to GitHub** unless the user says so in your session.
- **Download only scripts and HTML reports** to the Mac, except the browser tracks, which the user
  asked for.
- **No jobs on `acompile`** (downloads only). Analysis jobs run on **amem**:
  `-p amem --qos=mem-normal --mem=240G --cpus-per-task=12`. Do not lower time walls. The user saw the
  peak memory use per step (most need < 30 G, DiffBind 116 G with `DB_CORES=4`) and chose to keep this.
- **Discuss with the user before doing anything new.** Propose first, act after a yes.
- **Do not edit, commit in, move or delete `/pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun`**
  while the chain below is queued: its jobs source that repo's `study_config/config.sh` at start.
  Do not cancel or resubmit any chain job unless one fails.

## Where things are
- This folder: Alpine `/pl/active/ZhuLab_BDC/GitHub/ISL1_project/cutnrun_ISL1_KO_histone/`,
  Mac `~/Documents/git_repository/ISL1_project/cutnrun_ISL1_KO_histone/`. `README.md` here is
  the lab notebook (every step, command, job id, result table).
- Old repo, still used by the running jobs: `/pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun`
  (at 8a92179, same code as here before the path change 6f3e7b1).
- Scratch (unchanged): `S=/scratch/alpine/hzhu10@xsede.org/ISL1_KO_histone_cutNrun`
  - `results/` the re-run on the filtered fastqs; `results_with10x/` the old analysis
  - `clean_fastq/`, `tenx_filter/<sample>.tsv` (step 00f), `logs/` (one file per job:
    `<jobname>.<jobid>[_<task>].out`; pipeline steps 01/02 log under `logs/pipeline/`)
- SSH: `ssh -i ~/.ssh/id_rsa hzhu10@xsede.org@login.rc.colorado.edu`, filter output with
  `grep -v "consoletype\|tty (1)"`. Wrap `squeue`/`sacct` in `timeout 60`: the Slurm controller
  sometimes times out (`Socket timed out`); just retry.

## The 10x filter (step 00f, done)
Drop a pair if R1 has >= 10 T at bases 29-38 or >= 24 T in bases 29-58. All 62 libraries:
2,542,696,439 pairs in, 3,380,087 removed (0.133%), 794,594 of them with a 10x cell barcode.
Per library 0.088-0.379%. The ~0.09% floor is the expected loss of genuine pairs; the batch 3
rep1 H3CTR libraries are worst (D9 rep1 ISL1 0.38%, D18 rep1 H3K27ac 0.37%, D9 rep1 H3K27me3 0.33%).
Per-library table: `$S/tenx_filter/*.tsv` (sample, pairs_in, removed, removed_with_10x_barcode,
pairs_out, pct_removed). Not yet in the README.

## The job chain (submitted 2026-11-10 by `scripts/rerun_after_10x_filter.sh`)
All on amem, 12 CPUs, 240 G. Old tracks job 35562672 was cancelled (user's choice) and the old
results moved to `results_with10x/` before submission.

| Step | Job id | Name | What it does | State at 2026-11-11 ~11:00 |
|---|---|---|---|---|
| 01 | 35568628 | cutnrun_align | align 62 libraries (bowtie2), filtered BAMs | done, 62/62, max 2 h 01 m |
| 02 | 35568629 | cutnrun_igg_bg | dedup matched IgG backgrounds | done, 58 m |
| 03 | 35568630 | b3_cnr_peaks | MACS3 peaks, 44 tasks | 40 done, 4 running |
| 04 | 35568631 | b3_relaxed | relaxed calls + pseudoreps, 40 tasks | running |
| 05 | 35568632 | b3_repro | reproducible sets (IDR / overlap) | pending |
| 13 | 35568633 | ctr_matched | earlier-batch control ISL1 matched calls | running |
| 14 | 35568634 | ctr_relaxed | earlier-batch control ISL1 relaxed calls | pending |
| 15 | 35568635 | ctr_multirep | three-replicate sets, `ctr_trio` (D15 1,3,4; D18 1,3,4; D22 2,4,5) | pending |
| 18 | 35568636 | isl1_final | kept ISL1 sets (MACS3 >= 2 of 3, IDR any pair) | pending |
| 16 | 35568637 / 38 | hist_counts | H3K27ac union / H3K27me3 3 kb domains, counts | pending |
| 17 | 35568639 / 40 | hist_deseq2 | DESeq2 (drop-one D18 H3CTR rep2 ac / D22 H3CTR rep2 me3) | pending |
| 19 | 35568641 / 42 | hist_sens | sensitive windows (500 bp ac / 5 kb me3) | pending |
| 21 | 35568643 / 44 | hist_diffbind | DiffBind (me3: macs3_broad peaks, broad calls) | pending |
| 23a | 35568645 | tracks_prep | DESeq2 size factors on analysis regions | pending |
| 23b | 35568646 | tracks_bw | 70 bigWigs, 20 at a time | pending |
| 23c | 35568647 | tracks_fin | replicate means, beds, IGV XMLs, copy to /pl/active | pending |
| 24 | 35568648 | tenx_cmp | before/after comparison vs `results_with10x` | pending |

Pending jobs held with `QOSMaxCpuPerUserLimit` are normal: the QOS caps a user at 256 CPUs (21 jobs).
Not re-run (exploratory; old outputs stay in `results_with10x/`): 06-12, 20, 22, ctr_multi.

## How it was watched
A background loop on the Mac polled every 15 min:
`sacct -n -X -P -o JobID,JobName,State,Elapsed -j $(seq -s, 35568629 35568648)` plus `squeue`,
stopping on any FAILED / TIMEOUT / OUT_OF_MEMORY / CANCELLED / NODE_FAIL or a
`DependencyNeverSatisfied` reason, or when 35568647 is COMPLETED. That watcher dies with the old
session; start your own the same way (background Bash, not short sleeps).
If a step fails: read `$S/logs/<name>.<jobid>*.out`, fix, and resubmit that step and everything
downstream of it (dependencies are in `scripts/rerun_after_10x_filter.sh`). **Resubmit from the
ISL1_project copy** (its scripts source the new config path; same scratch outputs). Tell the user first.

## Where outputs land (important)
Steps 18, 23c and 24 copy their keepers to `${STUDY_DIR}/results/...`, and `STUDY_DIR` comes from the
config the job sourced, i.e. the **old repo**:
- `/pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun/results/isl1_ctr_final/` (18)
- `/pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun/results/tracks/` (23c)
- `/pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun/results/tenx_check/` (24)
The `results/` copied into ISL1_project on 2026-11-11 holds the **pre-filter** versions. After the
chain finishes, copy these three folders (and anything else the user wants) into
`ISL1_project/cutnrun_ISL1_KO_histone/results/`, replacing the old ones, after asking the user.
All scratch outputs are also in `$S/results/`.

## What to check when each part finishes
- **18**: `results/isl1_ctr_final/` sets per day. Compare sizes with the old ones:
  macs3_ge2 3,299 / 10,945 / 9,804; idr_any_pair 3,365 / 7,873 / 9,577 (D15 / D18 / D22).
- **21**: DiffBind counts vs old: H3K27ac KO vs Ctr 42 / 20 / 1,381 / shared 3,234;
  H3K27me3 802 / 4,690 / 4,936 / 6,609.
- **23c**: 70 bigWigs in `results/tracks/replicates/`, group means in `groups/`, beds in `beds/`,
  three `igv_session_*.xml`, and the copy on /pl/active. Then download the tracks (the user asked):
  ```
  rsync -a --exclude 'sf/*.counts.txt' -e "ssh -i ~/.ssh/id_rsa" \
    hzhu10@xsede.org@login.rc.colorado.edu:/pl/active/ZhuLab_BDC/GitHub/ISL1_KO_histone_cutNrun/results/tracks/ \
    ~/Documents/git_repository/ISL1_project/cutnrun_ISL1_KO_histone/results/tracks/
  ```
  (or from the ISL1_project copy once moved there). `results/` is gitignored.
- **24**: `results/tenx_check/{libraries,sets,diff}.tsv`, `lost/*.bed`. Report to the user: 10x
  fragments that were in each old BAM, peaks and DiffBind regions lost/gained, and whether lost
  peaks sat on 10x fragments.

## After the chain
1. Report steps 18 / 21 / 24 to the user as above; download the tracks.
2. Add a README section here: the 10x filter (rule, per-library table from `tenx_filter/`), re-run job
   ids (table above), tracks, before/after results. Commit in ISL1_project with
   `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Push only if the user says so.
3. Once nothing reads the old repo any more, the user's planned cleanup of
   `ISL1_KO_histone_cutNrun` can go ahead (ask the user; the other session is handling it).

## Mac-side report builds
None were built by this session. Earlier HTML reports live under `reports/` (copied into this
folder). Build or download any new ones into
`~/Documents/git_repository/ISL1_project/cutnrun_ISL1_KO_histone/`, not into
`~/Documents/Claude/ISL1_KO_histone_cutNrun` (old Mac clone, retired).

## Open decisions waiting on the user (do not act without asking)
- Normalization and fold-change cut for the final differential tables.
- Porting the tiered day-specific labels to DiffBind.
- A final annotated output step.
- Copying the histone peak sets to /pl/active.
- Integration analyses; the spike-in question; an earlier unanswered scRNA question.
- Optional offer still open: compare ISL1 peaks with 3' ends of highly expressed genes from the
  scRNA data, as a check for residual 10x signal.
