# ISL1 CUT&RUN time course (control, D15 / D18 / D22)

ISL1 binding over islet differentiation in control cells: three ISL1 CUT&RUN libraries per day, each
with its own IgG, analysed from raw fastqs in one Slurm chain. Built on the lab pipeline
[cutnrun-pipeline](https://github.com/ZhuLab-CUAnschutz/cutnrun-pipeline), vendored unchanged in
`pipeline/` (commit `1afe307`).

Repo (Alpine): `/pl/active/ZhuLab_BDC/GitHub/ISL1_project/cutnrun_ISL1_timecourse` -- part of
`Han-Z-BDC/ISL1_project`. Scratch: `/scratch/alpine/hzhu10@xsede.org/ISL1_cutnrun_timecourse`.

This folder replaces, for the paper, the exploratory analyses in `cutnrun_ISL1_KO_histone/` (batch 3,
with histone marks) and `cutnrun_ISL1_batch1-2/` (earlier ISL1 batches), which record how the libraries,
the replicate choice and the methods below were chosen. Both are kept at git tag `cutnrun-exploratory`
(`git checkout cutnrun-exploratory -- <folder>` restores one); the peak-methods report that chose
macs3_ge2 is `reports/isl1_peak_methods/report.html` in the Mac copy.

## Design (decided 2026-11-12)

- **Control ISL1 only**, D15, D18 and D22; no KO, no histone marks, no D9.
- **Three libraries per day:** one batch 3 library (rep1) and two earlier libraries (rep2-3), the
  replicates that agree best with each other (chosen in `cutnrun_ISL1_KO_histone/`).
- **Question:** what changes over time. Contrasts D18 vs D15, D22 vs D18, D22 vs D15.

| Sample | Day | Batch (model) | Delivery | Original library |
|---|---|---|---|---|
| D15_Ctr_rep1 | D15 | batch3 | 7535-0003 (2026-11) | D15_H3CTR_rep1 |
| D15_Ctr_rep2 | D15 | earlier | 7535-0001 (2026-08) | D12_Ctr_rep1 (Ctr-GFP, electroporated D12) |
| D15_Ctr_rep3 | D15 | earlier | 7535-0002 (2026-10) | D12_Ctr_rep2 (Ctr-GFP, electroporated D12) |
| D18_Ctr_rep1 | D18 | batch3 | 7535-0003 | D18_H3CTR_rep1 |
| D18_Ctr_rep2 | D18 | earlier | 7535-0001 | D15_Ctr_rep1 (Ctr-GFP, electroporated D15) |
| D18_Ctr_rep3 | D18 | earlier | 7535-0003 | D15_CtrGFP_rep2 (prep batch 2, sequenced with batch 3) |
| D22_Ctr_rep1 | D22 | batch3 | 7535-0003 | D22_H3CTR_rep2 |
| D22_Ctr_rep2 | D22 | earlier | 7535-0002 | D22_Ctr_rep2 |
| D22_Ctr_rep3 | D22 | earlier | 7535-0002 | D22_Ctr_rep3 |

Each ISL1 library has an IgG of the same name (`<sample>_IgG`). Days are collection days; the earlier
libraries' original names use the electroporation day. Sample sheets: `study_config/samples.tsv`
(pipeline format) and `study_config/sample_metadata.tsv`.

## Steps

| Step | Script | What |
|---|---|---|
| 01 | `scripts/01_filter_10x.sb` | Remove 10x 3' GEX reads that leaked into the batch 3 delivery (R1 >= 10 T at bases 29-38 or >= 24 T in 29-58; both mates). The 8 libraries of delivery 7535-0003 only. |
| 02 | `pipeline/scripts/01_align.sb` | cutadapt, bowtie2 (hg38, very-sensitive-local, 10-700 bp), MAPQ 20, properly paired; duplicates marked, kept |
| 03 | `pipeline/scripts/02_build_igg_background.sb` | each IgG deduplicated (matched IgG) |
| 04 | `scripts/04_call_peaks.sb` | MACS3 per ISL1 library vs its own IgG: BAMPE, q 0.01, `--call-summits`, duplicates kept, blacklist filtered; FRiP |
| 05 | `scripts/05_consensus.sb` | per day, regions with peaks in >= 2 of 3 libraries (**macs3_ge2**, primary) and 3 of 3; shuffled-chance estimate |
| 06 | `scripts/06_diff_sites.sb` | 401 bp sites: each macs3_ge2 region centred on its strongest summit, summit +/- 200 bp, overlaps across days resolved by the stronger summit |
| 07 | `scripts/07_diffbind.sb` | IgG greylist (GreyListChIP, p 0.999, 9 IgGs); **DiffBind (DESeq2)**, `~ batch + day` (batch 3 vs earlier), background (15 kb bins) and RLE normalization side by side; check: DESeq2 on featureCounts counts, same design; FDR < 0.05 |
| 08a-c | `scripts/08{a,b,c}_*.sb` | bigWigs (deepTools, fragments, 10 bp): ISL1 scaled to CPM at the mean library size with the DiffBind background size factors, IgG plain CPM; one track per day = mean of its 3 libraries; beds; IGV session |
| 09 | `scripts/09_report.sb` | `report.html` and copy of the keepers to `results/` on /pl/active |

Analysis jobs run on `amem` (12 CPUs, 240 GB). Tracks for figures: load `results/tracks/igv_session.xml`
(hg38); the three day tracks share one autoscale group, so heights compare across days.

## Run (Alpine)

```
cd /pl/active/ZhuLab_BDC/GitHub/ISL1_project/cutnrun_ISL1_timecourse
bash study_config/setup_links.sh
bash scripts/run_all.sh
```

## Outputs

On scratch under `results/`; steps 08c and 09 copy the keepers to `results/` in this folder on
/pl/active (not tracked by git): `report.html`, `consensus/`, `diff/` (all-site tables per method,
normalization and contrast; up/down beds; summaries), `tracks/`, `qc/`, `tenx_filter/`.
