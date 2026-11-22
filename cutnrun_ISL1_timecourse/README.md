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
| 07 | `scripts/07_diffbind.sb` | IgG greylist (GreyListChIP, p 0.999, 9 IgGs); **DiffBind (DESeq2)**, `~ batch + day` (batch 3 vs earlier), background (15 kb bins) and RLE normalization side by side; check: DESeq2 on featureCounts counts, same design; FDR < 0.05; also writes `counts_diffbind.tsv` (DiffBind reads per site and library), which `integration/` uses with the DiffBind background size factors |
| 08a-c | `scripts/08{a,b,c}_*.sb` | bigWigs (deepTools, fragments, 10 bp): ISL1 scaled to CPM at the mean library size with the DiffBind background size factors, IgG plain CPM; one track per day = mean of its 3 libraries; beds; IGV session |
| 10a-c | `scripts/10{a,b,c}_*.sb` | HOMER motifs (v5.1) on each day's macs3_ge2 consensus peaks, 200 bp windows centred on the strongest summit, known + de novo (8/10/12 bp), repeats masked (hg38 soft-masked with RepeatMasker); two backgrounds per day: genome (GC-matched random windows) and the other two days' peaks minus those overlapping the target day; 10c: genomic annotation (annotatePeaks.pl) of each day's peaks and of 50,000 random genomic windows |
| 11a-b | `scripts/11{a,b}_*.sb` | HOMER with backgrounds matched to the peaks' genomic composition (about half the peaks are within 1 kb of a TSS): A, promoter peaks vs windows near random TSSs and distal peaks vs genome windows > 1 kb from a TSS; B, all peaks vs 3 x random windows of each genomic class (promoter, intron, intergenic, exon, TTS) |
| 12 | `scripts/12_peak_overlap.sb` (acpu, 2 CPUs, 7 GB) | overlap of the per-day peak calls on the step 07 sites: each site's days in the macs3_ge2 consensus (7 overlap classes), each library's own MACS3 calls of the sites, and peak presence checked against the DiffBind counts (background normalization); `peak_overlap/` |
| 13 | `scripts/13_abc_targets.sb` (acpu, 4 CPUs, 14 GB; `scripts/abc_targets.R`) | ISL1 target genes (user, 2026-11-22; replaces the GREAT / bulk links of `integration/` for the paper): the D15 / D18 / D22 macs3_ge2 regions pooled and merged (14,300 sites, days kept); a site links to a gene by **TSS** (within +/- 1 kb of any transcript start, 10x GRCh38 2020-A GTF as the scRNA-seq) or **ABC** (overlaps an ABC element with ABC >= 0.02, own-promoter links included, in ENP2, ENP3 or early_SC_alpha of the diff-multiome ABC run, `/pl/active/ZhuLab_BDC/Han_analysis/diff_multiome_2026/out/abc/results`; PP2 left out); several genes per site allowed; `abc_targets/`: `links.tsv`, `genes.tsv`, `sites.tsv`, `summary.tsv`, `unmatched.tsv` |
| 09 | `scripts/09_report.sb` (acpu, 2 CPUs, 7 GB) | `report.html` and copy of the keepers to `results/` on /pl/active |

Analysis jobs run on `amem` (12 CPUs, 240 GB), except steps 12 and 09 (small, acpu). Tracks for figures: load `results/tracks/igv_session.xml`
(hg38); the three day tracks share one autoscale group, so heights compare across days.

## Peak-call overlap (step 12)

Each of the 14,431 sites carries the days whose macs3_ge2 consensus peaks contain it; each site is also overlapped with
every library's own peaks. Peak presence is a threshold call without a test, so it is compared with the DiffBind counts
(background normalization).

| Called on | Sites | Median log2FC D18 vs D15 / D22 vs D15 / D22 vs D18 |
|---|---|---|
| all three days | 2,450 | 0.38 / 0.24 / -0.07 |
| D18 + D22 | 4,130 | 0.41 / 0.35 / -0.01 |
| D18 only | 3,958 | 0.41 / 0.12 / -0.14 |
| D22 only | 3,153 | 0.12 / 0.34 / 0.13 |
| D15 + D18 | 368 | 0.17 / -0.20 / -0.22 |
| D15 only | 230 | -0.19 / -0.28 / -0.05 |
| D15 + D22 | 142 | 0.02 / 0.17 / 0.08 |

- **The day calls mostly follow the earlier libraries.** The batch 3 libraries call 6.5-34% of the sites (the earlier
  libraries 28-82%) although they have as many peaks of their own; the earlier D15 libraries are the weakest (7,809 and
  14,019 peaks), so part of the small D15 set is library strength.
- **Presence agrees with the direction of the counts, rarely with significance.** Sites called on D18 but not D15: 95%
  have a positive log2FC, 7.6% FDR < 0.05; sites called on both days rise almost as much (median 0.36 vs 0.41), so most
  "gained" calls are the overall rise after D15 crossing the threshold. D22 vs D18: direction agrees (85% of D22-not-D18
  sites up, 12% of D18-not-D22 sites), no site significant; about half of the sites missing on D18 or D22 are called by
  1 of that day's 3 libraries.
- The overlap classes describe the same pattern as the counts; they are not independent evidence and carry the batch
  caveat. Their target genes against the bulk KO effect are in `integration/` (step 09).

## ISL1 target genes (step 13)

TSS or ABC link (see Steps). 10,978 of 14,300 sites link to >= 1 gene (TSS 7,440; ABC 10,694; ABC only 3,538); 17,908
genes, 14,753 protein coding (76% of the GTF's 19,352): 5,068 by TSS only, 4,269 by TSS and enhancer, 8,571 by ABC
enhancer only. TSS links alone: 9,337 genes (7,354 protein coding, 38%). Most ABC-only links are weak (median ABC
0.034); ABC >= 0.1 with TSS: 10,847 genes. PAX6 is TSS-linked here (the GREAT domains of `integration/` gave its
intronic promoter site to the lncRNA PAUPAR). The 0.02 rule is broad; a stricter ABC cut is a choice still open.

## Run (Alpine)

```
cd /pl/active/ZhuLab_BDC/GitHub/ISL1_project/cutnrun_ISL1_timecourse
bash study_config/setup_links.sh
bash scripts/run_all.sh
```

## Outputs

On scratch under `results/`; steps 08c and 09 copy the keepers to `results/` in this folder on
/pl/active (not tracked by git): `report.html`, `consensus/`, `diff/` (all-site tables per method,
normalization and contrast; up/down beds; summaries), `homer/`, `peak_overlap/`, `tracks/`, `qc/`, `tenx_filter/`.
