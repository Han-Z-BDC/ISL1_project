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
`acpu` runs the alignment array. `acompile` (8 cores per user, 12 h) suits single small jobs such
as FastQC.

```
cd pipeline
sbatch --partition=acompile --qos=compile scripts/00_fastqc.sb
sbatch --array=0-61 --cpus-per-task=8 --mem=30G --time=06:00:00 scripts/01_align.sb
```

Peak calling (planned, ENCODE-style):
- **ISL1 and Flag:** MACS3 narrow peaks vs IgG.
- **H3K27ac:** MACS3 narrow peaks.
- **H3K27me3:** MACS3 `--broad`.
- **Histone reproducibility:** replicated peaks by overlap (ENCODE's histone rule) rather than IDR.

## Local review copy (Mac)

`~/Documents/Claude/ISL1_KO_histone_cutNrun` is a clone (remote `alpine`) holding code plus
HTML reports only. Refresh with `bash scripts/fetch_reports.sh`.
