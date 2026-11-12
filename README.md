# ISL1_project

All analyses for the ISL1 stem-cell islet differentiation project, one folder per assay. Each folder keeps the
history of the repo it came from and copies in the lab pipeline it is built on (`pipeline/`, at a recorded commit).
Data and large outputs stay on Alpine storage (scratch / `/projects`); this repo holds code, configs and small tables.

| Folder | Assay | Built on | Came from |
|---|---|---|---|
| `bulkRNAseq_ISL1_KO_timecourse/` | ISL1 KO time-course bulk RNA-seq, run YJZ5LS (Ctr vs ISL1 KO, D15 / D18 / D22; D9 not used) | own `pipeline/` (bulk-rnaseq-pipeline to be built) | new |
| `cutnrun_ISL1_timecourse/` | ISL1 CUT&RUN time course for the paper: control, D15 / D18 / D22, 3 libraries per day (from batches 1-3), peaks, consensus, differential binding, tracks | `cutnrun-pipeline` | new (2026-11-12) |
| `cutnrun_ISL1_KO_histone/` | CUT&RUN for ISL1, Flag, H3K27ac, H3K27me3 (+ IgG), D9-D22, H3 control vs ISL1 KO | `cutnrun-pipeline` | `ISL1_KO_histone_cutNrun` |
| `cutnrun_ISL1_batch1-2/` | ISL1 CUT&RUN (+ IgG), D12 / D15 / D22, control vs ISL1 KO (earlier batches) | `cutnrun-pipeline` | `ISL1_cutNrun` |
| `scRNAseq/` | scRNA-seq, 10x GEM-X On-Chip Multiplexing | `single-cell-pipeline` | `ISL1_scRNAseq` |
| `integration/` | ISL1 CUT&RUN x bulk RNA-seq: direct ISL1 targets and their time-dependent KO effect | | new |

Each assay folder has its own README with its design, steps and outputs.

## Working in this repo

Scripts are edited in the Alpine copy (`/pl/active/ZhuLab_BDC/GitHub/ISL1_project`), committed and pushed to
`Han-Z-BDC/ISL1_project` (private). This repo moves to the lab organization and becomes public at publication.

Until the switch-over, the assay jobs still run from the original repos under `/pl/active/ZhuLab_BDC/GitHub/`; their
configs point at those paths. New commits made there are brought in with
`git fetch <name> main && git merge -X subtree=<folder>/ <name>/main` (remote names:
`ISL1_KO_histone_cutNrun`, `ISL1_cutNrun`, `ISL1_scRNAseq`).
