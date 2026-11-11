# ISL1 scRNA-seq (10x OCM)

Single-cell RNA-seq for the ISL1 stem-cell islet differentiation project, run on the
10x GEM-X **On-Chip Multiplexing (OCM)** platform: up to 4 samples per chip channel,
each tagged with an OCM barcode (OB1-OB4) and demultiplexed by Cell Ranger.

Data: Admera Health project 260924-0015, delivered 2026-10-08 to
`/pl/active/ZhuLab_BDC/raw_data/ISL1_scRNAseq/7535-0002-10082026_215649/`.
That delivery also contains ISL1 CUT&RUN fastqs; this repo only uses the `*_scRNA_*` files.

| Channel (fastq_id) | OCM samples |
|---|---|
| `HZ_Diff34_D40_scRNA` | OB1 Diff34_D40_ctr_1, OB2 Diff34_D40_KO_1, OB3 Diff34_D40_ctr_2, OB4 Diff34_D40_KO_2 (order confirmed from data 2026-10-10: OB2/OB4 show the ISL1-KO signature -- INS/GCG/SST/ARX lost, TPH1/FEV up; OB1~OB3, OB2~OB4) |

## Setup

- `config.sh` -- all paths (raw data, reference, Cell Ranger binary, list of channels).
- `configs/<run_id>.multi.csv` -- one `cellranger multi` config per OCM channel. The
  `[samples]` section maps OCM barcodes to sample names; fill it in before running.
- Reference: `refdata-cellranger-arc-GRCh38-2020-A-2.0.0`, same as the diff-multiome
  data, so genes match across datasets.
- Cell Ranger 9.0.1 (OCM needs >= 9.0), at `/projects/hzhu10@xsede.org/software/cellranger-9.0.1`.

## Pipeline

0a. `scripts/00a_admera_download.sb` -- FTP download from Admera (credentials via env vars, not in git).
0. `scripts/00_md5_check.sb` -- verify the Admera delivery against its `md5.txt`.
1. `scripts/01_cellranger_multi.sb` -- `cellranger multi` per channel on `amem`.
   ```
   sbatch scripts/01_cellranger_multi.sb
   ```
   Outputs (on scratch): `/scratch/alpine/hzhu10@xsede.org/ISL1_scRNAseq/cellranger/<run_id>/outs/per_sample_outs/<sample_id>/`
   (`web_summary.html`, `count/sample_filtered_feature_bc_matrix.h5`, BAM).
