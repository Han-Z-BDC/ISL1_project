# STAR index choice

`/pl/active/ZhuLab_BDC/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/` has two STAR
indices. This pipeline uses **`STARindex_2.7.5c/`**, not `star/`.

| | `STARindex_2.7.5c/` (used) | `star/` (not used) |
|---|---|---|
| Built by | This lab, on this cluster | 10x Genomics, on their own infra (path in `genomeParameters.txt` is `/mnt/deck7/rudy/arc_refs/...`) |
| Purpose | General-purpose standalone STAR index | Bundled inside the Cell Ranger ARC reference tarball, for ARC's own internal RNA step |
| STAR genome format version | `2.7.4a` | `2.7.1a` |
| `--sjdbOverhang` | `99` (built assuming ~100bp reads) | not set |

Both are built from the same underlying GRCh38-2020-A fasta/GTF, so gene models are
identical either way — the difference is provenance/intent, not annotation content.
`STARindex_2.7.5c` is the one actually built for standalone alignment.

STAR binary: use the `bulkRNAseq` conda env's `STAR 2.7.5c` (pinned to match the index's
build version exactly). Also available on this cluster: `schicar` env (STAR 2.7.5c) and
`pairedtag_env` (STAR 2.7.10b, newer, also index-compatible).
