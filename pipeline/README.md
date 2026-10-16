# CUT&RUN pipeline

FASTQ to peaks for CUT&RUN: QC → adapter trimming → alignment → an IgG
background of your choice → peak calling → CPM bigwig tracks.

This repo is self-contained — clone it into your own folder (`/projects/$USER/...`
or anywhere else you have write access) and run it there. Nothing in it
points at any specific person's or lab's storage; every path you actually
need is set once, in your own `config.sh`, which you create from the
template and never commit.

## 1. Set up the environment

```bash
acompile --ntasks=4 --time=01:00:00      # build environments on a compute node, not the login node
module load miniforge
conda env create -f environment.yml      # builds an env named "cutnrun"
conda activate cutnrun
```

`config.sh` (next step) points `TOOLS` at this environment's `bin/` directory
so every pipeline script calls tools by full path rather than relying on
`conda activate` inside a batch job.

## 2. Configure

```bash
cp config.example.sh config.sh
nano config.sh
```

Fill in, at minimum: `RAW_DIR`, `OUT_DIR`, `BOWTIE2_INDEX`, `BLACKLIST`, and
`TOOLS` (the path `conda env create` just built — find it with
`conda env list`). `config.sh` is gitignored on purpose; it's specific to
your account and your project.

## 3. Describe your samples

```bash
cp samples.example.tsv samples.tsv
```

One row per library. `role` is `TF` (your target — ISL1, or whatever
transcription factor/mark you're profiling) or `IgG`. Every `TF` row's
`control_sample` must name an `IgG` row in the same file.

## 4. Pick an IgG background method

Set `IGG_METHOD` in `config.sh` to one of:

| | `matched` | `pooled` |
|---|---|---|
| What it does | Each target sample uses its own paired IgG control | Several IgG replicates are deduplicated, downsampled to a common depth, and merged into one shared background |
| When to use it | Default — works with a single decent IgG replicate per target | You have ≥2 usable IgG replicates and want a more robust, less noise-prone background |
| Extra config needed | None | `POOL_IGG_SAMPLES`, `POOL_IGG_SEEDS`, `POOL_TARGET_FRAGMENTS` |

Either way, whatever IgG BAM ends up as the MACS3 control is deduplicated —
target (TF) samples always keep their duplicates, all the way through to
peak calling (`--keep-dup all`). CUT&RUN libraries start from little input
material, so at a true binding site, seeing the same fragment repeatedly is
real signal, not just PCR noise — removing duplicates from the target risks
flattening your strongest peaks.

## 5. Run it, in order

```bash
sbatch scripts/00_fastqc.sb              # wait for it to finish, check qc/raw_fastqc/multiqc_report.html
bash scripts/submit_align.sh             # array job, one task per sample -- don't sbatch 01_align.sb directly
sbatch scripts/02_build_igg_background.sb  # builds whichever background IGG_METHOD picked
sbatch scripts/03_call_macs3.sb          # peaks: target (dup-retained) vs background (deduplicated)
sbatch scripts/04_tracks.sb              # CPM bigwigs for every sample, for a genome browser
```

Chain them with `--dependency=afterok:<jobid>` instead of babysitting the
queue between steps — see the lab manual's Dependencies section.

## Output layout

```
$OUT_DIR/
├── qc/            fastqc, multiqc, per-sample cutadapt/bowtie2/markdup logs
├── trimmed/       adapter-trimmed FASTQs
├── bam/           filtered, duplicate-marked BAMs -- every sample
├── igg_background/  deduplicated IgG BAM(s), method depends on IGG_METHOD
├── macs3/         narrowPeak / summits, raw and blacklist-filtered
└── tracks/bigwig/ CPM-normalized bigwigs, every sample
```

## Before trusting your peaks

- Check `qc/raw_fastqc/multiqc_report.html` first — adapter content and
  duplication-rate flags are expected for CUT&RUN, but check quality scores
  and overall alignment rate look healthy.
- Compare target vs. IgG signal at a few known-positive and known-negative
  loci in a genome browser before doing anything downstream with the calls.
- `matched` and `pooled` will not give identical peak sets — that's
  expected, not a bug. If you switch methods on the same data, treat it as
  a real methodological change worth noting in your analysis, not a rerun
  of "the same" pipeline.
