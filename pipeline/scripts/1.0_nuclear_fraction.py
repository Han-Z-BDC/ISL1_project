#!/usr/bin/env python
"""
RNA QC step 1.0: nuclear fraction per called cell (DropletQC-style).

nuclear_fraction = intronic / (intronic + exonic) reads per cell barcode,
from Cell Ranger's per-read region tag (RE: E = exonic, N = intronic) in the
BAM. Counts primary, uniquely mapped (MAPQ 255), non-duplicate reads
(Cell Ranger flags UMI duplicates 0x400), so each count is one UMI.
Same measure as DropletQC::nuclear_fraction_tags (Muskovic & Powell 2021).

Why: droplets holding a piece of cytoplasm without a nucleus have almost no
intronic reads (NF near 0), and they can pass count/gene/%MT filters. Step
1.1 drops cells below min_nuclear_fraction (configs/qc_RNA.yaml) when that
is set. Very high NF (nuclei that lost their cytoplasm) is not filtered here.

Needs the BAM (Cell Ranger create-bam true). By default it is found next to
the filtered matrix after resolving symlinks: sample_alignments.bam (multi)
or possorted_genome_bam.bam (count).

Writes {outdir}/{sample}.nuclear_fraction.tsv (barcode, exonic, intronic,
nuclear_fraction) for every barcode in the filtered matrix.
"""
import argparse
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pysam


def find_bam(filtered_h5: str) -> Path:
    d = Path(os.path.realpath(filtered_h5)).parent
    for cand in (d / "sample_alignments.bam", d / "possorted_genome_bam.bam",
                 d.parent / "possorted_genome_bam.bam"):
        if cand.exists():
            return cand
    raise SystemExit(f"ERROR: no BAM found next to {d} -- pass --bam (Cell Ranger needs create-bam true)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--h5-path", required=True, help="Cell Ranger filtered_feature_bc_matrix.h5 (the called cells)")
    ap.add_argument("--bam", help="BAM with CB and RE tags; default: found next to the filtered matrix")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    with h5py.File(args.h5_path, "r") as f:
        bcs = [b.decode() for b in f["matrix"]["barcodes"][:]]
    idx = {b: i for i, b in enumerate(bcs)}
    bam = Path(args.bam) if args.bam else find_bam(args.h5_path)

    ex = np.zeros(len(bcs), np.int64)
    intr = np.zeros(len(bcs), np.int64)
    n = 0
    with pysam.AlignmentFile(str(bam), threads=args.threads) as f:
        for r in f.fetch(until_eof=True):
            n += 1
            if r.flag & 0xD04 or r.mapping_quality != 255:   # unmapped, secondary, duplicate, supplementary
                continue
            try:
                i = idx[r.get_tag("CB")]
                region = r.get_tag("RE")
            except KeyError:
                continue
            if region == "E":
                ex[i] += 1
            elif region == "N":
                intr[i] += 1

    df = pd.DataFrame({"barcode": bcs, "exonic": ex, "intronic": intr})
    tot = df["intronic"] + df["exonic"]
    df["nuclear_fraction"] = np.where(tot > 0, df["intronic"] / tot.clip(lower=1), np.nan).round(4)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / f"{args.sample_id}.nuclear_fraction.tsv", sep="\t", index=False)
    print(f"[{args.sample_id}] {bam}: {n} reads scanned, {len(bcs)} cells, "
          f"median nuclear fraction {df['nuclear_fraction'].median():.3f}")


if __name__ == "__main__":
    main()
