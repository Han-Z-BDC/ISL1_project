#!/usr/bin/env python3
"""Library table for the browser tracks (scripts/23a; agreed 2026-11-09).

Usage: python3 tracks_manifest.py <sample_metadata.tsv> <ctr_replicates.tsv> <results_dir> <out.tsv>

One row per BAM: name, kind (ISL1 / H3K27ac / H3K27me3 / IgG), bam, bin size, groups (comma list of
the averaged tracks it goes into). Libraries:
  - batch 3 H3CTR / ISL1KO H3K27ac and H3K27me3 (filtered BAMs, duplicates kept, as for the calls);
    group <day>_<line>_<mark>
  - batch 3 H3CTR / ISL1KO block IgGs (dedup BAMs, the controls of the calls); group <day>_<line>_IgG
  - the 15 control ISL1 libraries of ctr_replicates.tsv, named by collection day (Ctr rep1-5); the
    chosen replicates go into <day>_Ctr_ISL1 (D15/D18 rep1,3,4; D22 rep2,4,5; D9 rep1,2) and their
    IgGs into <day>_Ctr_ISL1_IgG. A BAM used twice (batch 3 IgGs) gets one row with both groups.
Bins: 10 bp for ISL1, H3K27ac and IgG; 50 bp for H3K27me3.
"""
import csv
import os
import sys

meta, reps, res, out = sys.argv[1:5]
LINES = {"H3CTR", "ISL1KO"}
CHOSEN = {"D9": {"1", "2"}, "D15": {"1", "3", "4"}, "D18": {"1", "3", "4"}, "D22": {"2", "4", "5"}}
BIN = {"ISL1": 10, "H3K27ac": 10, "H3K27me3": 50, "IgG": 10}
rows = {}  # bam -> dict


def add(name, kind, bam, group):
    r = rows.setdefault(bam, {"name": name, "kind": kind, "bam": bam, "bin": BIN[kind], "groups": []})
    if group and group not in r["groups"]:
        r["groups"].append(group)


with open(meta) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["line"] not in LINES:
            continue
        g = f"{r['day']}_{r['line']}"
        if r["antibody"] in ("H3K27ac", "H3K27me3"):
            add(r["sample_name"], r["antibody"], os.path.join(res, "bam", r["sample_name"] + ".filtered.bam"), f"{g}_{r['antibody']}")
        elif r["antibody"] == "IgG":
            add(r["sample_name"], "IgG", os.path.join(res, "igg_background", r["sample_name"] + ".dedup.bam"), f"{g}_IgG")

with open(reps) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        chosen = r["replicate"] in CHOSEN.get(r["day"], set())
        add(r["sample"], "ISL1", r["isl1_bam"], f"{r['day']}_Ctr_ISL1" if chosen else "")
        add(r["sample"].replace("_ISL1", "_IgG"), "IgG", r["igg_bam"], f"{r['day']}_Ctr_ISL1_IgG" if chosen else "")

with open(out, "w") as fh:
    fh.write("name\tkind\tbam\tbin\tgroups\n")
    order = {"ISL1": 0, "H3K27ac": 1, "H3K27me3": 2, "IgG": 3}
    for r in sorted(rows.values(), key=lambda x: (order[x["kind"]], x["name"])):
        fh.write(f"{r['name']}\t{r['kind']}\t{r['bam']}\t{r['bin']}\t{','.join(r['groups'])}\n")
print(f"{len(rows)} libraries -> {out}")
