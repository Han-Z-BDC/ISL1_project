#!/usr/bin/env python3
"""Build samples.tsv (pipeline format) and sample_metadata.tsv from the Admera delivery folder.

Library IDs are the fastq prefixes (<library>_R1_001.fastq.gz). Each sample block
(day, line, replicate) has one IgG, which is the control_sample of every other antibody in
that block. Two antibody names in the delivery have typos and are normalized:
H327Kac -> H3K27ac (D15 ISL1-KO), H3K27Me -> H3K27me3 (D22 H3-CTR).

usage: make_samples.py <delivery_dir> <batch_folder_name> <samples.tsv> <sample_metadata.tsv>
"""
import re
import sys
from pathlib import Path

LINES = {"H3-CTR": "H3CTR", "ISL1-KO": "ISL1KO", "Ctr-GFP": "CtrGFP", "KO-ISL1": "KOISL1"}
ANTIBODIES = {"ISL1": "ISL1", "IgG": "IgG", "Flag": "Flag", "H3K27ac": "H3K27ac", "H327Kac": "H3K27ac",
              "H3K27me3": "H3K27me3", "H3K27Me": "H3K27me3"}
PAT = re.compile(r"^(D\d+)-(.+)-rep-(\d+)_([^_]+)_(S\d+_L\d+)$")


def main(delivery, batch, samples_tsv, metadata_tsv):
    rows = []
    for f in sorted(Path(delivery).glob("*_R1_001.fastq.gz")):
        lib = f.name[: -len("_R1_001.fastq.gz")]
        m = PAT.match(lib)
        if not m:
            sys.exit(f"ERROR: unexpected library name {lib}")
        day, line, rep, ab, _ = m.groups()
        if line not in LINES or ab not in ANTIBODIES:
            sys.exit(f"ERROR: unknown line or antibody in {lib}")
        if not (Path(delivery) / f"{lib}_R2_001.fastq.gz").exists():
            sys.exit(f"ERROR: no R2 for {lib}")
        block = f"{day}_{LINES[line]}_rep{rep}"
        rows.append({"library_id": lib, "block": block, "day": day, "line": LINES[line], "replicate": rep,
                     "antibody": ANTIBODIES[ab], "sample_name": f"{block}_{ANTIBODIES[ab]}"})

    igg = {r["block"]: r["sample_name"] for r in rows if r["antibody"] == "IgG"}
    names = [r["sample_name"] for r in rows]
    if len(set(names)) != len(names):
        sys.exit("ERROR: duplicate sample names after normalization")
    for r in rows:
        if r["antibody"] != "IgG" and r["block"] not in igg:
            sys.exit(f"ERROR: no IgG in block {r['block']}")

    order = {"IgG": 0, "ISL1": 1, "Flag": 2, "H3K27ac": 3, "H3K27me3": 4}
    rows.sort(key=lambda r: (int(r["day"][1:]), r["line"], int(r["replicate"]), order[r["antibody"]]))
    with open(samples_tsv, "w") as f:
        f.write("library_id\tsample_name\trole\tcontrol_sample\treplicate\n")
        for r in rows:
            is_igg = r["antibody"] == "IgG"
            f.write(f"{r['library_id']}\t{r['sample_name']}\t{'IgG' if is_igg else 'TF'}\t"
                    f"{'NA' if is_igg else igg[r['block']]}\t{r['replicate']}\n")
    with open(metadata_tsv, "w") as f:
        f.write("sample_name\tlibrary_id\tday\tline\treplicate\tantibody\tbatch_folder\n")
        for r in rows:
            f.write(f"{r['sample_name']}\t{r['library_id']}\t{r['day']}\t{r['line']}\t{r['replicate']}\t"
                    f"{r['antibody']}\t{batch}\n")
    print(f"{len(rows)} libraries, {len(igg)} blocks -> {samples_tsv}, {metadata_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
