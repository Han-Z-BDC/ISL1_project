#!/usr/bin/env python3
"""Step 10a: HOMER inputs per day, from the macs3_ge2 consensus peaks (step 05).

Target: every macs3_ge2 region of the day, centred on its strongest summit (max -log10 p) among the
day's three libraries' step 04 MACS3 peaks inside the region (same rule as step 06), summit +/- half
width (200 bp for half width 100). chr1-22/X/Y; windows touching the blacklist are dropped.
Other-days background: the target windows of the other two days, minus every window that overlaps a
consensus region of the target day (peaks shared with the target day would cancel out); windows of
the two days that overlap each other are kept once (the stronger summit; ties by position).
  python homer_inputs.py <OUT_DIR> <sample_metadata.tsv> <bedtools> <blacklist> <genome.txt> <half_width> <out_dir>
Output: <out_dir>/<day>.target.bed, <out_dir>/<day>.other_days.bed (BED6, strand +), <out_dir>/qc.tsv
"""
import bisect
import csv
import os
import subprocess
import sys

OUT, META, BT, BLACKLIST, GENOME, HALF, OD = sys.argv[1:8]
HALF = int(HALF)
DAYS = ["D15", "D18", "D22"]
CANON = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
os.makedirs(OD, exist_ok=True)

chrom_len = {}
with open(GENOME) as f:
    for line in f:
        c, n = line.split()[:2]
        chrom_len[c] = int(n)

with open(META, newline="") as f:
    libs = [r for r in csv.DictReader(f, delimiter="\t") if r["antibody"] == "ISL1"]


def read_regions(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] in CANON:
                out.setdefault(p[0], []).append((int(p[1]), int(p[2]), p[3]))
    for c in out:
        out[c].sort()
    return out


def overlaps(regs, st, en):
    """regs: sorted (start, end, name) of one chromosome, non-overlapping (merged consensus)."""
    i = bisect.bisect_left(regs, (en,))
    return any(rs < en and re > st for rs, re, _ in regs[max(0, i - 50):i])


qc = []
regions, windows = {}, {}
for day in DAYS:
    regions[day] = read_regions(os.path.join(OUT, "consensus", f"{day}_Ctr_ISL1.macs3_ge2.bed"))
    summits = {}
    for r in libs:
        if r["day"] == day:
            with open(os.path.join(OUT, "macs3", f"{r['sample_name']}_peaks.filtered.narrowPeak")) as f:
                for line in f:
                    p = line.split("\t")
                    if p[0] in CANON:
                        summits.setdefault(p[0], []).append((int(p[1]) + int(p[9]), float(p[7])))
    for c in summits:
        summits[c].sort()
    win, n_reg, n_skip = [], 0, 0
    for c, regs in regions[day].items():
        pos = [s[0] for s in summits.get(c, [])]
        for start, end, name in regs:
            n_reg += 1
            i, j = bisect.bisect_left(pos, start), bisect.bisect_left(pos, end)
            if i == j:
                n_skip += 1
                continue
            # strongest summit; ties go to the leftmost position (deterministic)
            s, score = max(summits[c][i:j], key=lambda x: (x[1], -x[0]))
            if s - HALF < 0 or s + HALF > chrom_len.get(c, 0):
                n_skip += 1
                continue
            win.append((c, s - HALF, s + HALF, name, score))
    win.sort(key=lambda x: (x[0], x[1]))
    tmp = os.path.join(OD, f"{day}.tmp.bed")
    with open(tmp, "w") as f:
        for c, st, en, name, score in win:
            f.write(f"{c}\t{st}\t{en}\t{name}\t{score:.4f}\t+\n")
    r = subprocess.run([BT, "intersect", "-v", "-a", tmp, "-b", BLACKLIST], capture_output=True, text=True, check=True)
    os.remove(tmp)
    kept = [l.split("\t") for l in r.stdout.splitlines() if l.strip()]
    windows[day] = [(p[0], int(p[1]), int(p[2]), p[3], float(p[4])) for p in kept]
    with open(os.path.join(OD, f"{day}.target.bed"), "w") as f:
        f.write(r.stdout)
    qc += [(f"{day} macs3_ge2 regions", n_reg), (f"{day} no summit / past chromosome end", n_skip),
           (f"{day} blacklist", len(win) - len(kept)), (f"{day} target windows", len(kept))]

for day in DAYS:
    others = [w for d in DAYS if d != day for w in windows[d]]
    n_all = len(others)
    others = [w for w in others if not overlaps(regions[day].get(w[0], []), w[1], w[2])]
    n_shared = n_all - len(others)
    # strongest first; overlapping windows of the two other days kept once
    others.sort(key=lambda x: (-x[4], x[0], x[1], x[3]))
    acc, bg = {}, []
    for w in others:
        a = acc.setdefault(w[0], [])
        i = bisect.bisect_right(a, w[1] - 2 * HALF)
        if i < len(a) and a[i] < w[2]:
            continue
        bisect.insort(a, w[1])
        bg.append(w)
    bg.sort(key=lambda x: (x[0], x[1]))
    with open(os.path.join(OD, f"{day}.other_days.bed"), "w") as f:
        for c, st, en, name, score in bg:
            f.write(f"{c}\t{st}\t{en}\t{name}\t{score:.4f}\t+\n")
    qc += [(f"{day} other-days windows", n_all), (f"{day} other-days, overlap a {day} region (dropped)", n_shared),
           (f"{day} other-days, duplicate between the two days (dropped)", len(others) - len(bg)),
           (f"{day} other-days background", len(bg))]

with open(os.path.join(OD, "qc.tsv"), "w") as f:
    f.write("item\tn\n")
    for k, v in qc:
        f.write(f"{k}\t{v}\n")
for k, v in qc:
    print(f"{k}\t{v}")
