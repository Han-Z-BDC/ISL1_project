#!/usr/bin/env python3
"""Step 06: fixed-width ISL1 site set for differential binding over time.

Candidates: every macs3_ge2 region of D15, D18 and D22 (step 05), centred on the strongest summit
(max -log10 p) among that day's three libraries' step 04 MACS3 peaks inside the region. Sites are summit - 200 .. summit + 200 (401 bp). Overlapping or book-ended
candidates across days are resolved by keeping the stronger summit, iteratively (Corces et al.
2018). chr1-22/X/Y only; sites touching the blacklist are dropped.
  python diff_sites.py <OUT_DIR> <sample_metadata.tsv> <bedtools> <blacklist> <genome.txt> <half_width> <out_dir>
Output: <out_dir>/sites.bed (chr, start, end, name, -log10 p, strand, summit_day, days_called,
n_days), <out_dir>/qc.tsv
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


qc = []
cands = []
day_regions = {}
for day in DAYS:
    regions = read_regions(os.path.join(OUT, "consensus", f"{day}_Ctr_ISL1.macs3_ge2.bed"))
    day_regions[day] = regions
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
    n_reg = n_nosummit = 0
    for c, regs in regions.items():
        pos = [s[0] for s in summits.get(c, [])]
        for start, end, name in regs:
            n_reg += 1
            i, j = bisect.bisect_left(pos, start), bisect.bisect_left(pos, end)
            if i == j:
                n_nosummit += 1
                continue
            # strongest summit; ties go to the leftmost position (deterministic)
            s, score = max(summits[c][i:j], key=lambda x: (x[1], -x[0]))
            cands.append((score, c, s, day, name))
    qc.append((f"{day} macs3_ge2 regions", n_reg))
    if n_nosummit:
        qc.append((f"{day} regions without a summit (skipped)", n_nosummit))

qc.append(("candidates, all days", len(cands)))
# strongest first; ties broken by chromosome, position, day (deterministic)
cands.sort(key=lambda x: (-x[0], x[1], x[2], x[3]))
accepted = {}
kept = []
n_edge = n_overlap = 0
for score, c, s, day, name in cands:
    if s - HALF < 0 or s + HALF + 1 > chrom_len.get(c, 0):
        n_edge += 1
        continue
    a = accepted.setdefault(c, [])
    i = bisect.bisect_left(a, s - 2 * HALF - 1)
    # conflict: an accepted summit within 401 bp (overlapping or book-ended 401 bp windows)
    if i < len(a) and a[i] <= s + 2 * HALF + 1:
        n_overlap += 1
        continue
    bisect.insort(a, s)
    kept.append((c, s - HALF, s + HALF + 1, score, day))
qc.append(("removed: overlap a stronger summit", n_overlap))
if n_edge:
    qc.append(("removed: window past chromosome end", n_edge))

kept.sort(key=lambda x: (x[0], x[1]))
tmp = os.path.join(OD, "sites.tmp.bed")
with open(tmp, "w") as f:
    for c, st, en, score, day in kept:
        f.write(f"{c}\t{st}\t{en}\t{score}\t{day}\n")
r = subprocess.run([BT, "intersect", "-v", "-a", tmp, "-b", BLACKLIST], capture_output=True, text=True, check=True)
final = [l.split("\t") for l in r.stdout.splitlines() if l.strip()]
os.remove(tmp)
qc.append(("removed: blacklist", len(kept) - len(final)))


def called(day, c, st, en):
    regs = day_regions[day].get(c, [])
    i = bisect.bisect_left(regs, (en,))
    return any(rs < en and re > st for rs, re, _ in regs[max(0, i - 50):i])


patterns = {}
with open(os.path.join(OD, "sites.bed"), "w") as f:
    for k, (c, st, en, score, day) in enumerate(final, 1):
        st, en = int(st), int(en)
        days = [d for d in DAYS if called(d, c, st, en)]
        key = ",".join(days)
        patterns[key] = patterns.get(key, 0) + 1
        f.write(f"{c}\t{st}\t{en}\tisl1_site_{k:05d}\t{float(score):.4f}\t.\t{day.strip()}\t{key}\t{len(days)}\n")
qc.append(("sites, final", len(final)))
for key in sorted(patterns, key=lambda k: (-len(k.split(",")), k)):
    qc.append((f"sites called on {key}", patterns[key]))

with open(os.path.join(OD, "qc.tsv"), "w") as f:
    f.write("item\tn\n")
    for k, v in qc:
        f.write(f"{k}\t{v}\n")
for k, v in qc:
    print(f"{k}\t{v}")
