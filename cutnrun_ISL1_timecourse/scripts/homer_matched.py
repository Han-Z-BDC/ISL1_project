#!/usr/bin/env python3
"""Step 11a helper: HOMER inputs with backgrounds matched to the peaks' genomic composition.

  python homer_matched.py tss <genes.gtf> <genome.txt> <n> <half_width> <out.bed>
      Promoter pool: n windows of 2 x half_width bp, each centred at a random offset in [-1000, +1000]
      bp of a random TSS (unique transcript TSSs of the GTF, chr1-22/X/Y; seed 3). Filtering
      (blacklist, N, peaks) is done by the caller.
  python homer_matched.py select <homer_dir> <pool_promoter.bed> <pool_genome.bed> <pool_genome.annotatePeaks.txt> <out_dir>
      Per day, from the step 10a target windows and their step 10c annotation:
        <day>.promoter.bed   targets within 1 kb of a TSS (HOMER "Distance to TSS")
        <day>.distal.bed     the other targets
        <day>.bg_promoter.bed  promoter pool (A: promoter peaks vs promoters)
        <day>.bg_distal.bed    genome pool windows > 1 kb from a TSS, up to 100,000 (A: distal vs genome)
        <day>.bg_matched.bed   B: 3 x the targets of each class, drawn from the pools: promoter (within
                               1 kb of a TSS) from the promoter pool, the rest by HOMER category
                               (intron, Intergenic, exon, TTS, other) from the genome pool
      and qc.tsv (counts per class for targets and backgrounds). Draws use seed 11 (deterministic).
"""
import os
import random
import sys

DAYS = ["D15", "D18", "D22"]
CANON = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
CLASSES = ["promoter", "intron", "Intergenic", "exon", "TTS", "other"]
FOLD = 3
MAX_DISTAL = 100000


def write_bed(path, rows):
    rows = sorted(rows, key=lambda r: (r[0], r[1], r[3]))
    with open(path, "w") as f:
        for c, s, e, n in rows:
            f.write(f"{c}\t{s}\t{e}\t{n}\t0\t+\n")


def read_bed(path):
    with open(path) as f:
        return [(p[0], int(p[1]), int(p[2]), p[3]) for p in (l.rstrip("\n").split("\t") for l in f) if p[0]]


def read_annot(path):
    """PeakID -> (class, |distance to TSS|)."""
    out = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        ia, idist = head.index("Annotation"), head.index("Distance to TSS")
        for line in f:
            p = line.rstrip("\n").split("\t")
            d = abs(int(p[idist])) if p[idist] not in ("", "NA") else 10 ** 9
            cat = p[ia].split(" (")[0].strip()
            if d <= 1000:
                cls = "promoter"
            elif cat in ("intron", "Intergenic", "exon", "TTS"):
                cls = cat
            else:
                cls = "other"
            out[p[0]] = (cls, d)
    return out


if sys.argv[1] == "tss":
    gtf, genome, n, half, out = sys.argv[2:7]
    n, half = int(n), int(half)
    clen = {}
    with open(genome) as f:
        for line in f:
            c, l = line.split()[:2]
            clen[c] = int(l)
    tss = set()
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if p[2] == "transcript" and p[0] in CANON:
                tss.add((p[0], int(p[3]) - 1 if p[6] == "+" else int(p[4]) - 1))
    tss = sorted(tss)
    rng = random.Random(3)
    rows, seen = [], set()
    while len(rows) < n:
        c, t = tss[rng.randrange(len(tss))]
        centre = t + rng.randint(-1000, 1000)
        s, e = centre - half, centre + half
        if s < 0 or e > clen.get(c, 0) or (c, s) in seen:
            continue
        seen.add((c, s))
        rows.append((c, s, e, f"tsspool_{len(rows) + 1}"))
    write_bed(out, rows)
    print(f"{len(tss)} unique TSSs; {len(rows)} promoter-pool windows")
    sys.exit(0)

H, POOL_P, POOL_G, POOL_G_ANN, OD = sys.argv[2:7]
os.makedirs(OD, exist_ok=True)
rng = random.Random(11)
pool_p = read_bed(POOL_P)
g_ann = read_annot(POOL_G_ANN)
pool_g = {c: [] for c in CLASSES}
for r in read_bed(POOL_G):
    cls = g_ann.get(r[3], ("other", 0))[0]
    pool_g[cls].append(r)
qc = [("pool promoter windows", len(pool_p))] + [(f"pool genome windows, {c}", len(pool_g[c])) for c in CLASSES]

distal_all = [r for c in CLASSES if c != "promoter" for r in pool_g[c]]
bg_distal = sorted(rng.sample(distal_all, min(MAX_DISTAL, len(distal_all))))
qc.append(("background distal (A, all days)", len(bg_distal)))

for day in DAYS:
    tg = read_bed(os.path.join(H, "inputs", f"{day}.target.bed"))
    ann = read_annot(os.path.join(H, "annotation", f"{day}.annotatePeaks.txt"))
    by = {c: [] for c in CLASSES}
    for r in tg:
        by[ann[r[3]][0]].append(r)
    write_bed(os.path.join(OD, f"{day}.promoter.bed"), by["promoter"])
    write_bed(os.path.join(OD, f"{day}.distal.bed"), [r for c in CLASSES if c != "promoter" for r in by[c]])
    write_bed(os.path.join(OD, f"{day}.bg_promoter.bed"), pool_p)
    write_bed(os.path.join(OD, f"{day}.bg_distal.bed"), bg_distal)
    matched = []
    for c in CLASSES:
        src = pool_p if c == "promoter" else pool_g[c]
        want = FOLD * len(by[c])
        got = rng.sample(src, min(want, len(src)))
        matched += got
        qc.append((f"{day} targets, {c}", len(by[c])))
        qc.append((f"{day} matched background, {c}", len(got)))
        if len(got) < want:
            qc.append((f"{day} matched background, {c}: pool short by", want - len(got)))
    write_bed(os.path.join(OD, f"{day}.bg_matched.bed"), matched)
    qc += [(f"{day} targets, promoter", len(by["promoter"])),
           (f"{day} targets, distal", len(tg) - len(by["promoter"])),
           (f"{day} matched background, total", len(matched))]

with open(os.path.join(OD, "qc.tsv"), "w") as f:
    f.write("item\tn\n")
    seen = set()
    for k, v in qc:
        if k in seen:
            continue
        seen.add(k)
        f.write(f"{k}\t{v}\n")
for k, v in qc:
    print(f"{k}\t{v}")
