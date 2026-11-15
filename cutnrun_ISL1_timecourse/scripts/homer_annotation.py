#!/usr/bin/env python3
"""Step 10c: summary of HOMER annotatePeaks.pl output per set (D15, D18, D22 consensus peaks; genome
baseline of random windows): count and percent per HOMER category, and percent within 1 kb of a TSS.
  python homer_annotation.py <annotation_dir>
Output: <annotation_dir>/summary.tsv (set, category, n, pct)
"""
import csv
import os
import sys

A = sys.argv[1]
SETS = ["D15", "D18", "D22", "genome_random"]
CATS = ["promoter-TSS", "5' UTR", "exon", "intron", "3' UTR", "TTS", "non-coding", "Intergenic"]
rows = []
for s in SETS:
    with open(os.path.join(A, f"{s}.annotatePeaks.txt"), newline="") as f:
        r = csv.reader(f, delimiter="\t")
        head = next(r)
        ia, idist = head.index("Annotation"), head.index("Distance to TSS")
        counts, n, near = {}, 0, 0
        for p in r:
            n += 1
            cat = p[ia].split(" (")[0].strip() if p[ia] else "NA"
            counts[cat] = counts.get(cat, 0) + 1
            if p[idist] not in ("", "NA") and abs(int(p[idist])) <= 1000:
                near += 1
    for c in CATS + sorted(k for k in counts if k not in CATS):
        if c in counts or c in CATS:
            rows.append((s, c, counts.get(c, 0), 100 * counts.get(c, 0) / n))
    rows.append((s, "within 1 kb of a TSS", near, 100 * near / n))
    rows.append((s, "total", n, 100.0))
with open(os.path.join(A, "summary.tsv"), "w") as f:
    f.write("set\tcategory\tn\tpct\n")
    for s, c, k, p in rows:
        f.write(f"{s}\t{c}\t{k}\t{p:.2f}\n")
cats = [c for s, c, _, _ in rows if s == "D15"]
print("category\t" + "\t".join(SETS))
for c in cats:
    print(c + "\t" + "\t".join(f"{next(p for s2, c2, _, p in rows if s2 == s and c2 == c):.1f}"
                                if any(s2 == s and c2 == c for s2, c2, _, _ in rows) else "-" for s in SETS))
