#!/usr/bin/env python3
"""Step 01: GREAT regulatory domains from the bulk GTF, and ISL1 site -> gene links.

GREAT "basal plus extension" (McLean et al. 2010), default settings:
  basal domain   5 kb upstream to 1 kb downstream of the TSS (strand-aware)
  extension      in both directions up to the nearest other gene's basal domain, at most 1 Mb from the TSS;
                 a domain never shrinks below its own basal domain
Genes: every gene record in the GTF on chr1-22/X/Y (all gene types), TSS = gene start (+) or end (-).
A site belongs to a gene when the site centre (the summit, start + 200 for the 401 bp sites) lies in the
gene's domain; a site can belong to more than one gene.

usage: python 01_great_site_to_gene.py <genes.gtf> <chrom_sizes> <sites.bed> <bulk_LRT_genotype.tsv> <out_dir>
Output in <out_dir>:
  gene_regions.bed  chr, start, end (domain, 0-based half-open), gene_id, gene_name, strand, tss, basal_start,
                    basal_end, gene_type, tested_in_bulk
  site_gene.tsv     site, chr, summit, gene_id, gene_name, gene_type, distance_to_tss (+ = downstream of the TSS
                    in the gene's direction), in_basal, tested_in_bulk
  summary.tsv       item, n
"""
import bisect
import csv
import heapq
import os
import re
import statistics
import sys

GTF, SIZES, SITES, BULK, OD = sys.argv[1:6]
UP, DOWN, EXT = 5000, 1000, 1_000_000
CANON = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
os.makedirs(OD, exist_ok=True)

chrom_len = {}
with open(SIZES) as f:
    for line in f:
        c, n = line.split()[:2]
        chrom_len[c] = int(n)

tested = set()
with open(BULK, newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        tested.add(r["gene_id"])

# ---- genes and basal domains
attr = re.compile(r'(gene_id|gene_name|gene_type) "([^"]*)"')
genes = {}
with open(GTF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[2] != "gene" or p[0] not in CANON:
            continue
        a = dict(attr.findall(p[8]))
        start0, end0, strand = int(p[3]) - 1, int(p[4]), p[6]
        tss = start0 if strand == "+" else end0 - 1
        if strand == "+":
            b0, b1 = tss - UP, tss + DOWN + 1
        else:
            b0, b1 = tss - DOWN, tss + UP + 1
        L = chrom_len[p[0]]
        genes.setdefault(p[0], []).append({
            "chr": p[0], "gene_id": a["gene_id"], "gene_name": a.get("gene_name", a["gene_id"]),
            "gene_type": a.get("gene_type", ""), "strand": strand, "tss": tss,
            "b0": max(0, b0), "b1": min(L, b1)})

# ---- extension
for c, gl in genes.items():
    L = chrom_len[c]
    by_start = sorted(gl, key=lambda g: (g["b0"], g["b1"]))
    starts = [g["b0"] for g in by_start]
    by_end = sorted(gl, key=lambda g: (g["b1"], g["b0"]))
    ends = [g["b1"] for g in by_end]
    for g in gl:
        # nearest basal domain ending at or before this basal start
        i = bisect.bisect_right(ends, g["b0"]) - 1
        left = ends[i] if i >= 0 else 0
        # nearest basal domain starting at or after this basal end
        j = bisect.bisect_left(starts, g["b1"])
        right = starts[j] if j < len(starts) else L
        g["d0"] = min(g["b0"], max(left, g["tss"] - EXT, 0))
        g["d1"] = max(g["b1"], min(right, g["tss"] + EXT + 1, L))

with open(os.path.join(OD, "gene_regions.bed"), "w") as f:
    for c in sorted(genes):
        for g in sorted(genes[c], key=lambda g: (g["d0"], g["d1"], g["gene_id"])):
            f.write("\t".join(map(str, [c, g["d0"], g["d1"], g["gene_id"], g["gene_name"], g["strand"], g["tss"],
                                         g["b0"], g["b1"], g["gene_type"], int(g["gene_id"] in tested)])) + "\n")

# ---- sites -> genes (sweep along each chromosome)
sites = {}
with open(SITES) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        s0, s1 = int(p[1]), int(p[2])
        sites.setdefault(p[0], []).append((s0 + (s1 - s0) // 2, p[3]))

rows = []
n_sites = n_assigned = 0
per_site = []
for c in sorted(sites):
    doms = sorted(genes.get(c, []), key=lambda g: g["d0"])
    k = 0
    active = []  # heap of (d1, gene_id, gene)
    for pos, name in sorted(sites[c]):
        n_sites += 1
        while k < len(doms) and doms[k]["d0"] <= pos:
            heapq.heappush(active, (doms[k]["d1"], doms[k]["gene_id"], doms[k]))
            k += 1
        while active and active[0][0] <= pos:
            heapq.heappop(active)
        hits = sorted((g for _, _, g in active), key=lambda g: g["gene_id"])
        per_site.append(len(hits))
        if hits:
            n_assigned += 1
        for g in hits:
            d = pos - g["tss"] if g["strand"] == "+" else g["tss"] - pos
            rows.append([name, c, pos, g["gene_id"], g["gene_name"], g["gene_type"], d,
                         int(g["b0"] <= pos < g["b1"]), int(g["gene_id"] in tested)])

rows.sort(key=lambda r: (r[0], r[3]))
with open(os.path.join(OD, "site_gene.tsv"), "w") as f:
    f.write("site\tchr\tsummit\tgene_id\tgene_name\tgene_type\tdistance_to_tss\tin_basal\ttested_in_bulk\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")

# ---- summary
n_genes = sum(len(v) for v in genes.values())
tested_in_gtf = {g["gene_id"] for v in genes.values() for g in v} & tested
sites_per_tested = {}
for r in rows:
    if r[8]:
        sites_per_tested[r[3]] = sites_per_tested.get(r[3], 0) + 1
summ = [
    ("genes in GTF (chr1-22/X/Y)", n_genes),
    ("genes tested in bulk", len(tested)),
    ("tested genes with a GREAT domain", len(tested_in_gtf)),
    ("tested genes with >= 1 site", len(sites_per_tested)),
    ("tested genes with >= 1 site, %", round(100 * len(sites_per_tested) / len(tested_in_gtf), 1)),
    ("sites per bound tested gene, median", statistics.median(sites_per_tested.values())),
    ("sites", n_sites),
    ("sites linked to >= 1 gene", n_assigned),
    ("sites linked to 1 gene", sum(1 for x in per_site if x == 1)),
    ("sites linked to 2 genes", sum(1 for x in per_site if x == 2)),
    ("sites linked to >= 3 genes", sum(1 for x in per_site if x >= 3)),
    ("site-gene links", len(rows)),
    ("site-gene links in a basal domain", sum(r[7] for r in rows)),
    ("|distance to TSS|, median bp", int(statistics.median(abs(r[6]) for r in rows))),
]
with open(os.path.join(OD, "summary.tsv"), "w") as f:
    f.write("item\tn\n")
    for k_, v in summ:
        f.write(f"{k_}\t{v}\n")
for k_, v in summ:
    print(f"{k_}\t{v}")
