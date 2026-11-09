#!/usr/bin/env python3
"""QC for the histone differential counts of scripts/16: PCA and correlation.

Usage: python3 diff_qc.py <results/diff/<mark>>

Reads samples.tsv, counts.txt (+ .summary) and filter.tsv. Uses regions passing the pre-filter.
Values are log2(CPM + 1), where CPM uses all fragments featureCounts saw for that library
(assigned + unassigned). Writes qc/{library_stats.tsv,pca.tsv,correlation.tsv,pca.png,
correlation.png,report.html}.
"""
import base64
import html
import io
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage

D = sys.argv[1]
Q = os.path.join(D, "qc")
TOP_VAR = 5000
DAY_COLOR = {"D9": "#4C72B0", "D15": "#DD8452", "D18": "#55A868", "D22": "#C44E52"}
MARKER = {"H3CTR": "o", "ISL1KO": "^"}

with open(os.path.join(D, "samples.tsv")) as fh:
    rows = [l.rstrip("\n").split("\t") for l in fh][1:]
meta = {os.path.basename(r[5]).replace(".filtered.bam", ""): dict(day=r[1], genotype=r[2], rep=r[3]) for r in rows}

with open(os.path.join(D, "counts.txt")) as fh:
    lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
names = [os.path.basename(c).replace(".filtered.bam", "") for c in lines[0][6:]]
ids = [r[0] for r in lines[1:]]
counts = np.array([[float(x) for x in r[6:]] for r in lines[1:]])
with open(os.path.join(D, "counts.txt.summary")) as fh:
    summ = [l.rstrip("\n").split("\t") for l in fh][1:]
total = np.array([sum(float(r[i]) for r in summ) for i in range(1, len(names) + 1)])
assigned = np.array([float(r[i]) for r in summ if r[0] == "Assigned" for i in range(1, len(names) + 1)])
with open(os.path.join(D, "filter.tsv")) as fh:
    keep = {l.split("\t")[0] for l in list(fh)[1:] if l.rstrip("\n").split("\t")[2] == "1"}
mask = np.array([i in keep for i in ids])

with open(os.path.join(Q, "library_stats.tsv"), "w") as out:
    out.write("sample\tday\tgenotype\treplicate\tfragments\tin_regions\tfrac_in_regions\tin_filtered_regions\n")
    for j, s in enumerate(names):
        m = meta[s]
        out.write(f"{s}\t{m['day']}\t{m['genotype']}\t{m['rep']}\t{int(total[j])}\t{int(assigned[j])}\t"
                  f"{assigned[j] / total[j]:.4f}\t{int(counts[mask, j].sum())}\n")

x = np.log2(counts[mask] / total * 1e6 + 1)
var = x.var(axis=1)
top = x[np.argsort(var)[::-1][:TOP_VAR]]
c = top - top.mean(axis=1, keepdims=True)
u, sv, vt = np.linalg.svd(c, full_matrices=False)
pcs = vt[:3].T * sv[:3]
expl = sv ** 2 / (sv ** 2).sum()
with open(os.path.join(Q, "pca.tsv"), "w") as out:
    out.write("sample\tday\tgenotype\treplicate\tPC1\tPC2\tPC3\n")
    for j, s in enumerate(names):
        m = meta[s]
        out.write(f"{s}\t{m['day']}\t{m['genotype']}\t{m['rep']}\t" + "\t".join(f"{v:.3f}" for v in pcs[j]) + "\n")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
for ax, (a, b) in zip(axes, [(0, 1), (0, 2)]):
    for j, s in enumerate(names):
        m = meta[s]
        ax.scatter(pcs[j, a], pcs[j, b], c=DAY_COLOR[m["day"]], marker=MARKER.get(m["genotype"], "s"),
                   s=70, edgecolors="black", linewidths=0.5)
        ax.annotate(f"r{m['rep']}", (pcs[j, a], pcs[j, b]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel(f"PC{a + 1} ({100 * expl[a]:.1f}%)")
    ax.set_ylabel(f"PC{b + 1} ({100 * expl[b]:.1f}%)")
handles = [plt.Line2D([], [], marker="o", ls="", color=v, label=k) for k, v in DAY_COLOR.items()
           if k in {m["day"] for m in meta.values()}]
handles += [plt.Line2D([], [], marker=v, ls="", color="grey", label=k) for k, v in MARKER.items()]
axes[1].legend(handles=handles, loc="best", fontsize=8)
fig.suptitle(f"PCA, top {TOP_VAR} variable regions, log2(CPM + 1)")
fig.tight_layout()
fig.savefig(os.path.join(Q, "pca.png"), dpi=130)
pca_png = fig

r = np.corrcoef(x.T)
order = leaves_list(linkage(x.T, "average", metric="correlation"))
with open(os.path.join(Q, "correlation.tsv"), "w") as out:
    out.write("sample\t" + "\t".join(names) + "\n")
    for i, s in enumerate(names):
        out.write(s + "\t" + "\t".join(f"{v:.3f}" for v in r[i]) + "\n")
fig2, ax = plt.subplots(figsize=(7.5, 6.5))
im = ax.imshow(r[np.ix_(order, order)], cmap="viridis")
lab = [names[i].replace(f"_{os.path.basename(D)}", "") for i in order]
ax.set_xticks(range(len(lab)), lab, rotation=90, fontsize=8)
ax.set_yticks(range(len(lab)), lab, fontsize=8)
for i in range(len(lab)):
    for j in range(len(lab)):
        ax.text(j, i, f"{r[order[i], order[j]]:.2f}", ha="center", va="center", fontsize=6, color="white")
fig2.colorbar(im, ax=ax, label="Pearson r")
ax.set_title("Pearson correlation, all filtered regions (clustered)")
fig2.tight_layout()
fig2.savefig(os.path.join(Q, "correlation.png"), dpi=130)


def b64(f):
    buf = io.BytesIO()
    f.savefig(buf, format="png", dpi=110)
    return base64.b64encode(buf.getvalue()).decode()


with open(os.path.join(Q, "library_stats.tsv")) as fh:
    st = [l.rstrip("\n").split("\t") for l in fh]
table = "<table><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in st[0]) + "</tr>" + "".join(
    "<tr>" + "".join(f"<td>{html.escape(v)}</td>" for v in row) + "</tr>" for row in st[1:]) + "</table>"
mark = os.path.basename(D)
with open(os.path.join(Q, "report.html"), "w") as out:
    out.write(f"""<!doctype html><html><head><meta charset="utf-8"><title>{mark} counts QC</title>
<style>body{{font-family:sans-serif;max-width:1100px;margin:24px auto;padding:0 16px}}
table{{border-collapse:collapse;font-size:13px}}td,th{{border:1px solid #ccc;padding:3px 6px;text-align:right}}
img{{max-width:100%}}</style></head><body>
<h1>{mark}: counts QC before differential testing</h1>
<p>{len(ids):,} union regions counted; {int(mask.sum()):,} pass the pre-filter and are used here.
Values: log2(CPM + 1), CPM over all fragments of the library.</p>
<h2>Libraries</h2>{table}
<h2>PCA</h2><img src="data:image/png;base64,{b64(pca_png)}">
<h2>Correlation</h2><img src="data:image/png;base64,{b64(fig2)}">
</body></html>""")
print(f"QC: {int(mask.sum())} filtered regions; PC1 {100 * expl[0]:.1f}%, PC2 {100 * expl[1]:.1f}%, PC3 {100 * expl[2]:.1f}%")
for j, s in enumerate(names):
    print(f"PCA {s}\t" + "\t".join(f"{v:.2f}" for v in pcs[j]))
