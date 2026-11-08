#!/usr/bin/env python3
"""ISL1 replicate correlation and PCA from the featureCounts tables of scripts/09.

Usage: python3 replicate_correlation.py <results/qc/replicate_correlation/isl1>

Reads <dir>/samples.tsv (ISL1 sample, antibody, matched IgG) and <dir>/counts/<set>.counts.txt
for each day and "all". Values are log2(CPM + 1), where CPM uses all fragments featureCounts saw
for that library (assigned + unassigned, from the .summary file). Writes correlation tables, PNG
figures and report.html (figures embedded) to <dir>.
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
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr

D = sys.argv[1]
DAYS = ["D9", "D15", "D18", "D22"]
DAY_COLOR = {"D9": "#4C72B0", "D15": "#DD8452", "D18": "#55A868", "D22": "#C44E52"}
TOP_VAR = 5000


def load(set_name):
    """-> (peak ids, sample names, counts [peaks x samples], total fragments per sample)."""
    path = os.path.join(D, "counts", f"{set_name}.counts.txt")
    with open(path) as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#")]
    head, rows = lines[0], lines[1:]
    names = [os.path.basename(c).replace(".filtered.bam", "") for c in head[6:]]
    counts = np.array([[float(x) for x in r[6:]] for r in rows])
    with open(path + ".summary") as fh:
        summ = [l.rstrip("\n").split("\t") for l in fh][1:]
    total = np.array([sum(float(r[i]) for r in summ) for i in range(1, len(names) + 1)])
    return [r[0] for r in rows], names, counts, total


def log_cpm(counts, total):
    return np.log2(counts / total * 1e6 + 1)


def short(name):
    # D18_H3CTR_rep1_ISL1 -> D18 rep1 ISL1
    p = name.split("_")
    return f"{p[0]} {p[2]} {p[3]}"


def fig_to_b64(fig, png_path):
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def heatmap(mat, labels, title, png, cluster=True, vmin=None, vmax=1.0):
    order = list(range(len(labels)))
    if cluster and len(labels) > 2:
        dist = np.clip(1 - mat, 0, None)
        np.fill_diagonal(dist, 0)
        order = list(leaves_list(linkage(squareform(dist, checks=False), "average")))
    m = mat[np.ix_(order, order)]
    lab = [labels[i] for i in order]
    fig, ax = plt.subplots(figsize=(0.55 * len(lab) + 2.5, 0.55 * len(lab) + 1.5))
    im = ax.imshow(m, cmap="viridis", vmin=vmin if vmin is not None else np.nanmin(m), vmax=vmax)
    ax.set_xticks(range(len(lab)), lab, rotation=90, fontsize=8)
    ax.set_yticks(range(len(lab)), lab, fontsize=8)
    for i in range(len(lab)):
        for j in range(len(lab)):
            ax.text(j, i, f"{m[i, j]:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if m[i, j] < (np.nanmin(m) + vmax) / 2 else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r (log2 CPM)")
    ax.set_title(title, fontsize=10)
    return fig_to_b64(fig, png)


def pca(x):
    """x: samples x features -> (scores, explained variance fraction)."""
    xc = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(xc, full_matrices=False)
    var = s ** 2 / np.sum(s ** 2)
    return u * s, var


def pca_plot(scores, var, names, title, png, pcs=(0, 1)):
    a, b = pcs
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for i, n in enumerate(names):
        day, rep, ab = n.split("_")[0], n.split("_")[2], n.split("_")[3]
        ax.scatter(scores[i, a], scores[i, b], s=70, color=DAY_COLOR[day],
                   marker="o" if ab == "ISL1" else "X", edgecolor="black", linewidth=0.5)
        ax.annotate(f"{day} {rep}" + ("" if ab == "ISL1" else " IgG"), (scores[i, a], scores[i, b]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel(f"PC{a + 1} ({100 * var[a]:.1f}%)")
    ax.set_ylabel(f"PC{b + 1} ({100 * var[b]:.1f}%)")
    ax.set_title(title, fontsize=10)
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=d) for d, c in DAY_COLOR.items()]
    handles += [plt.Line2D([], [], marker="o", ls="", color="grey", label="ISL1"),
                plt.Line2D([], [], marker="X", ls="", color="grey", label="IgG")]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="best")
    return fig_to_b64(fig, png)


figs, tables = [], []

# 1. Per day: union of rep1 + rep2 peaks.
per_day = []
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for ax, day in zip(axes, DAYS):
    _, names, counts, total = load(day)
    lc = log_cpm(counts, total)
    isl = [i for i, n in enumerate(names) if n.startswith(day + "_") and n.endswith("_ISL1")]
    igg = [i for i, n in enumerate(names) if n.startswith(day + "_") and n.endswith("_IgG")]
    r1, r2 = isl
    g1, g2 = igg
    row = {
        "day": day, "peaks": counts.shape[0],
        "rep1": names[r1], "rep2": names[r2],
        "pearson_rep1_rep2": pearsonr(lc[:, r1], lc[:, r2])[0],
        "spearman_rep1_rep2": spearmanr(counts[:, r1], counts[:, r2])[0],
        "pearson_rep1_IgG1": pearsonr(lc[:, r1], lc[:, g1])[0],
        "pearson_rep2_IgG2": pearsonr(lc[:, r2], lc[:, g2])[0],
        "pearson_IgG1_IgG2": pearsonr(lc[:, g1], lc[:, g2])[0],
        "frac_in_peaks_rep1": counts[:, r1].sum() / total[r1],
        "frac_in_peaks_rep2": counts[:, r2].sum() / total[r2],
    }
    per_day.append(row)
    ax.scatter(lc[:, r1], lc[:, r2], s=2, alpha=0.25, color=DAY_COLOR[day], rasterized=True)
    lim = max(lc[:, r1].max(), lc[:, r2].max())
    ax.plot([0, lim], [0, lim], color="grey", lw=0.8, ls="--")
    ax.set_xlabel(f"{short(names[r1])} log2(CPM+1)")
    ax.set_ylabel(f"{short(names[r2])} log2(CPM+1)")
    ax.set_title(f"{day}: {counts.shape[0]:,} union peaks\nPearson {row['pearson_rep1_rep2']:.2f}, "
                 f"Spearman {row['spearman_rep1_rep2']:.2f}", fontsize=9)
figs.append(("Replicate scatter per day (union of the day's two replicates' peaks)",
             fig_to_b64(fig, os.path.join(D, "per_day_scatter.png"))))

cols = ["day", "peaks", "rep1", "rep2", "pearson_rep1_rep2", "spearman_rep1_rep2", "pearson_rep1_IgG1",
        "pearson_rep2_IgG2", "pearson_IgG1_IgG2", "frac_in_peaks_rep1", "frac_in_peaks_rep2"]
with open(os.path.join(D, "per_day_correlation.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in per_day:
        fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
tables.append(("Per-day correlation (per-day union peaks)", cols, per_day))

# Per-day summary heatmap: days x comparisons.
comp = ["pearson_rep1_rep2", "spearman_rep1_rep2", "pearson_rep1_IgG1", "pearson_rep2_IgG2", "pearson_IgG1_IgG2"]
m = np.array([[r[c] for c in comp] for r in per_day])
fig, ax = plt.subplots(figsize=(7, 3.2))
im = ax.imshow(m, cmap="viridis", vmin=min(0, m.min()), vmax=1)
ax.set_xticks(range(len(comp)), ["rep1 vs rep2\nPearson", "rep1 vs rep2\nSpearman", "rep1 vs IgG1\nPearson",
                                 "rep2 vs IgG2\nPearson", "IgG1 vs IgG2\nPearson"], fontsize=8)
ax.set_yticks(range(len(DAYS)), DAYS)
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        ax.text(j, i, f"{m[i, j]:.2f}", ha="center", va="center", fontsize=8,
                color="white" if m[i, j] < 0.5 else "black")
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.03)
ax.set_title("Replicate correlation by day (per-day union peaks)", fontsize=10)
figs.append(("Replicate correlation by day", fig_to_b64(fig, os.path.join(D, "per_day_heatmap.png"))))

# 2. All samples: union of all 8 ISL1 replicates' peaks.
_, names, counts, total = load("all")
lc = log_cpm(counts, total)
isl = [i for i, n in enumerate(names) if n.endswith("_ISL1")]
labels = [short(n) for n in names]
corr_all = np.corrcoef(lc.T)
with open(os.path.join(D, "all_samples_pearson.tsv"), "w") as fh:
    fh.write("\t" + "\t".join(names) + "\n")
    for n, row in zip(names, corr_all):
        fh.write(n + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")
n_all = counts.shape[0]
figs.append((f"ISL1 libraries, all-sample union ({n_all:,} peaks)",
             heatmap(corr_all[np.ix_(isl, isl)], [labels[i] for i in isl],
                     f"ISL1, Pearson log2 CPM, {n_all:,} peaks", os.path.join(D, "heatmap_isl1.png"))))
figs.append((f"ISL1 and IgG libraries, all-sample union ({n_all:,} peaks)",
             heatmap(corr_all, labels, f"ISL1 + IgG, Pearson log2 CPM, {n_all:,} peaks",
                     os.path.join(D, "heatmap_isl1_igg.png"))))

# PCA: ISL1 only and ISL1 + IgG; all peaks and the most variable ones.
for tag, idx in (("ISL1", isl), ("ISL1 + IgG", list(range(len(names))))):
    x = lc[:, idx].T
    for sub, feats in (("all peaks", np.arange(x.shape[1])),
                       (f"top {TOP_VAR:,} variable", np.argsort(x.var(axis=0))[::-1][:TOP_VAR])):
        scores, var = pca(x[:, feats])
        stem = f"pca_{tag.replace(' + ', '_').lower()}_{'all' if sub == 'all peaks' else 'top'}"
        figs.append((f"PCA, {tag}, {sub}",
                     pca_plot(scores, var, [names[i] for i in idx], f"PCA {tag}, {sub} (PC1 vs PC2)",
                              os.path.join(D, stem + ".png"))))
        if tag == "ISL1":
            figs.append((f"PCA, {tag}, {sub}, PC2 vs PC3",
                         pca_plot(scores, var, [names[i] for i in idx], f"PCA {tag}, {sub} (PC2 vs PC3)",
                                  os.path.join(D, stem + "_pc23.png"), pcs=(1, 2))))

# Report.
css = """
:root { --bg: #ffffff; --fg: #1f2328; --muted: #59636e; --line: #d1d9e0; }
@media (prefers-color-scheme: dark) { :root { --bg: #0d1117; --fg: #e6edf3; --muted: #9198a1; --line: #3d444d; } }
body { background: var(--bg); color: var(--fg); font: 14px/1.5 -apple-system, Helvetica, Arial, sans-serif;
       max-width: 1100px; margin: 0 auto; padding: 16px; }
h1 { font-size: 22px; } h2 { font-size: 17px; margin-top: 32px; border-bottom: 1px solid var(--line); }
p.note { color: var(--muted); }
table { border-collapse: collapse; font-size: 12px; overflow-x: auto; display: block; }
td, th { border: 1px solid var(--line); padding: 4px 8px; text-align: right; }
th:first-child, td:first-child { text-align: left; }
img { max-width: 100%; background: #fff; border-radius: 4px; }
"""
out = [f"<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>"
       f"<title>ISL1 Replicate Correlation</title><style>{css}</style></head><body>",
       "<h1>ISL1 replicate correlation and PCA (H3CTR, matched IgG as background)</h1>",
       "<p class='note'>Peaks: per-replicate MACS3 calls (q 0.01, blacklist filtered, chr1-22/X/Y). "
       "Counts: fragments per peak (featureCounts -p --countReadPairs, duplicates kept). "
       "Values: log2(CPM + 1), CPM over all fragments in the library. Pearson on log2 CPM, Spearman on counts. "
       f"PCA centered, not scaled; 'top variable' = the {TOP_VAR:,} peaks with the highest variance across the "
       "samples in that PCA. Heatmaps are ordered by average-linkage clustering on 1 - r.</p>"]
for title, cols_, rows in tables:
    out.append(f"<h2>{html.escape(title)}</h2><table><tr>" + "".join(f"<th>{c}</th>" for c in cols_) + "</tr>")
    for r in rows:
        out.append("<tr>" + "".join(f"<td>{r[c]:.3f}</td>" if isinstance(r[c], float) else
                                     f"<td>{html.escape(str(r[c]))}</td>" for c in cols_) + "</tr>")
    out.append("</table>")
for title, b64 in figs:
    out.append(f"<h2>{html.escape(title)}</h2><img alt='{html.escape(title)}' src='data:image/png;base64,{b64}'>")
out.append("</body></html>")
with open(os.path.join(D, "report.html"), "w") as fh:
    fh.write("\n".join(out))
for r in per_day:
    print("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols))
print("written", os.path.join(D, "report.html"))
