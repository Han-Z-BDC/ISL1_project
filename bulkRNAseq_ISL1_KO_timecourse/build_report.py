#!/usr/bin/env python3
"""Build the DESeq2 time-course report for the ISL1 KO bulk RNA-seq (YJZ5LS).

Run on the Mac after fetching from Alpine:
  results/DE_DESeq2/*            (pipeline/06_deseq2_timecourse.R)
  results/qc/multiqc_*.txt, results/qc/multiqc_full_report.html   (pipeline/04b_multiqc_full.sb)

Usage:  ~/anaconda3/bin/python3 build_report.py
Output: results/report_DESeq2_timecourse.html
        results/DE_DESeq2/timecourse_clusters_LRT_interaction.tsv

Time-course clusters: genes with LRT_interaction padj < 0.05 (the KO effect changes with day),
z-scored mean VST per group (Ctr/KO x D15/D18/D22), k-means with N_CLUSTERS clusters
(fixed seed, best of N_INIT starts; clusters numbered by size).
Uses plotly only (matplotlib does not import in this env).
"""
import html
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
DE = os.path.join(RES, "DE_DESeq2")
QC = os.path.join(RES, "qc")
OUT = os.path.join(RES, "report_DESeq2_timecourse.html")

ALPHA = 0.05
N_CLUSTERS = 6
N_INIT = 25
SEED = 12345
N_VOLCANO_LABELS = 6
MARKERS = ["ISL1", "ARX", "GCG", "SST", "INS", "MAFA", "PDX1", "NKX6-1", "NEUROG3", "CHGA"]

DAYS = ["D15", "D18", "D22"]
GENOS = ["Ctr", "KO"]
GROUPS = [f"{g}_{d}" for g in GENOS for d in DAYS]

# Colors (reference palette, light mode). Genotype identity: slot 1 / slot 2.
C = {"Ctr": "#2a78d6", "KO": "#eb6834", "ns": "#c9c8c2", "text": "#0b0b0b",
     "text2": "#52514e", "grid": "#e6e5e0", "surface": "#fcfcfb"}
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DAY_SYMBOL = {"D15": "circle", "D18": "square", "D22": "diamond"}


def read(name, folder=DE):
    return pd.read_csv(os.path.join(folder, name), sep="\t", keep_default_na=False, na_values=["NA"])


def layout(fig, height, **kw):
    fig.update_layout(
        height=height, margin=dict(l=60, r=20, t=40, b=50),
        paper_bgcolor=C["surface"], plot_bgcolor=C["surface"],
        font=dict(family="Inter, -apple-system, Segoe UI, Helvetica, Arial, sans-serif", size=12, color=C["text"]),
        hoverlabel=dict(bgcolor="white", font_size=12), **kw)
    fig.update_xaxes(gridcolor=C["grid"], zerolinecolor=C["grid"], linecolor=C["text2"], ticks="outside",
                     tickcolor=C["grid"], title_font=dict(color=C["text2"]))
    fig.update_yaxes(gridcolor=C["grid"], zerolinecolor=C["grid"], linecolor=C["text2"], ticks="outside",
                     tickcolor=C["grid"], title_font=dict(color=C["text2"]))
    return fig


_first_fig = [True]


def fig_html(fig):
    inc = _first_fig[0]
    _first_fig[0] = False
    return fig.to_html(full_html=False, include_plotlyjs=inc, config={"displaylogo": False, "responsive": True})


def table_html(df, cls="tbl", fmt=None):
    fmt = fmt or {}
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            s = fmt[c](v) if c in fmt else ("" if pd.isna(v) else str(v))
            num = isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool)
            cells.append(f'<td class="{"num" if num else ""}">{html.escape(s)}</td>')
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f'<div class="tblwrap"><table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def short(sample):
    p = sample.split("_")
    g = "Ctr" if p[2].startswith("Ctr") else "KO"
    return f"{g} {p[3]} r{p[4][-1]}"


# ---------------------------------------------------------------- inputs
meta = pd.read_csv(os.path.join(HERE, "sample_sheet.csv"))
meta["group"] = meta["genotype"] + "_" + meta["day"]
meta["label"] = meta["sample"].map(short)
order = meta.sort_values(["genotype", "day", "replicate"], kind="stable")["sample"].tolist()

summary = read("summary.tsv")
pca = read("pca_top500.tsv")
pvar = read("pca_variance.tsv")
cor = read("sample_correlation_vst.tsv").set_index("sample")
vst = read("vst.tsv")
norm = read("normalized_counts.tsv")
contr = {d: read(f"KO_vs_Ctr_{d}.tsv") for d in DAYS}
lrt_int = read("LRT_interaction.tsv")
lrt_geno = read("LRT_genotype.tsv")

# ---------------------------------------------------------------- QC table
star = read("multiqc_star.txt", QC).set_index("Sample")
dd = read("multiqc_umitools_dedup.txt", QC).set_index("Sample")
fc = read("multiqc_featurecounts.txt", QC)
fc["Sample"] = fc["Sample"].str.replace(r"\.dedup$", "", regex=True)
fc = fc.set_index("Sample")
ie = read("multiqc_rseqc_infer_experiment.txt", QC)
ie["Sample"] = ie["Sample"].str.replace(r"\.infer_experiment$", "", regex=True)
ie = ie.set_index("Sample")
qc = pd.DataFrame({
    "Sample": [short(s) for s in order],
    "Reads (M)": [star.loc[s, "total_reads"] / 1e6 for s in order],
    "Unique map %": [star.loc[s, "uniquely_mapped_percent"] for s in order],
    "Multi map %": [star.loc[s, "multimapped_percent"] for s in order],
    "UMI duplicates %": [100 - dd.loc[s, "dedup_percent_passing"] for s in order],
    "Assigned to genes (M)": [fc.loc[s, "Assigned"] / 1e6 for s in order],
    "Sense strand %": [100 * ie.loc[s, "se_sense"] for s in order],
})
f1 = lambda v: f"{v:.1f}"
qc_tbl = table_html(qc, fmt={c: f1 for c in qc.columns[1:]})

# ---------------------------------------------------------------- DE summary table
lab = {f"KO_vs_Ctr_{d}": f"KO vs Ctr, {d}" for d in DAYS}
lab.update({f"{g}_{b}_vs_{a}": f"{g}: {b} vs {a}" for g in GENOS for a, b in [("D15", "D18"), ("D15", "D22"), ("D18", "D22")]})
lab.update({"LRT_interaction": "KO effect changes with day (LRT)", "LRT_genotype": "KO effect at any day (LRT)"})
s = summary.copy()
s.insert(0, "Test", s["test"].map(lab))
s = s.drop(columns="test").rename(columns={"n_tested": "Genes tested", "n_padj_lt_0.05": "padj < 0.05",
                                           "n_up": "Up", "n_down": "Down"})
intfmt = lambda v: "" if pd.isna(v) else f"{int(v):,}"
sum_tbl = table_html(s, fmt={c: intfmt for c in s.columns[1:]})

# ---------------------------------------------------------------- PCA
pv = dict(zip(pvar["PC"], pvar["pct_variance"]))
fig = go.Figure()
for g in GENOS:
    for d in DAYS:
        sub = pca[(pca["genotype"] == g) & (pca["day"] == d)]
        fig.add_trace(go.Scatter(
            x=sub["PC1"], y=sub["PC2"], mode="markers", name=f"{g} {d}",
            legendgroup=g, legendgrouptitle_text=g if d == "D15" else None,
            marker=dict(color=C[g], symbol=DAY_SYMBOL[d], size=13, line=dict(color=C["surface"], width=2)),
            customdata=sub["sample"], hovertemplate="%{customdata}<br>PC1 %{x:.1f}, PC2 %{y:.1f}<extra></extra>"))
        fig.add_annotation(x=sub["PC1"].mean(), y=sub["PC2"].max(), text=f"{g} {d}", showarrow=False, yshift=16,
                           font=dict(size=11, color=C["text2"]))
fig.update_xaxes(title=f"PC1 ({pv['PC1']:.1f}% of variance)")
fig.update_yaxes(title=f"PC2 ({pv['PC2']:.1f}% of variance)")
pca_html = fig_html(layout(fig, 460, legend=dict(groupclick="toggleitem")))

# ---------------------------------------------------------------- correlation heatmap
cm = cor.loc[order, order]
labels = [short(x) for x in order]
fig = go.Figure(go.Heatmap(
    z=cm.values, x=labels, y=labels, colorscale=[[i / (len(SEQ_BLUE) - 1), c] for i, c in enumerate(SEQ_BLUE)],
    zmin=float(np.floor(cm.values.min() * 100) / 100), zmax=1, xgap=2, ygap=2,
    colorbar=dict(title=dict(text="Pearson r", font=dict(color=C["text2"])), thickness=12),
    hovertemplate="%{y} vs %{x}<br>r = %{z:.4f}<extra></extra>"))
fig.update_yaxes(autorange="reversed", showgrid=False)
fig.update_xaxes(showgrid=False, tickangle=-45)
cor_html = fig_html(layout(fig, 520))

# ---------------------------------------------------------------- volcanoes
fig = make_subplots(rows=1, cols=3, shared_yaxes=True, horizontal_spacing=0.03,
                    subplot_titles=[f"KO vs Ctr, {d}" for d in DAYS])
YCAP = 300
xmax = max(np.nanmax(np.abs(contr[d]["log2FoldChange"].to_numpy(dtype=float))) for d in DAYS)
xlim = float(min(np.ceil(xmax), 14))
for j, d in enumerate(DAYS, start=1):
    df = contr[d].dropna(subset=["padj"]).copy()
    df["y"] = np.minimum(-np.log10(df["padj"].clip(lower=1e-300)), YCAP)
    df["x"] = df["log2FoldChange"].clip(-xlim, xlim)
    sig = df["padj"] < ALPHA
    cats = [("ns", ~sig, "Not significant"), ("KO", sig & (df["log2FoldChange"] > 0), "Higher in KO"),
            ("Ctr", sig & (df["log2FoldChange"] < 0), "Higher in Ctr")]
    for key, m, nm in cats:
        sub = df[m]
        fig.add_trace(go.Scattergl(
            x=sub["x"], y=sub["y"], mode="markers", name=f"{nm}", legendgroup=key, showlegend=(j == 1),
            marker=dict(color=C[key], size=5 if key == "ns" else 6, opacity=0.6 if key == "ns" else 0.85),
            customdata=np.stack([sub["gene_name"], sub["log2FoldChange"], sub["padj"], sub["baseMean"]], axis=1),
            hovertemplate="<b>%{customdata[0]}</b><br>log2FC %{customdata[1]:.2f}<br>padj %{customdata[2]:.2e}"
                          "<br>mean count %{customdata[3]:.0f}<extra></extra>"), row=1, col=j)
    top = df[sig].sort_values(["padj", "pvalue", "gene_id"], kind="stable").head(N_VOLCANO_LABELS)
    lab_df = pd.concat([top, df[df["gene_name"] == "ISL1"]]).drop_duplicates("gene_id")
    for i, (_, r) in enumerate(lab_df.iterrows()):
        # staggered offsets so neighbouring labels do not sit on each other
        side = 1 if r["x"] > 0 else -1
        fig.add_annotation(x=r["x"], y=r["y"], text=r["gene_name"], showarrow=True, arrowhead=0, arrowwidth=1,
                           arrowcolor=C["text2"], ax=side * (30 + 25 * (i % 2)), ay=[-10, 14, -30, 30][i % 4],
                           font=dict(size=10, color=C["text"]), bgcolor="rgba(252,252,251,0.85)", row=1, col=j)
    fig.add_hline(y=-np.log10(ALPHA), line=dict(color=C["text2"], width=1, dash="dot"), row=1, col=j)
    fig.update_xaxes(title="log2 fold change (KO / Ctr)", range=[-xlim - 0.5, xlim + 0.5], row=1, col=j)
fig.update_yaxes(title="-log10 adjusted p", row=1, col=1)
for a in fig.layout.annotations[:3]:
    a.font = dict(size=13, color=C["text"])
volc_html = fig_html(layout(fig, 520, legend=dict(orientation="h", y=-0.18)))


# ---------------------------------------------------------------- time-course clusters
def kmeans(X, k, seed, n_init, iters=300):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_init):
        # k-means++ initialisation
        cent = [X[rng.integers(len(X))]]
        for _ in range(1, k):
            d2 = np.min(((X[:, None, :] - np.array(cent)[None]) ** 2).sum(-1), axis=1)
            cent.append(X[rng.choice(len(X), p=d2 / d2.sum())])
        cent = np.array(cent)
        lab_ = None
        for _ in range(iters):
            dist = ((X[:, None, :] - cent[None]) ** 2).sum(-1)
            new = np.argmin(dist, axis=1)
            if lab_ is not None and np.array_equal(new, lab_):
                break
            lab_ = new
            for c in range(k):
                if np.any(lab_ == c):
                    cent[c] = X[lab_ == c].mean(0)
        inertia = ((X - cent[lab_]) ** 2).sum()
        if best is None or inertia < best[0]:
            best = (inertia, lab_.copy())
    return best[1]


sig_int = lrt_int[lrt_int["padj"] < ALPHA][["gene_id", "gene_name", "padj"]]
v = vst.set_index("gene_id").loc[sig_int["gene_id"]]
gm = pd.DataFrame({g: v[meta.loc[meta["group"] == g, "sample"]].mean(axis=1) for g in GROUPS})
z = gm.sub(gm.mean(axis=1), axis=0).div(gm.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
raw = kmeans(z.to_numpy(), N_CLUSTERS, SEED, N_INIT)
sizes = pd.Series(raw).value_counts()
first = pd.Series(np.arange(len(raw))).groupby(raw).min()
rank = sorted(sizes.index, key=lambda c: (-sizes[c], first[c]))
remap = {c: i + 1 for i, c in enumerate(rank)}
clus = np.array([remap[c] for c in raw])

ct = sig_int.reset_index(drop=True).copy()
ct.insert(2, "cluster", clus)
for g in GROUPS:
    ct[f"meanVST_{g}"] = gm[g].round(4).to_numpy()
for g in GROUPS:
    ct[f"z_{g}"] = z[g].round(4).to_numpy()
ct = ct.sort_values(["cluster", "padj", "gene_id"], kind="stable")
ct.to_csv(os.path.join(DE, "timecourse_clusters_LRT_interaction.tsv"), sep="\t", index=False)

ncol = 3
nrow = int(np.ceil(N_CLUSTERS / ncol))
titles = [f"Cluster {c} (n = {int((clus == c).sum()):,})" for c in range(1, N_CLUSTERS + 1)]
fig = make_subplots(rows=nrow, cols=ncol, shared_yaxes=True, shared_xaxes=True, subplot_titles=titles,
                    vertical_spacing=0.12, horizontal_spacing=0.04)
zz = z.to_numpy()
for c in range(1, N_CLUSTERS + 1):
    r_, c_ = (c - 1) // ncol + 1, (c - 1) % ncol + 1
    m = clus == c
    for g in GENOS:
        cols = [GROUPS.index(f"{g}_{d}") for d in DAYS]
        q1, med, q3 = (np.percentile(zz[m][:, cols], p, axis=0) for p in (25, 50, 75))
        fig.add_trace(go.Scatter(x=DAYS + DAYS[::-1], y=list(q3) + list(q1[::-1]), fill="toself",
                                 fillcolor=C[g], opacity=0.15, line=dict(width=0), hoverinfo="skip",
                                 showlegend=False, legendgroup=g), row=r_, col=c_)
        fig.add_trace(go.Scatter(x=DAYS, y=med, mode="lines+markers", name=g, legendgroup=g,
                                 showlegend=(c == 1), line=dict(color=C[g], width=2),
                                 marker=dict(size=8, color=C[g], line=dict(color=C["surface"], width=2)),
                                 hovertemplate=f"{g} %{{x}}<br>median z %{{y:.2f}}<extra>Cluster {c}</extra>"),
                      row=r_, col=c_)
for r_ in range(1, nrow + 1):
    fig.update_yaxes(title="z-score (mean VST)", row=r_, col=1)
for a in fig.layout.annotations:
    a.font = dict(size=13, color=C["text"])
clus_html = fig_html(layout(fig, 300 * nrow + 60, legend=dict(orientation="h", y=1.08)))

clus_tops = []
for c in range(1, N_CLUSTERS + 1):
    t = ct[ct["cluster"] == c].head(25)["gene_name"].tolist()
    clus_tops.append(f"<li><b>Cluster {c}</b> ({int((clus == c).sum()):,} genes): {html.escape(', '.join(t))}</li>")

# ---------------------------------------------------------------- heatmap of the clustered genes
# Rows: the clustered genes, grouped by cluster; within a cluster ordered by average-linkage
# hierarchical clustering (Euclidean, optimal leaf order) on the same group-mean z-scores used
# for k-means. Columns: the 12 samples. Colour: VST z-scored per gene across the 12 samples.
from scipy.cluster.hierarchy import leaves_list, linkage

HEAT_LABELS = ["INS", "ISL1", "MAFA", "MAFB", "PAX6", "NEUROD1", "NKX2-2", "GCG", "ARX", "SST", "IRX2", "TTR",
               "LOXL4", "ENTPD3", "FFAR1", "HADH", "NKX6-1", "TAC1", "IFI6", "CXCL1", "NEUROG3", "TPH1"]
ZCLIP = 2.5
gid_all = sig_int["gene_id"].to_numpy()
row_order = []
for c in range(1, N_CLUSTERS + 1):
    idx = np.where(clus == c)[0]
    if len(idx) > 2:
        lk = linkage(z.to_numpy()[idx], method="average", metric="euclidean", optimal_ordering=True)
        idx = idx[leaves_list(lk)]
    row_order.extend(idx.tolist())
row_order = np.array(row_order)
hv = v.loc[gid_all[row_order], order].astype(float)
hz = hv.sub(hv.mean(axis=1), axis=0).div(hv.std(axis=1, ddof=1).replace(0, np.nan), axis=0).fillna(0)
hnames = sig_int["gene_name"].to_numpy()[row_order]
hclus = clus[row_order]
n_rows = len(row_order)

DIVERGING = [[0.0, "#184f95"], [0.2, "#3987e5"], [0.4, "#b7d3f6"], [0.5, "#f0efec"],
             [0.6, "#f4b9b8"], [0.8, "#e34948"], [1.0, "#a3292a"]]
xlab = [short(s_) for s_ in order]
fig = go.Figure(go.Heatmap(
    z=hz.clip(-ZCLIP, ZCLIP).to_numpy(), x=xlab, y=np.arange(n_rows), colorscale=DIVERGING, zmin=-ZCLIP, zmax=ZCLIP,
    customdata=np.repeat(np.stack([hnames, hclus.astype(str)], axis=1)[:, None, :], len(order), axis=1),
    hovertemplate="<b>%{customdata[0]}</b> (cluster %{customdata[1]})<br>%{x}<br>z = %{z:.2f}<extra></extra>",
    colorbar=dict(title=dict(text="z-score<br>(VST)", font=dict(color=C["text2"])), thickness=12, len=0.5, x=1.16)))
# cluster bands: separators and labels
bounds = np.r_[0, np.cumsum([int((clus == c).sum()) for c in range(1, N_CLUSTERS + 1)])]
for b in bounds[1:-1]:
    fig.add_shape(type="line", xref="paper", x0=0, x1=1, yref="y", y0=b - 0.5, y1=b - 0.5,
                  line=dict(color=C["surface"], width=4))
fig.add_shape(type="line", xref="x", x0=5.5, x1=5.5, yref="paper", y0=0, y1=1, line=dict(color=C["surface"], width=4))
centers = (bounds[:-1] + bounds[1:]) / 2 - 0.5
fig.update_yaxes(autorange="reversed", showgrid=False, tickmode="array", tickvals=centers,
                 ticktext=[f"Cluster {c}<br>n = {int((clus == c).sum())}" for c in range(1, N_CLUSTERS + 1)],
                 ticks="", zeroline=False)
fig.update_xaxes(showgrid=False, tickangle=-45, side="bottom")
# gene labels on the right, spread out so they do not overlap, with leader lines to their rows
HEAT_H, MT, MB = 1000, 40, 90
px_per_row = (HEAT_H - MT - MB) / n_rows
pos = {g: int(np.where(hnames == g)[0][0]) for g in HEAT_LABELS if g in set(hnames)}
items = sorted(pos.items(), key=lambda kv: kv[1])
min_gap = 16 / px_per_row  # 16 px between labels
placed = []
for g, r in items:
    y_lab = r if not placed else max(r, placed[-1][2] + min_gap)
    placed.append((g, r, y_lab))
overflow = placed[-1][2] - (n_rows - 1) if placed else 0
if overflow > 0:  # shift the stack up if it runs past the bottom
    placed = [(g, r, yl - overflow) for g, r, yl in placed]
    for i in range(len(placed) - 2, -1, -1):
        if placed[i][2] > placed[i + 1][2] - min_gap:
            placed[i] = (placed[i][0], placed[i][1], placed[i + 1][2] - min_gap)
for g, r, yl in placed:
    fig.add_annotation(xref="paper", x=1.0, yref="y", y=r, text=g, showarrow=True, arrowhead=0, arrowwidth=1,
                       arrowcolor=C["text2"], ax=55, ay=(yl - r) * px_per_row, xanchor="left",
                       font=dict(size=10, color=C["text"]))
layout(fig, HEAT_H)
fig.update_layout(margin=dict(l=90, r=170, t=MT, b=MB))
heat_html = fig_html(fig)

# ---------------------------------------------------------------- marker genes
nm = norm.set_index("gene_name")
mk = [g for g in MARKERS if g in nm.index]
ncol = 5
nrow = int(np.ceil(len(mk) / ncol))
fig = make_subplots(rows=nrow, cols=ncol, subplot_titles=mk, vertical_spacing=0.16, horizontal_spacing=0.05)
for i, gene in enumerate(mk):
    r_, c_ = i // ncol + 1, i % ncol + 1
    row = nm.loc[gene]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    for g in GENOS:
        xs, ys, hs = [], [], []
        for d in DAYS:
            ss = meta[(meta["genotype"] == g) & (meta["day"] == d)].sort_values("replicate", kind="stable")["sample"]
            for smp in ss:
                xs.append(d); ys.append(float(row[smp]) + 1); hs.append(short(smp))
        means = [np.mean([ys[k] for k in range(len(xs)) if xs[k] == d]) for d in DAYS]
        fig.add_trace(go.Scatter(x=DAYS, y=means, mode="lines", line=dict(color=C[g], width=2), name=g,
                                 legendgroup=g, showlegend=(i == 0), hoverinfo="skip"), row=r_, col=c_)
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", marker=dict(color=C[g], size=8, line=dict(color=C["surface"], width=2)),
                                 legendgroup=g, showlegend=False, customdata=hs,
                                 hovertemplate=f"{gene}<br>%{{customdata}}<br>%{{y:.0f}} (normalized count + 1)<extra></extra>"),
                      row=r_, col=c_)
    fig.update_yaxes(type="log", row=r_, col=c_)
for r_ in range(1, nrow + 1):
    fig.update_yaxes(title="normalized count + 1", row=r_, col=1)
for a in fig.layout.annotations:
    a.font = dict(size=13, color=C["text"])
mark_html = fig_html(layout(fig, 280 * nrow + 60, legend=dict(orientation="h", y=1.1)))

marker_rows = []
for gene in mk:
    rr = {"Gene": gene}
    for d in DAYS:
        hit = contr[d][contr[d]["gene_name"] == gene]
        if len(hit):
            h = hit.iloc[0]
            star_ = "" if pd.isna(h["padj"]) or h["padj"] >= ALPHA else " *"
            rr[f"log2FC {d}"] = f"{h['log2FoldChange']:.2f}{star_}"
            rr[f"padj {d}"] = "NA" if pd.isna(h["padj"]) else f"{h['padj']:.1e}"
    marker_rows.append(rr)
marker_tbl = table_html(pd.DataFrame(marker_rows))

# ---------------------------------------------------------------- top genes per day
top_blocks = []
for d in DAYS:
    df = contr[d].dropna(subset=["padj"])
    df = df[df["padj"] < ALPHA]
    for direction, m, title in [("up", df["log2FoldChange"] > 0, "higher in KO"), ("down", df["log2FoldChange"] < 0, "higher in Ctr")]:
        t = df[m].sort_values(["padj", "pvalue", "gene_id"], kind="stable").head(20)
        t = pd.DataFrame({"Gene": t["gene_name"], "log2FC": t["log2FoldChange"].map(lambda x: f"{x:.2f}"),
                          "shrunk": t["log2FC_shrunk"].map(lambda x: f"{x:.2f}"),
                          "padj": t["padj"].map(lambda x: f"{x:.1e}"), "mean count": t["baseMean"].map(lambda x: f"{x:,.0f}")})
        top_blocks.append(f'<div class="half"><h4>{d}: {title}</h4>{table_html(t, cls="tbl small")}</div>')

# ---------------------------------------------------------------- page
n_int = int((lrt_int["padj"] < ALPHA).sum())
n_geno = int((lrt_geno["padj"] < ALPHA).sum())
de_n = {d: int((contr[d]["padj"] < ALPHA).sum()) for d in DAYS}
isl1 = nm.loc["ISL1"]
isl1_ctr = [float(isl1[s]) for s in meta.loc[meta["genotype"] == "Ctr", "sample"]]
isl1_ko = [float(isl1[s]) for s in meta.loc[meta["genotype"] == "KO", "sample"]]

page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>ISL1 KO Time Course</title>
<style>
:root {{ color-scheme: light; --surface:{C['surface']}; --text:{C['text']}; --text2:{C['text2']}; --line:{C['grid']}; }}
body {{ margin:0; background:var(--surface); color:var(--text);
       font:15px/1.55 Inter, -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }}
main {{ max-width:1180px; margin:0 auto; padding:24px 16px 64px; }}
h1 {{ font-size:26px; margin:0 0 4px; }} h2 {{ font-size:20px; margin:40px 0 8px; border-top:1px solid var(--line); padding-top:24px; }}
h4 {{ margin:12px 0 6px; font-size:14px; }}
p, li {{ color:var(--text2); max-width:900px; }} .lead {{ font-size:16px; }}
code {{ font-size:13px; background:#f0efec; padding:1px 4px; border-radius:3px; }}
.tblwrap {{ overflow-x:auto; }}
table.tbl {{ border-collapse:collapse; font-size:13px; margin:8px 0; }}
.tbl th, .tbl td {{ padding:4px 10px; border-bottom:1px solid var(--line); text-align:left; white-space:nowrap; }}
.tbl th {{ color:var(--text2); font-weight:600; }} .tbl td.num {{ text-align:right; font-variant-numeric:tabular-nums; }}
.tbl.small {{ font-size:12px; }}
.grid {{ display:flex; flex-wrap:wrap; gap:0 24px; }} .half {{ flex:1 1 340px; min-width:0; }}
.tiles {{ display:flex; flex-wrap:wrap; gap:12px; margin:16px 0; }}
.tile {{ flex:1 1 160px; border:1px solid var(--line); border-radius:8px; padding:12px 14px; }}
.tile .v {{ font-size:26px; font-weight:600; font-variant-numeric:tabular-nums; }} .tile .l {{ color:var(--text2); font-size:13px; }}
</style></head><body><main>
<h1>ISL1 KO time course: bulk RNA-seq</h1>
<p class="lead">Run YJZ5LS (Plasmidsaurus 3'-end, UMI). Ctr (Ctr-H3) vs ISL1 KO at D15, D18 and D22, 2 replicates each
(D9 samples not used). Reads aligned with STAR, UMI-deduplicated with umi_tools, counted per gene with featureCounts
(forward strand) and tested with DESeq2. Significance: adjusted p &lt; {ALPHA}.</p>

<div class="tiles">
<div class="tile"><div class="v">{de_n['D15']:,}</div><div class="l">genes changed in KO at D15</div></div>
<div class="tile"><div class="v">{de_n['D18']:,}</div><div class="l">genes changed in KO at D18</div></div>
<div class="tile"><div class="v">{de_n['D22']:,}</div><div class="l">genes changed in KO at D22</div></div>
<div class="tile"><div class="v">{n_int:,}</div><div class="l">genes whose KO effect changes with day</div></div>
</div>

<h2>Design and tests</h2>
<p>Model <code>~ genotype + day + genotype:day</code> (fitted as 6 groups), all 12 samples in one fit; day is a factor.
No replicate term (rep1/rep2 are not differentiation batches). Genes kept: &ge; 10 counts in &ge; 2 samples.
KO vs Ctr and day comparisons are Wald tests; adjusted p comes from the unshrunk fold change, and the apeglm-shrunk fold
change is given for ranking. The two likelihood-ratio tests compare the full model with <code>~ genotype + day</code>
(KO effect changes with day) and with <code>~ day</code> (KO effect at any day). Tables: <code>results/DE_DESeq2/</code>.</p>
{sum_tbl}

<h2>Sample QC</h2>
<p>Per-sample read processing. UMI duplicates are reads removed by umi_tools as PCR copies; assigned reads are deduplicated
reads counted on a gene. Sense strand % is from rseqc infer_experiment. Full report:
<a href="qc/multiqc_full_report.html">MultiQC</a>.</p>
{qc_tbl}

<h2>PCA</h2>
<p>VST counts, 500 most variable genes. PC1 ({pv['PC1']:.1f}%) orders the samples by day; PC2 ({pv['PC2']:.1f}%) separates
the genotypes. Replicates sit together.</p>
{pca_html}

<h2>Sample correlation</h2>
<p>Pearson correlation of VST counts across all kept genes, samples ordered by genotype and day.</p>
{cor_html}

<h2>KO vs Ctr at each day</h2>
<p>Each point is a gene. Orange: higher in KO; blue: higher in Ctr; gray: not significant. Labeled: the {N_VOLCANO_LABELS}
genes with the smallest adjusted p, and ISL1. Fold changes beyond &plusmn;{xlim:g} and -log10 p beyond {YCAP} are drawn
at the edge. Dotted line: adjusted p = {ALPHA}.</p>
{volc_html}
<div class="grid">{''.join(top_blocks)}</div>

<h2>Time-course patterns of the interaction genes</h2>
<p>The {n_int:,} genes whose KO effect changes with day (interaction LRT), grouped by their pattern across the six groups:
mean VST per group, z-scored per gene, k-means with {N_CLUSTERS} clusters. Lines: median z-score; bands: 25th-75th
percentile. Gene lists: <code>results/DE_DESeq2/timecourse_clusters_LRT_interaction.tsv</code>.
The cluster number is a choice, not a result; the patterns are the summary.</p>
{clus_html}
<p>Most significant genes per cluster (by interaction p):</p>
<ul>{''.join(clus_tops)}</ul>
<h4>Heatmap of the clustered genes</h4>
<p>All {n_int:,} clustered genes (rows) in the 12 samples (columns, Ctr then KO, by day). Colour: VST z-scored per gene
across the 12 samples, clipped at &plusmn;{ZCLIP:g}; red above the gene's mean, blue below. Rows are grouped by cluster and
ordered within each cluster by hierarchical clustering (average linkage) of the group-mean z-scores used for k-means.
Labeled on the right: selected islet and marker genes.</p>
{heat_html}

<h2>Islet marker genes</h2>
<p>Normalized counts (log scale); points are replicates, lines join group means. ISL1 transcripts are reduced in KO
(normalized counts {min(isl1_ko):.0f}-{max(isl1_ko):.0f}, vs {min(isl1_ctr):.0f}-{max(isl1_ctr):.0f} in Ctr) but
not absent; the 3'-end library cannot show where the gene was edited. * adjusted p &lt; {ALPHA}.</p>
{mark_html}
{marker_tbl}
</main></body></html>
"""
with open(OUT, "w") as fh:
    fh.write(page)
print("wrote", OUT)
print("clusters:", {c: int((clus == c).sum()) for c in range(1, N_CLUSTERS + 1)})
