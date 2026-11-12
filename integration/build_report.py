#!/usr/bin/env python3
"""Build the ISL1 CUT&RUN x bulk RNA-seq integration report.

Run on the Mac from integration/ after `bash scripts/fetch_results.sh`:
  ~/anaconda3/bin/python3 build_report.py
Reads the step 01-05 tables in results/ and, for section 06, the bulk normalized counts
(../bulkRNAseq_ISL1_KO_timecourse/results/DE_DESeq2/normalized_counts.tsv).
Output: results/report_integration.html
"""
import html
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
BULK_COUNTS = os.path.join(HERE, "..", "bulkRNAseq_ISL1_KO_timecourse", "results", "DE_DESeq2", "normalized_counts.tsv")
# Section 06 marker genes (human islet cell-type markers; pan-endocrine genes left out)
MARKERS = {
    "alpha": ["GCG", "ARX", "IRX1", "IRX2", "TTR", "GC", "LOXL4", "DPP4", "FAP", "TM4SF4", "CRYBA2", "PTPRT", "KCTD12",
              "MAFB", "POU3F4", "PLCE1", "TMEM236", "SMIM24", "SLC7A2"],
    "beta": ["INS", "IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "G6PC2", "HADH", "PCSK1"],
    "delta": ["SST", "HHEX", "LEPR", "RBP4"],
}
OUT = os.path.join(RES, "report_integration.html")

GROUPS = ["no_KO_effect", "cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "stays", "other_KO"]
LABEL = {
    "no_KO_effect": "No KO effect",
    "cluster1": "Cluster 1: down, strongest at D18",
    "cluster2": "Cluster 2: up early, down by D22",
    "cluster3": "Cluster 3: up, emerging from D18",
    "cluster4": "Cluster 4: down, emerging and growing",
    "cluster5": "Cluster 5: up, late (D22)",
    "cluster6": "Cluster 6: up early, down by D22",
    "stays": "Stays: one sign at all days",
    "other_KO": "Other KO effect",
}
CONTRASTS = ["D18_vs_D15", "D22_vs_D18", "D22_vs_D15"]
CT_LABEL = {"D18_vs_D15": "D18 vs D15", "D22_vs_D18": "D22 vs D18", "D22_vs_D15": "D22 vs D15"}
C_GROUP = "#2a78d6"   # series slot 1
C_OTHER = "#898781"   # muted
GRID = "#e1e0d9"
AXIS = "#c3c2b7"


def rd(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=["NA"])


def fig_html(fig, first):
    return fig.to_html(full_html=False, include_plotlyjs=True if first else False,
                       config={"displaylogo": False, "responsive": True})


def cell(v):
    if isinstance(v, (float, np.floating)):
        if np.isnan(v):
            return ""
        if float(v).is_integer() and abs(v) >= 1:
            return f"{int(v):,}"
        return f"{v:.3g}"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}"
    return str(v)


def table_html(df):
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = "".join("<tr>" + "".join(f"<td>{html.escape(cell(v))}</td>" for v in r) + "</tr>"
                   for r in df.itertuples(index=False))
    return f'<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table></div>'


def style(fig, height):
    fig.update_layout(height=height, margin=dict(l=60, r=20, t=60, b=50), plot_bgcolor="#fcfcfb",
                      paper_bgcolor="#fcfcfb", font=dict(family="system-ui, sans-serif", size=12, color="#2b2b2b"),
                      hoverlabel=dict(bgcolor="white"))
    fig.update_xaxes(gridcolor=GRID, linecolor=AXIS, zeroline=False)
    fig.update_yaxes(gridcolor=GRID, linecolor=AXIS, zeroline=False)
    return fig


def axref(a, k):
    return f"{a} domain" if k == 1 else f"{a}{k} domain"


def ecdf(x, n=200):
    x = np.sort(np.asarray(x))
    # thin to <= n points for the line, keeping the ends
    idx = np.unique(np.linspace(0, len(x) - 1, min(n, len(x))).round().astype(int))
    return x[idx], (idx + 1) / len(x)


# ---- data
great = rd(os.path.join(RES, "01_great", "summary.tsv"))
enr = rd(os.path.join(RES, "02_direct_targets", "enrichment.tsv"))
cam = rd(os.path.join(RES, "03_binding_change", "camera.tsv"))
st = rd(os.path.join(RES, "03_binding_change", "site_stats.tsv")).set_index("site")
sg = rd(os.path.join(RES, "03_binding_change", "site_groups.tsv"))
members = {g: set(sg.loc[sg.group == g, "site"]) for g in GROUPS}
parts = []

# ---- step 02: bound fraction per group
e = enr[enr.assignment == "great"].set_index("group").loc[GROUPS]
ref = float(e.loc["no_KO_effect", "pct_bound"])
fig = go.Figure(go.Bar(
    x=[LABEL[g] for g in GROUPS], y=e.pct_bound, marker_color=[C_OTHER] + [C_GROUP] * (len(GROUPS) - 1),
    marker_line_width=0, width=0.6,
    customdata=np.stack([e.n_bound, e.n_genes, e.adj_odds_ratio.fillna(1), e.adj_fdr.fillna(np.nan)], axis=1),
    hovertemplate="%{x}<br>%{y:.1f}% bound (%{customdata[0]} of %{customdata[1]})"
                  "<br>adjusted OR %{customdata[2]}, FDR %{customdata[3]}<extra></extra>"))
fig.add_hline(y=ref, line=dict(color=C_OTHER, dash="dash", width=1),
              annotation_text=f"no KO effect {ref:.1f}%", annotation_position="top left")
fig.update_yaxes(title="% of genes with >= 1 ISL1 site (GREAT)", range=[0, 100])
fig.update_layout(title="Share of genes bound by ISL1, per KO-effect group", showlegend=False)
fig.update_xaxes(tickangle=-25)
parts.append(("02", fig_html(style(fig, 460), True)))
t02 = enr.copy()
t02["group"] = t02.group.map(LABEL)

# ---- step 03: CDF small multiples (one row per group, one column per contrast)
fig = make_subplots(rows=len(GROUPS), cols=len(CONTRASTS), shared_xaxes=True, vertical_spacing=0.025, horizontal_spacing=0.05,
                    column_titles=[CT_LABEL[c] for c in CONTRASTS])
for i, g in enumerate(GROUPS, start=1):
    inset = st.index.isin(members[g])
    for j, ct in enumerate(CONTRASTS, start=1):
        lfc = st[f"db_log2FC_{ct}"]
        row = cam[(cam.contrast == ct) & (cam.group == g)].iloc[0]
        for vals, name, color, dash in [(lfc[~inset], "All other sites", C_OTHER, "dot"),
                                        (lfc[inset], "Sites near the group's genes", C_GROUP, "solid")]:
            x, y = ecdf(vals)
            fig.add_trace(go.Scatter(x=x, y=y, mode="lines", line=dict(color=color, width=2, dash=dash), name=name,
                                     legendgroup=name, showlegend=(i == 1 and j == 1),
                                     hovertemplate=f"{name}<br>log2FC %{{x:.2f}}<br>cumulative share %{{y:.2f}}<extra></extra>"),
                          row=i, col=j)
        fig.add_annotation(text=f"{row.camera_direction}, FDR {row.camera_fdr:.2g} (n={row.n_sites})",
                           xref=axref("x", (i - 1) * len(CONTRASTS) + j), yref=axref("y", (i - 1) * len(CONTRASTS) + j),
                           x=0.02, y=0.95, showarrow=False, xanchor="left", font=dict(size=11, color="#555"))
    fig.update_yaxes(title_text=LABEL[g].split(":")[0], row=i, col=1, title_font=dict(size=11))
for j in range(1, len(CONTRASTS) + 1):
    fig.update_xaxes(title_text="ISL1 binding log2FC (DiffBind, background)", range=[-2, 2.5], row=len(GROUPS), col=j)
fig.update_layout(title="Binding change at sites near each group's genes vs all other sites (cumulative distribution)",
                  legend=dict(orientation="h", y=1.03, x=0, yanchor="bottom"))
fig = style(fig, 170 * len(GROUPS) + 160)
fig.update_layout(margin=dict(t=120), title=dict(y=0.995, yanchor="top"))
fig.update_annotations(selector=dict(text=CT_LABEL["D18_vs_D15"]), yshift=-2)
parts.append(("03", fig_html(fig, False)))
t03 = cam[["contrast", "group", "n_sites", "median_log2FC", "median_log2FC_other", "pct_gained", "pct_gained_other",
           "camera_direction", "camera_p", "camera_fdr", "est_correlation", "est_fdr", "wilcoxon_p"]].copy()
t03["group"] = t03.group.map(LABEL)

# ---- step 04: candidates
c_sum = rd(os.path.join(RES, "04_candidates", "summary.tsv"))
c_sum["group"] = c_sum.group.map(LABEL)
c_enr = rd(os.path.join(RES, "04_candidates", "gained_enrichment.tsv"))
c_enr["group"] = c_enr.group.map(LABEL)
c_gen = rd(os.path.join(RES, "04_candidates", "candidate_genes.tsv"))
c_gen["group"] = c_gen.group.map(LABEL)
for c in ["ko_lfc_D15", "ko_lfc_D18", "ko_lfc_D22"]:
    c_gen[c] = c_gen[c].round(2)

# ---- step 05: binding signal per day at each cluster's sites (boxplots)
SIG_GROUPS = ["no_KO_effect", "cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6"]
DAY_COL = {"D15": "#86b6ef", "D18": "#2a78d6", "D22": "#184f95"}   # one-hue sequential ramp, light -> dark
gsig = rd(os.path.join(RES, "05_binding_signal", "group_signal.tsv"))
wt = rd(os.path.join(RES, "05_binding_signal", "wilcoxon.tsv"))
xlab = {g: LABEL[g].replace(": ", "<br>") for g in SIG_GROUPS}


def stars(q):
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


def box_stats(v):
    v = np.asarray(v)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    lo, hi = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
    return q1, med, q3, v[v >= lo].min(), v[v <= hi].max()   # whiskers: furthest site within 1.5 x IQR


fig = go.Figure()
ymax = 0
for d, col in DAY_COL.items():
    bs = [box_stats(gsig.loc[gsig.group == g, f"signal_{d}"]) for g in SIG_GROUPS]
    q1, med, q3, lf, uf = (list(t) for t in zip(*bs))
    ymax = max(ymax, max(uf))
    fig.add_trace(go.Box(x=[xlab[g] for g in SIG_GROUPS], q1=q1, median=med, q3=q3, lowerfence=lf, upperfence=uf,
                         name=d, marker_color=col, line=dict(width=1.5), offsetgroup=d,
                         fillcolor="rgba({},{},{},0.45)".format(*(int(col[i:i + 2], 16) for i in (1, 3, 5))),
                         whiskerwidth=0.4, hoverinfo="y+name"))
for g in SIG_GROUPS:
    w = wt[wt.group == g].set_index("comparison")
    n = int(w.n_sites.iloc[0])
    txt = (f"D18/D15 {stars(w.loc['D18 vs D15', 'fdr'])}<br>D22/D18 {stars(w.loc['D22 vs D18', 'fdr'])}"
           f"<br><span style='color:#898781'>n = {n:,}</span>")
    fig.add_annotation(x=xlab[g], y=ymax * 1.02, text=txt, showarrow=False, yanchor="bottom", font=dict(size=11))
fig.update_layout(boxmode="group", boxgap=0.25, boxgroupgap=0.08,
                  title="ISL1 binding signal at the sites linked to each cluster's genes",
                  legend=dict(title="Day", orientation="h", y=1.08, x=1, xanchor="right", yanchor="bottom"))
fig.update_yaxes(title="ISL1 CUT&RUN signal (normalized reads per site)", range=[0, ymax * 1.18])
parts.append(("05", fig_html(style(fig, 560), False)))
t05 = wt.copy()
t05["group"] = t05.group.map(LABEL)

# ---- section 06: alpha, beta and delta programs (binding vs expression over time)
DAYS = ["D15", "D18", "D22"]
C_EXPR = "#eb6834"   # categorical slot 2
site_sig = rd(os.path.join(RES, "05_binding_signal", "site_signal.tsv")).set_index("site")
sgl = rd(os.path.join(RES, "01_great", "site_gene.tsv"))
gtab = rd(os.path.join(RES, "02_direct_targets", "gene_table.tsv")).set_index("gene_name")
bulk = rd(BULK_COUNTS)
LINE = {"Ctr": "_Ctr-H3_", "KO": "_ISL1-KO_"}
expr_cols = {f"{g}_{d}": [c for c in bulk.columns if f"{LINE[g]}{d}_" in c] for g in LINE for d in DAYS}
assert all(len(v) == 2 for v in expr_cols.values()), expr_cols
expr = pd.DataFrame({k: bulk[v].mean(axis=1) for k, v in expr_cols.items()})
expr.index = bulk.gene_name
expr = expr[~expr.index.duplicated()]
cand_names = set(c_gen.gene_name)

rows = []
for ctype, genes in MARKERS.items():
    for g in genes:
        sites = sgl.loc[sgl.gene_name == g, "site"].unique()
        b = site_sig.loc[sites, [f"signal_{d}" for d in DAYS]].sum() if len(sites) else pd.Series([np.nan] * 3)
        r = {"cell_type": ctype, "gene": g, "group": LABEL.get(gtab.group.get(g, ""), "not tested"), "n_sites": len(sites)}
        for d, v in zip(DAYS, b.to_numpy()):
            r[f"binding_{d}"] = round(float(v), 1) if len(sites) else np.nan
        for d in DAYS:
            r[f"ctr_expr_{d}"] = round(float(expr.loc[g, f"Ctr_{d}"])) if g in expr.index else np.nan
        for d in DAYS:
            r[f"ko_log2FC_{d}"] = round(float(gtab.loc[g, f"ko_lfc_{d}"]), 2) if g in gtab.index else np.nan
        r["candidate"] = g in cand_names
        rows.append(r)
mk = pd.DataFrame(rows)

# panel A: binding signal per day at the sites of each marker set, all sites as reference
sets = {"All sites": site_sig.index.to_numpy()}
for ctype, genes in MARKERS.items():
    sets[f"{ctype.capitalize()} genes"] = sgl.loc[sgl.gene_name.isin(genes), "site"].unique()
fig = go.Figure()
ymax = 0
for d, col in DAY_COL.items():
    bs = [box_stats(site_sig.loc[s, f"signal_{d}"]) for s in sets.values()]
    q1, med, q3, lf, uf = (list(t) for t in zip(*bs))
    ymax = max(ymax, max(uf))
    fig.add_trace(go.Box(x=list(sets), q1=q1, median=med, q3=q3, lowerfence=lf, upperfence=uf, name=d,
                         marker_color=col, line=dict(width=1.5), offsetgroup=d, whiskerwidth=0.4, hoverinfo="y+name",
                         fillcolor="rgba({},{},{},0.45)".format(*(int(col[i:i + 2], 16) for i in (1, 3, 5)))))
wrows = []
for name, s in sets.items():
    x = site_sig.loc[s]
    for a, b in [("D18", "D15"), ("D22", "D18"), ("D22", "D15")]:
        wrows.append({"set": name, "comparison": f"{a} vs {b}", "n_sites": len(s),
                      "median_a": round(x[f"signal_{a}"].median(), 1), "median_b": round(x[f"signal_{b}"].median(), 1),
                      "pct_sites_up": round(100 * (x[f"signal_{a}"] > x[f"signal_{b}"]).mean(), 1),
                      "p": wilcoxon(x[f"signal_{a}"], x[f"signal_{b}"]).pvalue})
w06 = pd.DataFrame(wrows)
order = np.argsort(w06.p.to_numpy())
q = np.empty(len(order))
q[order] = np.minimum.accumulate((w06.p.to_numpy()[order] * len(order) / np.arange(1, len(order) + 1))[::-1])[::-1]
w06["fdr"] = np.minimum(q, 1)
for name, s in sets.items():
    w = w06[w06.set == name].set_index("comparison")
    txt = (f"D18/D15 {stars(w.loc['D18 vs D15', 'fdr'])}<br>D22/D18 {stars(w.loc['D22 vs D18', 'fdr'])}"
           f"<br><span style='color:#898781'>n = {len(s):,}</span>")
    fig.add_annotation(x=name, y=ymax * 1.02, text=txt, showarrow=False, yanchor="bottom", font=dict(size=11))
fig.update_layout(boxmode="group", boxgap=0.3, boxgroupgap=0.08,
                  title="ISL1 binding signal at the sites linked to islet cell-type marker genes",
                  legend=dict(title="Day", orientation="h", y=1.08, x=1, xanchor="right", yanchor="bottom"))
fig.update_yaxes(title="ISL1 CUT&RUN signal (normalized reads per site)", range=[0, ymax * 1.2])
parts.append(("06a", fig_html(style(fig, 480), False)))
w06["p"] = w06.p.map(lambda v: float(f"{v:.3g}"))
w06["fdr"] = w06.fdr.map(lambda v: float(f"{v:.3g}"))

# panel B: per bound alpha gene, binding (sum over its sites) and control expression, log2 change from D15
bound_alpha = mk[(mk.cell_type == "alpha") & (mk.n_sites > 0)].reset_index(drop=True)
ncol = 4
nrow = int(np.ceil(len(bound_alpha) / ncol))
fig = make_subplots(rows=nrow, cols=ncol, subplot_titles=[f"{g} ({n} site{'s' if n > 1 else ''})"
                                                          for g, n in zip(bound_alpha.gene, bound_alpha.n_sites)],
                    shared_yaxes=True, vertical_spacing=0.12, horizontal_spacing=0.04)
for k, r in bound_alpha.iterrows():
    i, j = k // ncol + 1, k % ncol + 1
    bind = np.log2(np.array([r[f"binding_{d}"] for d in DAYS]) + 1)
    ex = np.log2(np.array([r[f"ctr_expr_{d}"] for d in DAYS]) + 1)
    for y, name, col in [(bind - bind[0], "ISL1 binding (sum of its sites)", C_GROUP),
                         (ex - ex[0], "Expression in control (RNA-seq)", C_EXPR)]:
        fig.add_trace(go.Scatter(x=DAYS, y=y, mode="lines+markers", name=name, legendgroup=name,
                                 showlegend=(k == 0), line=dict(color=col, width=2), marker=dict(size=8),
                                 hovertemplate=f"{r.gene} {name}<br>%{{x}}: log2 change %{{y:.2f}}<extra></extra>"),
                      row=i, col=j)
fig.update_yaxes(title_text="log2 change from D15", col=1)
fig.update_layout(title="Alpha genes: ISL1 binding near the gene vs its expression in control islets",
                  legend=dict(orientation="h", y=1.0, x=0, yanchor="bottom"))
fig = style(fig, 230 * nrow + 170)
fig.update_layout(margin=dict(t=150), title=dict(y=0.995, yanchor="top"))
fig.update_layout(legend=dict(y=1.0 + 75 / (230 * nrow + 170 - 150 - 50)))
parts.append(("06b", fig_html(fig, False)))
mk_show = mk.copy()
mk_show["candidate"] = mk_show.candidate.map({True: "yes", False: ""})

# ---- page
doc = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>ISL1 integration</title>
<style>
body{{font-family:system-ui,sans-serif;background:#fcfcfb;color:#2b2b2b;max-width:1150px;margin:0 auto;padding:16px;line-height:1.45}}
h1{{font-size:1.5em}} h2{{font-size:1.2em;margin-top:2em;border-bottom:1px solid #e1e0d9;padding-bottom:4px}}
.tw{{overflow-x:auto}} table{{border-collapse:collapse;font-size:12px;margin:8px 0}}
th,td{{border-bottom:1px solid #e1e0d9;padding:3px 8px;text-align:left;white-space:nowrap}} th{{background:#f1f0eb}}
p.note{{color:#555;font-size:13px}}
</style></head><body>
<h1>ISL1 CUT&amp;RUN x bulk RNA-seq: direct targets and time-dependent KO effect</h1>
<p class="note">Bulk RNA-seq: Ctr vs ISL1 KO at D15 / D18 / D22 (run YJZ5LS); KO-effect groups are the 6 clusters of the
LRT_interaction genes from the bulk report, genes whose KO effect stays, other KO genes, and genes without a KO effect.
ISL1 CUT&amp;RUN: control line, 3 libraries per day, 14,431 sites, linked to genes with GREAT default domains.</p>

<h2>01 GREAT domains, site to gene</h2>
{table_html(great)}

<h2>02 Share of genes bound by ISL1 per KO-effect group</h2>
<p class="note">Bound = at least one site in the gene's GREAT domain. Each group against genes without a KO effect: Fisher's
exact test and logistic regression with domain length (adjusted OR, FDR). Check: a site whose nearest TSS is the gene's,
within 50 kb.</p>
{parts[0][1]}
{table_html(t02)}

<h2>03 Binding change near each KO-effect group</h2>
<p class="note">Sites near a group's genes (GREAT) against all other sites. Lines: cumulative distribution of the step 27
DiffBind log2FC (background normalization); a curve shifted right means binding rises more than at other sites. Test:
limma voom on the 9 libraries (~ batch + day) and camera per group (competitive; main FDR with limma's default
inter-site correlation 0.01, est_fdr with the correlation estimated from the data). The Wilcoxon p treats sites as
independent and is descriptive only. Gained = FDR &lt; 0.05 up in D18 vs D15 or D22 vs D15.</p>
{parts[1][1]}
{table_html(t03)}

<h2>04 Candidate time-dependent direct targets</h2>
<p class="note">A candidate has a site in its GREAT domain whose binding change matches the timing of its KO effect:
cluster 1 (KO effect strongest at D18) a site gained D18 vs D15; clusters 3 and 4 (KO effect emerging from D18) a site
gained D18 or D22 vs D15; cluster 5 (late) a site gained D22 vs D15; clusters 2 and 6 (KO effect fading) a site lost.
Gained / lost: step 27 DiffBind, background normalization, FDR &lt; 0.05. ko_sig_on_change_day: the gene's own KO
effect is significant (padj &lt; 0.05) on the day of the matching change.</p>
{table_html(c_sum)}
<p class="note">Genes with at least one gained site per group against genes without a KO effect (logistic regression with
domain length and number of sites in the domain).</p>
{table_html(c_enr)}
<p class="note">Candidate genes, one row each (best matching site; full list of sites in results/04_candidates/candidates.tsv).</p>
{table_html(c_gen)}

<h2>05 ISL1 binding signal per day at each cluster's sites</h2>
<p class="note">Sites linked (GREAT) to the genes of each cluster, and to genes without a KO effect as the reference.
Signal: reads per site divided by the step 27 background size factors, mean of the 3 libraries per day; boxes show
the median and quartiles, whiskers 1.5 x IQR, outliers not drawn. Stars: paired Wilcoxon signed-rank test between days
within the group (BH; *** FDR &lt; 0.001, ** &lt; 0.01, * &lt; 0.05). Binding rises after D15 at most sites, so these
tests show that binding changes at a group's sites, not that it changes more than at other sites; that comparison is
the camera test in section 03 (cluster 4: up D18 vs D15, FDR 0.23; down D22 vs D18, FDR 0.023; both not
significant with the estimated inter-site correlation).</p>
{parts[2][1]}
{table_html(t05)}

<h2>06 Alpha, beta and delta programs: ISL1 binding and expression over time</h2>
<p class="note">Marker genes (human islet cell-type markers, pan-endocrine genes left out): alpha {", ".join(MARKERS["alpha"])};
beta {", ".join(MARKERS["beta"])}; delta {", ".join(MARKERS["delta"])}. Sites linked to them by GREAT; signal as in
section 05. Stars: paired Wilcoxon between days within the set (BH over this panel), descriptive as in section 05.</p>
{parts[3][1]}
{table_html(w06)}
<p class="note">Each alpha gene with at least one linked site: ISL1 binding (sum of the normalized signal of its sites)
and its mean normalized expression in control islets (Ctr-H3), both as log2 change from D15.</p>
{parts[4][1]}
<p class="note">All marker genes: KO-effect group, linked sites, ISL1 binding (sum over sites), control expression
(mean normalized counts) and the KO effect (shrunk log2FC KO vs Ctr) per day; candidate = in the section 04 list.</p>
{table_html(mk_show)}
</body></html>"""
with open(OUT, "w") as f:
    f.write(doc)
print("wrote", OUT)
