#!/usr/bin/env python3
"""Build the ISL1 CUT&RUN x bulk RNA-seq integration report.

Run on the Mac from integration/ after `bash scripts/fetch_results.sh`:
  ~/anaconda3/bin/python3 build_report.py
Reads results/01_great/summary.tsv, results/02_direct_targets/enrichment.tsv,
results/03_binding_change/{camera,site_stats,site_groups}.tsv.
Output: results/report_integration.html
"""
import html
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
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
CONTRASTS = ["D18_vs_D15", "D22_vs_D15"]
CT_LABEL = {"D18_vs_D15": "D18 vs D15", "D22_vs_D15": "D22 vs D15"}
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
fig = make_subplots(rows=len(GROUPS), cols=2, shared_xaxes=True, vertical_spacing=0.025, horizontal_spacing=0.06,
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
                           xref=f"x{(i - 1) * 2 + j} domain".replace("x1 domain", "x domain"),
                           yref=f"y{(i - 1) * 2 + j} domain".replace("y1 domain", "y domain"),
                           x=0.02, y=0.95, showarrow=False, xanchor="left", font=dict(size=11, color="#555"))
    fig.update_yaxes(title_text=LABEL[g].split(":")[0], row=i, col=1, title_font=dict(size=11))
for j in (1, 2):
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
</body></html>"""
with open(OUT, "w") as f:
    f.write(doc)
print("wrote", OUT)
