#!/usr/bin/env python3
"""Build the ISL1 CUT&RUN x bulk RNA-seq integration report.

Run on the Mac from integration/ after `bash scripts/fetch_results.sh`:
  ~/anaconda3/bin/python3 build_report.py
Reads the step 01-08 tables in results/ and, for section 06, the bulk normalized counts
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
# Section 06 marker genes: markers.tsv (cell_type, gene), shared with step 06
_mk = pd.read_csv(os.path.join(HERE, "markers.tsv"), sep="\t", comment="#")
MARKERS = {ct: list(g.gene) for ct, g in _mk.groupby("cell_type", sort=False)}

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

# ---- section 06: cluster 4 (alpha-leaning) vs cluster 2 (beta-leaning)
DAYS = ["D15", "D18", "D22"]
C_EXPR = "#eb6834"   # categorical slot 2
CLU = {"cluster4": "Cluster 4 (alpha-leaning)", "cluster2": "Cluster 2 (beta-leaning)"}
site_sig = rd(os.path.join(RES, "05_binding_signal", "site_signal.tsv")).set_index("site")
sgl = rd(os.path.join(RES, "01_great", "site_gene.tsv"))
gtab = rd(os.path.join(RES, "02_direct_targets", "gene_table.tsv")).set_index("gene_name")
gtab = gtab[~gtab.index.duplicated()]   # a few gene names occur twice; markers are unique
bulk = rd(BULK_COUNTS)
LINE = {"Ctr": "_Ctr-H3_", "KO": "_ISL1-KO_"}
expr_cols = {f"{g}_{d}": [c for c in bulk.columns if f"{LINE[g]}{d}_" in c] for g in LINE for d in DAYS}
assert all(len(v) == 2 for v in expr_cols.values()), expr_cols
expr = pd.DataFrame({k: bulk[v].mean(axis=1) for k, v in expr_cols.items()})
expr.index = bulk.gene_name
expr = expr[~expr.index.duplicated()]
cand_names = set(c_gen.gene_name)
lib = rd(os.path.join(RES, "06_cluster_contrast", "library_values.tsv"))
ctest = rd(os.path.join(RES, "06_cluster_contrast", "tests.tsv"))
set_n = rd(os.path.join(RES, "06_cluster_contrast", "set_sizes.tsv")).set_index("set").n_sites

# marker table: every marker gene with its group, sites, binding, control expression and KO effect
rows = []
for ctype, genes in MARKERS.items():
    for g in genes:
        sites = sgl.loc[sgl.gene_name == g, "site"].unique()
        b = site_sig.loc[sites, [f"signal_{d}" for d in DAYS]].sum() if len(sites) else None
        r = {"cell_type": ctype, "gene": g, "group": LABEL.get(gtab.group.get(g, ""), "not tested"),
             "n_sites": len(sites)}
        for k, d in enumerate(DAYS):
            r[f"binding_{d}"] = round(float(b.iloc[k]), 1) if b is not None else np.nan
        for d in DAYS:
            r[f"ctr_expr_{d}"] = round(float(expr.loc[g, f"Ctr_{d}"])) if g in expr.index else np.nan
        for d in DAYS:
            r[f"ko_log2FC_{d}"] = round(float(gtab.loc[g, f"ko_lfc_{d}"]), 2) if g in gtab.index else np.nan
        r["candidate"] = g in cand_names
        rows.append(r)
mk = pd.DataFrame(rows)
mk["cluster"] = mk.gene.map(gtab.group)

# panel A: binding signal per site at alpha vs beta marker sites, per day (boxplots as in section 05);
# annotation: the library-level test of the median per site (step 06, alpha_vs_beta, level_<day>)
CMP = {"cluster4_vs_cluster2": "Cluster 4 / cluster 2 sites",
       "alpha_vs_beta": "Alpha / beta marker sites",
       "alpha_in_c4_vs_beta_in_c2": "Alpha markers in cluster 4 / beta markers in cluster 2",
       "cluster4_vs_all": "Cluster 4 / all sites", "cluster2_vs_all": "Cluster 2 / all sites"}
C_PAIR = ["#2a78d6", "#eb6834"]   # categorical slots 1 and 2: alpha-leaning, beta-leaning
tt = ctest.set_index(["comparison", "metric", "term"])


def pair_boxplot(site_sets, comparison, title):
    """Signal per site at two site sets, per day; stars from the step 06 library-level test of the median per site."""
    fig = go.Figure()
    ymax = 0
    for (name, sites), col in zip(site_sets.items(), C_PAIR):
        bs = [box_stats(site_sig.loc[sites, f"signal_{d}"]) for d in DAYS]
        q1, med, q3, lf, uf = (list(t) for t in zip(*bs))
        ymax = max(ymax, max(uf))
        fig.add_trace(go.Box(x=DAYS, q1=q1, median=med, q3=q3, lowerfence=lf, upperfence=uf,
                             name=f"{name} ({len(sites):,} sites)", marker_color=col, line=dict(width=1.5),
                             offsetgroup=name, whiskerwidth=0.4, hoverinfo="y+name",
                             fillcolor="rgba({},{},{},0.45)".format(*(int(col[i:i + 2], 16) for i in (1, 3, 5)))))
    for d in DAYS:
        r = tt.loc[(comparison, "median", f"level_{d}")]
        fig.add_annotation(x=d, y=ymax * 1.02, showarrow=False, yanchor="bottom", font=dict(size=11),
                           text=f"{stars(r.p)}<br><span style='color:#898781'>log2 ratio {r.estimate:+.2f}, "
                                f"p {r.p:.2g}</span>")
    fig.update_layout(boxmode="group", boxgap=0.35, boxgroupgap=0.1, title=title,
                      legend=dict(orientation="h", y=1.08, x=1, xanchor="right", yanchor="bottom"))
    fig.update_yaxes(title="ISL1 CUT&RUN signal (normalized reads per site)", range=[0, ymax * 1.2])
    return fig_html(style(fig, 500), False)


marker_sites = {f"{ct.capitalize()} marker genes": sgl.loc[sgl.gene_name.isin(MARKERS[ct]), "site"].unique()
                for ct in ("alpha", "beta")}
parts.append(("06a", pair_boxplot(marker_sites, "alpha_vs_beta",
                                  "ISL1 binding signal at alpha and beta marker gene sites, per day")))
cluster_sites = {"Cluster 4, alpha-leaning": sg.loc[sg.group == "cluster4", "site"].unique(),
                 "Cluster 2, beta-leaning": sg.loc[sg.group == "cluster2", "site"].unique()}
parts.append(("06a2", pair_boxplot(cluster_sites, "cluster4_vs_cluster2",
                                   "ISL1 binding signal at cluster 4 and cluster 2 sites, per day")))
ct_show = ctest.copy()
ct_show["comparison"] = ct_show.comparison.map(CMP)
# panel B: marker genes of each cluster with >= 1 site: binding (sum of its sites) and control expression,
# log2 change from D15; alpha markers of cluster 4 first, then beta markers of cluster 2
panel = pd.concat([mk[(mk.cell_type == "alpha") & (mk.cluster == "cluster4") & (mk.n_sites > 0)],
                   mk[(mk.cell_type == "beta") & (mk.cluster == "cluster2") & (mk.n_sites > 0)]]).reset_index(drop=True)
ncol = 4
nrow = int(np.ceil(len(panel) / ncol))
titles = [f"{'alpha' if r.cell_type == 'alpha' else 'beta'} · {r.gene} ({r.n_sites} site{'s' if r.n_sites > 1 else ''})"
          for r in panel.itertuples()]
fig = make_subplots(rows=nrow, cols=ncol, subplot_titles=titles, shared_yaxes=True, vertical_spacing=0.1,
                    horizontal_spacing=0.04)
for k, r in panel.iterrows():
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
fig.update_layout(title="Marker genes: ISL1 binding near the gene vs its expression in control islets "
                        "(alpha markers in cluster 4, then beta markers in cluster 2)",
                  legend=dict(orientation="h", y=1.0, x=0, yanchor="bottom"))
fig = style(fig, 230 * nrow + 170)
fig.update_layout(margin=dict(t=150), title=dict(y=0.995, yanchor="top"))
fig.update_layout(legend=dict(y=1.0 + 75 / (230 * nrow + 170 - 150 - 50)))
parts.append(("06b", fig_html(fig, False)))
n_alpha4 = int(((mk.cell_type == "alpha") & (mk.cluster == "cluster4")).sum())
n_beta4 = int(((mk.cell_type == "beta") & (mk.cluster == "cluster4")).sum())
n_alpha2 = int(((mk.cell_type == "alpha") & (mk.cluster == "cluster2")).sum())
n_beta2 = int(((mk.cell_type == "beta") & (mk.cluster == "cluster2")).sum())
mk_show = mk.drop(columns="cluster")
mk_show["candidate"] = mk_show.candidate.map({True: "yes", False: ""})

# ---- section 07: target genes of the differentially bound sites and their KO effect
s07 = rd(os.path.join(RES, "07_diff_site_targets", "summary.tsv")).set_index("binding")
t07 = rd(os.path.join(RES, "07_diff_site_targets", "tests.tsv"))
d07 = rd(os.path.join(RES, "07_diff_site_targets", "diff_site_genes.tsv"))
CLS = {"gained": ("Gained site", "#2a78d6"), "stable": ("Only stable sites", "#eb6834"),
       "no site": ("No site", "#1baf7a")}   # categorical slots 1-3
fig = make_subplots(rows=1, cols=2, shared_yaxes=True, horizontal_spacing=0.06,
                    subplot_titles=["Lower in the KO (ISL1 activates)", "Higher in the KO (ISL1 represses)"])
for j, direction in enumerate(["down", "up"], start=1):
    for cls, (name, col) in CLS.items():
        y = [s07.loc[cls, f"pct_{direction}_{d}"] for d in DAYS]
        fig.add_trace(go.Bar(x=DAYS, y=y, name=f"{name} ({int(s07.loc[cls, 'n_genes']):,} genes)", marker_color=col,
                             marker_line_width=0, legendgroup=cls, showlegend=(j == 1), offsetgroup=cls,
                             hovertemplate=f"{name}<br>%{{x}}: %{{y:.1f}}% of genes<extra></extra>"), row=1, col=j)
fig.update_yaxes(title_text="% of genes with a KO effect that day (padj < 0.05)", col=1)
fig.update_layout(barmode="group", bargap=0.3, bargroupgap=0.08,
                  title="KO effect per day of the genes near sites that gain ISL1 binding",
                  legend=dict(orientation="h", y=1.0, x=1, xanchor="right", yanchor="bottom"))
fig = style(fig, 460)
fig.update_layout(margin=dict(t=120), title=dict(y=0.97, yanchor="top"), legend=dict(y=1.16))
parts.append(("07", fig_html(fig, False)))
s07_show = s07.reset_index()
# gained-site genes with a KO effect, one row per gene (most significant gained site)
gg = d07[d07.change.str.startswith("gained") & (d07.padj_genotype < 0.05)].copy()
gg["best_fdr"] = gg[["db_FDR_D18_vs_D15", "db_FDR_D22_vs_D15"]].min(axis=1)
gg = gg.sort_values(["gene_id", "best_fdr"], kind="stable").groupby("gene_id", sort=False).first().reset_index()
gg["n_gained_sites"] = gg.gene_id.map(d07[d07.change.str.startswith("gained")].groupby("gene_id").site.nunique())
gg["group"] = gg.group.map(LABEL)
gg = gg.sort_values(["ko_padj_D22", "gene_name"], kind="stable")
g07 = gg[["gene_name", "group", "n_gained_sites", "site", "change", "distance_to_tss", "db_log2FC_D18_vs_D15",
          "db_log2FC_D22_vs_D15", "ko_lfc_D15", "ko_lfc_D18", "ko_lfc_D22", "ko_padj_D22"]].copy()
for c in ["db_log2FC_D18_vs_D15", "db_log2FC_D22_vs_D15", "ko_lfc_D15", "ko_lfc_D18", "ko_lfc_D22"]:
    g07[c] = g07[c].round(2)
g07.columns = ["gene", "KO-effect group", "gained sites", "best site", "change", "distance to TSS",
               "site log2FC D18/D15", "site log2FC D22/D15", "KO log2FC D15", "KO log2FC D18", "KO log2FC D22",
               "KO padj D22"]

# boxplot: binding signal per day at the gained sites of the 343 genes (gained site + KO effect), split by the
# direction of the KO effect (on the latest day with a significant effect), with all sites as the reference
def ko_direction(r):
    for d in ("D22", "D18", "D15"):
        if r[f"ko_padj_{d}"] < 0.05:
            return "lower" if r[f"ko_lfc_{d}"] < 0 else "higher"
    return "no single day"


gdir = gg.set_index("gene_id").apply(ko_direction, axis=1)
gs = d07[d07.change.str.startswith("gained") & d07.gene_id.isin(gg.gene_id)].copy()
gs["direction"] = gs.gene_id.map(gdir)
sets07 = {"All sites": site_sig.index.to_numpy(),
          f"Gained sites, all {len(gg)} genes": gs.site.unique(),
          f"Genes lower in KO ({int((gdir == 'lower').sum())})": gs.loc[gs.direction == "lower", "site"].unique(),
          f"Genes higher in KO ({int((gdir == 'higher').sum())})": gs.loc[gs.direction == "higher", "site"].unique()}
fig = go.Figure()
ymax = 0
for d, col in DAY_COL.items():
    bs = [box_stats(site_sig.loc[v, f"signal_{d}"]) for v in sets07.values()]
    q1, med, q3, lf, uf = (list(t) for t in zip(*bs))
    ymax = max(ymax, max(uf))
    fig.add_trace(go.Box(x=list(sets07), q1=q1, median=med, q3=q3, lowerfence=lf, upperfence=uf, name=d,
                         marker_color=col, line=dict(width=1.5), offsetgroup=d, whiskerwidth=0.4, hoverinfo="y+name",
                         fillcolor="rgba({},{},{},0.45)".format(*(int(col[i:i + 2], 16) for i in (1, 3, 5)))))
w07 = []
for name, v in sets07.items():
    x = site_sig.loc[v]
    for a, b in [("D18", "D15"), ("D22", "D18"), ("D22", "D15")]:
        w07.append({"set": name, "comparison": f"{a} vs {b}", "n_sites": len(v),
                    "median_a": round(x[f"signal_{a}"].median(), 1), "median_b": round(x[f"signal_{b}"].median(), 1),
                    "pct_sites_up": round(100 * (x[f"signal_{a}"] > x[f"signal_{b}"]).mean(), 1),
                    "p": float(f"{wilcoxon(x[f'signal_{a}'], x[f'signal_{b}']).pvalue:.3g}")})
w07 = pd.DataFrame(w07)
for name, v in sets07.items():
    w = w07[w07.set == name].set_index("comparison")
    fig.add_annotation(x=name, y=ymax * 1.02, showarrow=False, yanchor="bottom", font=dict(size=11),
                       text=f"D18/D15 {stars(w.loc['D18 vs D15', 'p'])}<br>D22/D18 {stars(w.loc['D22 vs D18', 'p'])}"
                            f"<br><span style='color:#898781'>n = {len(v):,}</span>")
fig.update_layout(boxmode="group", boxgap=0.3, boxgroupgap=0.08,
                  title="ISL1 binding signal at the gained sites of the genes with a KO effect",
                  legend=dict(title="Day", orientation="h", y=1.08, x=1, xanchor="right", yanchor="bottom"))
fig.update_yaxes(title="ISL1 CUT&RUN signal (normalized reads per site)", range=[0, ymax * 1.2])
parts.append(("07b", fig_html(style(fig, 520), False)))

# ---- section 08: all sites grouped by binding pattern across days (k-means, k = 4)
D8 = os.path.join(RES, "08_site_clusters")
K8 = 4
SC_LABEL = {0: "Stable (LRT p ≥ 0.05)", 1: "1: rises by D18, stays", 2: "2: peaks at D18, falls by D22",
            3: "3: late, rises at D22", 4: "4: falls at D22"}
C_B3 = "#eb6834"   # series slot 2
sc8 = rd(os.path.join(D8, "site_clusters.tsv"))
cen8 = rd(os.path.join(D8, "centers.tsv"))
for c_ in cen8.columns:
    cen8[c_] = pd.to_numeric(cen8[c_], errors="coerce")
elb8 = rd(os.path.join(D8, "elbow.tsv"))
ho8 = rd(os.path.join(D8, "heldout.tsv"))
hoc8 = rd(os.path.join(D8, "heldout_clusters.tsv"))
lk8 = rd(os.path.join(D8, "links.tsv"))
zc = [f"z_{d}" for d in DAYS]

# 08a: z-scored profiles per group (median line, 25th-75th percentile band)
fig = make_subplots(rows=1, cols=K8, shared_yaxes=True, horizontal_spacing=0.04,
                    subplot_titles=[f"{SC_LABEL[c]}<br>n = {int((sc8[f'k{K8}'] == c).sum()):,}" for c in range(1, K8 + 1)])
for c in range(1, K8 + 1):
    z = sc8.loc[sc8[f"k{K8}"] == c, zc].to_numpy()
    q1, med, q3 = np.percentile(z, [25, 50, 75], axis=0)
    fig.add_trace(go.Scatter(x=DAYS + DAYS[::-1], y=np.r_[q3, q1[::-1]], mode="lines", fill="toself", fillcolor="rgba(42,120,214,0.18)",
                             line=dict(width=0), hoverinfo="skip", showlegend=False), row=1, col=c)
    fig.add_trace(go.Scatter(x=DAYS, y=med, mode="lines+markers", line=dict(color=C_GROUP, width=2), marker=dict(size=8),
                             showlegend=False, customdata=np.c_[q1, q3],
                             hovertemplate="%{x}: median z %{y:.2f} (IQR %{customdata[0]:.2f} to %{customdata[1]:.2f})<extra></extra>"),
                  row=1, col=c)
fig.update_yaxes(title_text="ISL1 binding, z-score across days", col=1)
fig.update_layout(title="ISL1 sites grouped by binding pattern (2,747 sites with LRT p < 0.05; k-means, k = 4)")
fig = style(fig, 400)
fig.update_layout(margin=dict(t=100))
parts.append(("08a", fig_html(fig, False)))

# 08b: held-out check, 6-library cluster centre vs mean batch 3 profile
h = hoc8[hoc8.k == K8].reset_index(drop=True)
fig = make_subplots(rows=1, cols=K8, shared_yaxes=True, horizontal_spacing=0.04,
                    subplot_titles=[f"6-library group {r.cluster}, n = {int(r.n_sites):,}<br>r = {r.pearson_r:.2f}"
                                    for r in h.itertuples()])
for i, r in h.iterrows():
    for cols_, name, col in [(["centre_D15", "centre_D18", "centre_D22"], "6 earlier libraries: group centre", C_GROUP),
                             (["batch3_D15", "batch3_D18", "batch3_D22"], "3 batch 3 libraries: mean of the group's sites", C_B3)]:
        fig.add_trace(go.Scatter(x=DAYS, y=[r[c] for c in cols_], mode="lines+markers", name=name, legendgroup=name,
                                 showlegend=(i == 0), line=dict(color=col, width=2), marker=dict(size=8),
                                 hovertemplate=f"{name}<br>%{{x}}: z %{{y:.2f}}<extra></extra>"), row=1, col=i + 1)
fig.update_yaxes(title_text="z-score across days", col=1)
fig.update_layout(title="Held-out check: groups found with the 6 earlier libraries, profile in the 3 batch 3 libraries",
                  legend=dict(orientation="h", y=1.2, x=0, xanchor="left", yanchor="bottom"))
fig = style(fig, 430)
fig.update_layout(margin=dict(t=140), title=dict(y=0.98, yanchor="top"))
parts.append(("08b", fig_html(fig, False)))

t8_cen = cen8[(cen8.k == K8) | (cen8.cluster == 0)].copy()
t8_cen.insert(0, "group", t8_cen.cluster.map(SC_LABEL))
t8_cen = t8_cen.drop(columns=["k", "cluster"])
FEAT = {"gained_step27": "Gained in step 27", "lost_step27": "Lost in step 27", "alpha_markers": "Alpha marker sites",
        "beta_markers": "Beta marker sites", **{f"KO_cluster{i}": f"KO-effect cluster {i}" for i in range(1, 7)},
        "KO_stays": "KO effect stays", "KO_other_KO": "Other KO effect", "KO_no_KO_effect": "No KO effect"}
l = lk8[lk8.k == K8].copy()
for c_ in ["odds_ratio_vs_stable", "p_vs_stable", "fdr_vs_stable"]:
    l[c_] = pd.to_numeric(l[c_].replace({"Inf": "inf"}), errors="coerce")
l.loc[l.feature.isin(["gained_step27", "lost_step27"]), ["odds_ratio_vs_stable", "p_vs_stable", "fdr_vs_stable"]] = np.nan
l["feature"] = l.feature.map(FEAT)
l["value"] = [f"{int(a):,} ({p:.1f}%)" + ("" if np.isnan(o) else f"; OR {o:.2g}, FDR {q:.2g}")
              for a, p, o, q in zip(l.n_with, l.pct_with, l.odds_ratio_vs_stable, l.fdr_vs_stable)]
t8_links = l.pivot(index="feature", columns="cluster", values="value").reindex(list(FEAT.values()))
t8_links.columns = [SC_LABEL[c] for c in t8_links.columns]
t8_links = t8_links.reset_index().rename(columns={"feature": "Sites"})
# gene level: genes linked to sites of each group vs genes whose linked sites are all stable
LINK8 = {"GREAT": "GREAT domain", "within_50kb": "TSS within ±50 kb"}
gl8 = rd(os.path.join(D8, "gene_links.tsv"))
gl8["value"] = [f"{int(n):,} ({p:.1f}%); OR {o:.2g}, adj. OR {ao:.2g}, adj. FDR {q:.2g}"
                for n, p, o, ao, q in zip(gl8.genes, gl8.pct_in_site_group, gl8.odds_ratio, gl8.adj_odds_ratio, gl8.adj_fdr)]
t8_genes = []
for ln, lab in LINK8.items():
    x = gl8[gl8.link == ln]
    pv = x.pivot(index="ko_group", columns="site_group", values="value").reindex([g for g in GROUPS if g in set(x.ko_group)])
    pv.columns = [f"Genes with a group {c} site" for c in pv.columns]
    st = x.groupby("ko_group").pct_stable_only.first()
    pv.insert(0, "Genes with only stable sites (%)", [f"{st[g]:.1f}%" for g in pv.index])
    pv.insert(0, "KO gene group", [LABEL[g] for g in pv.index])
    pv.insert(0, "Link", lab)
    t8_genes.append(pv.reset_index(drop=True))
t8_genes = pd.concat(t8_genes, ignore_index=True)
ng8 = {ln: {c: int(gl8[(gl8.link == ln) & (gl8.site_group == c)].genes_in_site_group.iloc[0]) for c in range(1, K8 + 1)}
       for ln in LINK8}

# gene list: KO clusters 2 and 4 genes with a varying site under either link
gs8 = rd(os.path.join(D8, "gene_site_groups.tsv"))
gs8["stable_only"] = gs8.stable_only.astype(str).str.upper() == "TRUE"
gv = gs8[~gs8.stable_only & gs8.group.isin(["cluster2", "cluster4"])]


def site_str(r):
    s = [f"{c}: {int(r[f'sites_group{c}'])}" for c in range(1, K8 + 1) if r[f"sites_group{c}"] > 0]
    return "; ".join(s) + f" (stable {int(r.sites_group0)})"


gv = gv.assign(sites=gv.apply(site_str, axis=1))
w = gv.pivot(index=["gene_id", "gene_name", "group", "marker"], columns="link", values="sites").reset_index()
for ln in LINK8:
    if ln not in w.columns:
        w[ln] = ""
w = w.fillna("").sort_values(["group", "gene_name"], kind="stable")
t8_list = pd.DataFrame({"Gene": w.gene_name, "KO gene group": w.group.map(LABEL), "Marker": w.marker,
                        "Varying sites by group, GREAT": w.GREAT, "Varying sites by group, ±50 kb": w.within_50kb})
n8_list = {g: int((w.group == g).sum()) for g in ["cluster2", "cluster4"]}

t8_ho = ho8.rename(columns={"pct_batch3_nearest_own": "% sites whose batch 3 profile is closest to own group",
                            "pct_expected_by_chance": "% expected by chance",
                            "pct_same_as_9lib_best_match": "% same group as the 9-library clustering"})

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
limma voom on the DiffBind counts of the 9 libraries (~ batch + day) and camera per group (competitive; main FDR with limma's default
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
Signal: DiffBind reads per site divided by the DiffBind background size factors (step 00), mean of the 3 libraries per day; boxes show
the median and quartiles, whiskers 1.5 x IQR, outliers not drawn. Stars: paired Wilcoxon signed-rank test between days
within the group (BH; *** FDR &lt; 0.001, ** &lt; 0.01, * &lt; 0.05). Binding rises after D15 at most sites, so these
tests show that binding changes at a group's sites, not that it changes more than at other sites; that comparison is
the camera test in section 03 (cluster 4: up D18 vs D15, FDR 0.21; down D22 vs D18, FDR 0.018; both not
significant with the estimated inter-site correlation).</p>
{parts[2][1]}
{table_html(t05)}

<h2>06 Alpha-leaning against beta-leaning sites, per day and over time</h2>
<p class="note">Marker genes (human islet cell-type markers): alpha {", ".join(MARKERS["alpha"])}; beta
{", ".join(MARKERS["beta"])}. Cluster 4 holds {n_alpha4} alpha and {n_beta4} beta markers (with SST, PYY and GAST);
cluster 2 holds {n_beta2} beta and {n_alpha2} alpha markers. Site sets (sites linked by GREAT): cluster 4
({set_n["cluster4"]:,}) and cluster 2 ({set_n["cluster2"]:,}); all alpha ({set_n["alpha"]}) and all beta
({set_n["beta"]}) marker sites; alpha markers in cluster 4 ({set_n["alpha_in_c4"]}) and beta markers in cluster 2
({set_n["beta_in_c2"]}). For each of the 9 libraries: signal per site = DiffBind reads / background size factor (step 00),
summarized per set as the total (depends on the number of sites) or the median per site (does not). log2 ratio
of the two sets fitted on ~ batch + day with the libraries as replicates (5 residual df): level_&lt;day&gt; tests whether
one set binds more on that day (between different site sets this also reflects site number and intrinsic site
strength), change_&lt;a&gt;_vs_&lt;b&gt; whether one set gains or loses relative to the other between days. p values are not
corrected for the number of tests. Figure: signal per site (mean of the 3 libraries per day) at the alpha and the
beta marker sites, and at the cluster 4 and the cluster 2 sites; boxes and whiskers as in section 05; stars and p
from the library-level test of the median per site of the same two sets (level_&lt;day&gt;; *** &lt; 0.001,
** &lt; 0.01, * &lt; 0.05, ns otherwise).</p>
{parts[3][1]}
{parts[4][1]}
{table_html(ct_show)}
<p class="note">Marker genes with at least one linked site: ISL1 binding (sum of the normalized signal of its sites)
and mean normalized expression in control islets (Ctr-H3), both as log2 change from D15.</p>
{parts[5][1]}
<p class="note">All marker genes: KO-effect group, linked sites, ISL1 binding (sum over sites), control expression
(mean normalized counts) and the KO effect (shrunk log2FC KO vs Ctr) per day; candidate = in the section 04 list.</p>
{table_html(mk_show)}

<h2>07 Target genes of the differentially bound sites and their KO effect</h2>
<p class="note">Differentially bound sites: step 27 DiffBind, background normalization, FDR &lt; 0.05 (gained: up in D18 vs
D15 or D22 vs D15; lost: down; no site changes between D18 and D22). Target genes: GREAT, bulk-tested genes only.
Gene classes: gained site, lost site (no gained site), only stable sites, no site. KO effect per day: KO vs Ctr padj
&lt; 0.05, lower in the KO (ISL1 activates) or higher (ISL1 represses). Tests compare genes with a gained site to
genes with only stable sites: Fisher's exact test, and logistic regression with domain length and number of sites
(adjusted OR, BH over the tests), since larger domains with more sites are more likely both to hold a gained site
and to have a KO effect.</p>
{parts[6][1]}
{table_html(s07_show)}
{table_html(t07)}
<p class="note">Genes with a gained site and a KO effect on any day (LRT_genotype padj &lt; 0.05), one row per gene with
its most significant gained site, ordered by the D22 KO padj. Every differential site and target gene:
results/07_diff_site_targets/diff_site_genes.tsv.</p>
{table_html(g07)}
<p class="note">Binding signal (normalized reads per site, mean of the 3 libraries per day) at the gained sites of these
genes, all together and split by the direction of the gene's KO effect on the latest day it is significant, with all
14,431 sites as the reference. Boxes and whiskers as in section 05. The sites were selected for gaining binding
against D15, so their rise from D15 is expected. D22 vs D18 is not part of the selection, but sites picked for a
high D18 reading also tend to read lower at the next measurement (regression to the mean), so part of their fall
to D22 is expected as well. 11 genes have a KO effect only in the LRT, on no single day, and are in the first box
only. Stars: paired Wilcoxon
between days within the set (sites as units, uncorrected; descriptive).</p>
{parts[7][1]}
{table_html(w07)}

<h2>08 All ISL1 sites grouped by binding pattern across days</h2>
<p class="note">All 14,431 sites, step 27 DiffBind counts and background size factors. DESeq2 LRT (~ batch + day against
~ batch) picks the sites whose binding varies across days (p &lt; 0.05, a filter against noise, not a significance call);
the rest are the stable group. Per site: DESeq2 VST, batch 3 effect removed, mean per day, z-score across the 3 days;
k-means on the varying sites (k = 4, from the elbow: the share of variance between groups is 0.71 / 0.81 / 0.87 / 0.91
for k = 3-6; more groups split the same shapes). Lines: median z per day; band: 25th-75th percentile.</p>
{dict(parts)["08a"]}
<p class="note"><b>Caveat: only the rise after D15 replicates in an independent batch.</b> The clustering repeated with the 6
earlier libraries finds the same shapes, but the 3 batch 3 libraries (D15 rep1, D18 rep1, D22 rep2) do not follow them
between D18 and D22: every group rises from D15 in batch 3, while its D18-vs-D22 shape is not reproduced (the
"peaks at D18" group rises further to D22, the "late" group peaks at D18, the "falls at D22" group rises). Single
sites match their own group no more often than chance. The earlier libraries come from another differentiation with
relabelled days (D12 as D15, D15 as D18), which may carry differences between D18 and D22 that batch 3 does not
show. The groups and their links below describe the 9-library data; the distinction between groups 1-4 beyond the
rise after D15 awaits the unified CUT&amp;RUN time course.</p>
{dict(parts)["08b"]}
{table_html(t8_ho)}
<p class="note">Groups: number of sites, z-score centre per day and median VST per day (batch-corrected).</p>
{table_html(t8_cen)}
<p class="note">Sites per group with each feature: number (share of the group); for the gene links (GREAT, step 01), odds
ratio and BH FDR of a Fisher test of the group against the stable sites. Gained and lost sites are part of what the
LRT selects, so they have no test.</p>
{table_html(t8_links)}
<h3>Genes linked to each group and their KO effect</h3>
<p class="note">Gene level, two ways of linking sites to bulk-tested genes: GREAT domains (step 01; main) and, as a check,
every gene whose TSS lies within ±50 kb of the site. A gene counts in every group in which it has at least one linked
site (GREAT: {ng8["GREAT"][1]:,} / {ng8["GREAT"][2]:,} / {ng8["GREAT"][3]:,} / {ng8["GREAT"][4]:,} genes in groups 1-4;
±50 kb: {ng8["within_50kb"][1]:,} / {ng8["within_50kb"][2]:,} / {ng8["within_50kb"][3]:,} / {ng8["within_50kb"][4]:,})
and is compared with the genes whose linked sites are all stable. Genes with a varying site have more sites (median 2
against 1), which alone raises their chance of a KO effect, so the main figure is the adjusted odds ratio: logistic
regression of KO group membership on group membership with log2(1 + number of sites) and, for GREAT, log10 domain
length as covariates; BH FDR over all rows of each link. Each cell: genes in the KO group (share of the site group's
genes); Fisher odds ratio; adjusted odds ratio and its FDR.</p>
<p class="note">KO cluster 4 (down, emerging and growing) is enriched near the varying sites after adjustment: groups 1
and 2 with GREAT (adj. OR 2.2 and 2.0), groups 2 and 4 with ±50 kb (2.9 and 4.0). KO cluster 2 (up early) is enriched
only before adjustment with GREAT, and after adjustment only for group 2 with ±50 kb (1.8). Group 3 (late rise) links
to no KO group under either link.</p>
{table_html(t8_genes)}
<details><summary>Gene list: {n8_list["cluster2"] + n8_list["cluster4"]} genes of KO clusters 2 ({n8_list["cluster2"]}) and 4
({n8_list["cluster4"]}) with a varying site under either link</summary>
<p class="note">Varying sites per group (group: number of sites), with the number of stable sites in brackets; blank: no
varying site under that link.</p>
{table_html(t8_list)}
</details>
</body></html>"""
with open(OUT, "w") as f:
    f.write(doc)
print("wrote", OUT)
