#!/usr/bin/env python3
"""Plot genes of choice the way the report plots the islet marker genes, with significance stars.

Per gene: DESeq2 normalized count + 1 (log scale) per sample, Ctr and KO, D15 / D18 / D22, lines through the
group means (same style as the "Islet marker genes" section of build_report.py).
Stars above each day: KO vs Ctr at that day (KO_vs_Ctr_<day>.tsv, Wald padj).
Stars in the panel title: time-dependent KO effect (LRT_interaction.tsv, padj).
Stars: * padj < 0.05, ** < 0.01, *** < 0.001; ns = not significant; NA = padj not computed by DESeq2.

Run on the Mac after fetching results/DE_DESeq2 from Alpine:
  ~/anaconda3/bin/python3 plot_genes.py INS [GENE ...] [--out NAME]
Output: results/gene_plots/<NAME>.html (default NAME = the gene names joined by "_")
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

HERE = os.path.dirname(os.path.abspath(__file__))
DE = os.path.join(HERE, "results", "DE_DESeq2")
OUT_DIR = os.path.join(HERE, "results", "gene_plots")

DAYS = ["D15", "D18", "D22"]
GENOS = ["Ctr", "KO"]
C = {"Ctr": "#2a78d6", "KO": "#eb6834", "text": "#0b0b0b", "text2": "#52514e", "muted": "#898781",
     "grid": "#e6e5e0", "surface": "#fcfcfb"}
NCOL = 5


def read(name):
    return pd.read_csv(os.path.join(DE, name), sep="\t", keep_default_na=False, na_values=["NA"])


def stars(p):
    if pd.isna(p):
        return "NA"
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def short(sample):
    p = sample.split("_")
    g = "Ctr" if p[2].startswith("Ctr") else "KO"
    return f"{g} {p[3]} r{p[4][-1]}"


def pick(df, gene):
    """The row of a gene by name; a name shared by several gene ids takes the first (stable) and warns."""
    hit = df[df["gene_name"] == gene]
    if len(hit) > 1:
        print(f"note: {gene} matches {len(hit)} gene ids; using {hit.iloc[0]['gene_id']}", file=sys.stderr)
    return hit.iloc[0] if len(hit) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes", nargs="+")
    ap.add_argument("--out", default=None, help="output name (default: genes joined by _)")
    a = ap.parse_args()

    meta = pd.read_csv(os.path.join(HERE, "sample_sheet.csv"))
    norm = read("normalized_counts.tsv")
    contr = {d: read(f"KO_vs_Ctr_{d}.tsv") for d in DAYS}
    lrt = read("LRT_interaction.tsv")

    genes = [g for g in a.genes if pick(norm, g) is not None]
    for g in set(a.genes) - set(genes):
        print(f"skipped: {g} is not among the tested genes", file=sys.stderr)
    if not genes:
        sys.exit("no gene to plot")

    titles = []
    for gene in genes:
        r = pick(lrt, gene)
        p = np.nan if r is None else r["padj"]
        titles.append(f"<b>{gene}</b><br><span style='font-size:11px'>KO effect over time: {stars(p)}</span>")

    ncol = min(NCOL, len(genes))
    nrow = int(np.ceil(len(genes) / ncol))
    fig = make_subplots(rows=nrow, cols=ncol, subplot_titles=titles, vertical_spacing=0.22 / max(nrow - 1, 1) + 0.08,
                        horizontal_spacing=0.06)
    rows = []
    for i, gene in enumerate(genes):
        r_, c_ = i // ncol + 1, i % ncol + 1
        row = pick(norm, gene)
        top = 0.0
        for g in GENOS:
            xs, ys, hs = [], [], []
            for d in DAYS:
                ss = meta[(meta["genotype"] == g) & (meta["day"] == d)].sort_values("replicate", kind="stable")["sample"]
                for smp in ss:
                    xs.append(d); ys.append(float(row[smp]) + 1); hs.append(short(smp))
            top = max(top, max(ys))
            means = [np.mean([ys[k] for k in range(len(xs)) if xs[k] == d]) for d in DAYS]
            fig.add_trace(go.Scatter(x=DAYS, y=means, mode="lines", line=dict(color=C[g], width=2), name=g,
                                     legendgroup=g, showlegend=(i == 0), hoverinfo="skip"), row=r_, col=c_)
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", legendgroup=g, showlegend=False, customdata=hs,
                                     marker=dict(color=C[g], size=8, line=dict(color=C["surface"], width=2)),
                                     hovertemplate=f"{gene}<br>%{{customdata}}<br>%{{y:.0f}} (normalized count + 1)"
                                                   "<extra></extra>"), row=r_, col=c_)
        # per-day stars, above the day's highest point
        lo = min(float(row[s]) + 1 for s in meta["sample"])
        ymax = np.log10(top) + 0.12 * max(np.log10(top) - np.log10(lo), 0.5)
        rec = {"gene": gene}
        for d in DAYS:
            h = pick(contr[d], gene)
            p = np.nan if h is None else h["padj"]
            dmax = max(float(row[s]) + 1 for s in meta.loc[meta["day"] == d, "sample"])
            s = stars(p)
            fig.add_annotation(x=d, y=np.log10(dmax), yshift=14, text=s, showarrow=False, row=r_, col=c_,
                               font=dict(size=14 if "*" in s else 11, color=C["text"] if "*" in s else C["muted"]))
            rec[f"log2FC_{d}"] = None if h is None else round(h["log2FoldChange"], 3)
            rec[f"padj_{d}"] = p
        rl = pick(lrt, gene)
        rec["padj_LRT_interaction"] = np.nan if rl is None else rl["padj"]
        rows.append(rec)
        fig.update_yaxes(type="log", range=[np.log10(lo) - 0.1 * max(np.log10(top) - np.log10(lo), 0.5), ymax],
                         row=r_, col=c_)
    for r_ in range(1, nrow + 1):
        fig.update_yaxes(title="normalized count + 1", row=r_, col=1)
    for an in fig.layout.annotations[:len(genes)]:
        an.font = dict(size=13, color=C["text"])
    w = max(260 * ncol + 150, 470)
    fig.update_layout(
        width=w, height=320 * nrow + 170, margin=dict(l=70, r=80, t=80, b=120),
        paper_bgcolor=C["surface"], plot_bgcolor=C["surface"],
        font=dict(family="Inter, -apple-system, Segoe UI, Helvetica, Arial, sans-serif", size=12, color=C["text"]),
        legend=dict(orientation="v", x=1.02, xanchor="left", y=1, yanchor="top"), hoverlabel=dict(bgcolor="white"))
    fig.add_annotation(text="Above each day: KO vs Ctr at that day (DESeq2 Wald padj)<br>"
                            "Title: KO effect over time (LRT_interaction padj)<br>"
                            "* &lt; 0.05, ** &lt; 0.01, *** &lt; 0.001, ns not significant",
                       xref="paper", yref="paper", x=0, y=0, yshift=-45, yanchor="top", showarrow=False,
                       xanchor="left", align="left", font=dict(size=10, color=C["text2"]))
    fig.update_xaxes(gridcolor=C["grid"], linecolor=C["text2"], ticks="outside", tickcolor=C["grid"])
    fig.update_yaxes(gridcolor=C["grid"], linecolor=C["text2"], ticks="outside", tickcolor=C["grid"],
                     title_font=dict(color=C["text2"]))

    os.makedirs(OUT_DIR, exist_ok=True)
    name = a.out or "_".join(genes)
    out = os.path.join(OUT_DIR, f"{name}.html")
    fig.write_html(out, include_plotlyjs=True, config={"displaylogo": False,
                   "toImageButtonOptions": {"format": "svg", "filename": name}})
    tbl = pd.DataFrame(rows)
    tbl.to_csv(os.path.join(OUT_DIR, f"{name}.tsv"), sep="\t", index=False, na_rep="NA")
    print(tbl.to_string(index=False))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
