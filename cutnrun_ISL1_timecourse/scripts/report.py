#!/usr/bin/env python3
"""Step 09: HTML report of the ISL1 time course: library QC (steps 01-04), consensus peaks (05),
differential binding (06-07). Standard library only; plots are inline SVG.
  python report.py <results> <sample_metadata.tsv> <out.html>
"""
import csv
import datetime
import html
import math
import os
import sys

OUT, META, HTML_OUT = sys.argv[1:4]
D = os.path.join(OUT, "diff")
SCRATCH = os.path.dirname(OUT)
TH = 0.05
CONTRASTS = ["D18_vs_D15", "D22_vs_D18", "D22_vs_D15"]
NORMS = ["background", "rle"]
DAY_COL = {"D15": "#4c78a8", "D18": "#f58518", "D22": "#54a24b"}


def tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x, dec=0):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "–"
    return "–" if v != v else f"{v:,.{dec}f}"


def table(head, rows, num_cols=()):
    h = "".join(f'<th class="{"n" if i in num_cols else ""}">{html.escape(c)}</th>' for i, c in enumerate(head))
    b = "".join("<tr>" + "".join(f'<td class="{"n" if i in num_cols else ""}">{c}</td>' for i, c in enumerate(r)) + "</tr>" for r in rows)
    return f'<div class="tw"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def scale(lo, hi, a, b):
    return lambda v: a + (v - lo) / (hi - lo) * (b - a) if hi > lo else (a + b) / 2


def ma_svg(rows, title, w=360, h=260):
    pts = []
    for r in rows:
        try:
            x, y, q = float(r["Conc"]), float(r["log2FC"]), float(r["FDR"])
        except ValueError:
            continue
        pts.append((x, y, q < TH))
    if not pts:
        return ""
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    lim = max(2.0, max(abs(min(ys)), abs(max(ys))))
    sx = scale(min(xs), max(xs), 40, w - 10); sy = scale(-lim, lim, h - 30, 22)
    ns = [p for p in pts if not p[2]]
    step = max(1, len(ns) // 4000)
    out = [f'<svg viewBox="0 0 {w} {h}" class="plot" role="img" aria-label="{html.escape(title)}">',
           f'<text x="{w/2}" y="14" class="pt">{html.escape(title)}</text>',
           f'<line x1="40" x2="{w-10}" y1="{sy(0):.1f}" y2="{sy(0):.1f}" class="ax"/>']
    for x, y, _ in ns[::step]:
        out.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="1.2" class="ns"/>')
    for x, y, s in pts:
        if s:
            out.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="1.6" class="{"up" if y > 0 else "dn"}"/>')
    out.append(f'<text x="{w/2}" y="{h-6}" class="al">mean log2 normalized reads</text>')
    out.append(f'<text x="10" y="{h/2}" class="al" transform="rotate(-90 10 {h/2})">log2FC</text>')
    for v in (-lim, lim):
        out.append(f'<text x="36" y="{sy(v)+4:.1f}" class="tk" text-anchor="end">{v:.0f}</text>')
    out.append("</svg>")
    return "".join(out)


def pca_svg(rows, var, title, w=360, h=280):
    xs = [float(r["PC1"]) for r in rows]; ys = [float(r["PC2"]) for r in rows]
    pad = 0.15 * max(max(xs) - min(xs), max(ys) - min(ys), 1)
    sx = scale(min(xs) - pad, max(xs) + pad, 40, w - 10); sy = scale(min(ys) - pad, max(ys) + pad, h - 30, 22)
    out = [f'<svg viewBox="0 0 {w} {h}" class="plot" role="img" aria-label="{html.escape(title)}">',
           f'<text x="{w/2}" y="14" class="pt">{html.escape(title)}</text>']
    for r, x, y in zip(rows, xs, ys):
        c = DAY_COL[r["day"]]
        if r["batch"] == "batch3":
            out.append(f'<rect x="{sx(x)-5:.1f}" y="{sy(y)-5:.1f}" width="10" height="10" fill="{c}"/>')
        else:
            out.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="5.5" fill="{c}"/>')
        out.append(f'<text x="{sx(x)+8:.1f}" y="{sy(y)+4:.1f}" class="tk">{html.escape(r["sample"].replace("_ISL1", ""))}</text>')
    out.append(f'<text x="{w/2}" y="{h-6}" class="al">PC1 ({var["PC1"]}%)</text>')
    out.append(f'<text x="10" y="{h/2}" class="al" transform="rotate(-90 10 {h/2})">PC2 ({var["PC2"]}%)</text>')
    out.append("</svg>")
    return "".join(out)


def read_text(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return ""


def grab(text, key):
    for line in text.splitlines():
        if line.startswith(key):
            try:
                return float(line.split(":", 1)[1].split()[0].rstrip("%"))
            except (ValueError, IndexError):
                return None
    return None


S = []
meta = tsv(META)
rows = []
for m in meta:
    s = m["sample_name"]
    q = os.path.join(OUT, "qc", "alignment", s)
    bt = read_text(q + ".bowtie2.log")
    aln = next((l.split("%")[0] for l in bt.splitlines() if "overall alignment rate" in l), None)
    md = read_text(q + ".markdup.metrics.txt")
    dup_t, exam = grab(md, "DUPLICATE TOTAL"), grab(md, "EXAMINED")
    tx = os.path.join(SCRATCH, "tenx_filter", s + ".tsv")
    tenx = tsv(tx)[0]["pct_removed"] + "%" if os.path.exists(tx) else "not filtered"
    ps = os.path.join(OUT, "qc", "peaks", s + ".peak_stats.tsv")
    p = tsv(ps)[0] if os.path.exists(ps) else {}
    rows.append([html.escape(s), m["day"], m["replicate"], m["batch"], html.escape(m["original_sample"]), tenx,
                 (aln + "%") if aln else "–", f"{100 * dup_t / exam:.1f}%" if dup_t is not None and exam else "–",
                 num(int(p["fragments"]) / 1e6, 1) + " M" if p else "–", num(p.get("regions")),
                 f'{100 * float(p["frip"]):.1f}%' if p else "–"])
S.append(("Libraries", "<p>10x removed: share of read pairs removed by step 01 (batch 3 delivery only). "
          "Align: bowtie2 overall alignment rate. Dup: duplicates marked by samtools markdup (kept for ISL1, "
          "removed from the IgG used as background). Fragments, peaks (distinct regions) and FRiP: step 04, "
          "MACS3 q 0.01 against the library's own IgG.</p>" +
          table(["Library", "Day", "Rep", "Batch", "Original library", "10x removed", "Align", "Dup", "Fragments", "Peaks", "FRiP"],
                rows, num_cols=(5, 6, 7, 8, 9, 10))))

cs = tsv(os.path.join(OUT, "consensus", "summary.tsv"))
S.append(("Consensus peaks per day (MACS3, >= 2 of 3)", "<p>A region is kept when it holds MACS3 peaks from "
          "at least 2 of the day's 3 libraries (macs3_ge2, the primary set); 3 of 3 is a high-confidence "
          "subset. Chance: each library's regions shuffled independently and counted alike. "
          "Batch 3 + earlier: supported by the batch 3 library and an earlier library.</p>" +
          table(["Day", "Set", "Regions", "Chance", "× chance", "Mb", "Median bp", "Batch 3 + earlier", "Earlier only"],
                [[c["day"], c["set"], f'<b>{num(c["regions"])}</b>', num(c["chance"]), c["over_chance"], c["Mb"],
                  num(c["median_bp"]), num(c["batch3_and_earlier"]), num(c["earlier_only"])] for c in cs],
                num_cols=(2, 3, 4, 5, 6, 7, 8))))

samples = tsv(os.path.join(D, "samples.tsv"))
nf = tsv(os.path.join(D, "norm_factors.tsv"))
S.append(("Differential binding: design", "<p>Control ISL1 CUT&amp;RUN only, three replicates per day (one batch 3 library, two earlier "
          "libraries). Model <code>~ batch + day</code>, batch = batch 3 vs earlier. Contrasts: D18 vs D15, D22 vs D18, "
          "D22 vs D15; log2FC &gt; 0 = higher on the later day. FDR &lt; 0.05. Main analysis DiffBind (DESeq2); "
          "check: DESeq2 on independent featureCounts counts. Two normalizations, side by side: <b>background</b> "
          "(15 kb genome bins) and <b>rle</b> (median of ratios on the sites).</p>" +
          table(["Sample", "Day", "Rep", "Batch", "Original library"],
                [[html.escape(s["sample"]), s["day"], s["replicate"], s["batch"], html.escape(s["original_sample"])] for s in samples])))

qc = tsv(os.path.join(D, "sites", "qc.tsv")) + tsv(os.path.join(D, "greylist_summary.tsv"))
S.append(("Differential binding: site set (401 bp)", "<p>macs3_ge2 regions of D15, D18 and D22, each centred on its strongest summit among "
          "that day's trio libraries; summit ± 200 bp; overlapping sites resolved by the stronger summit; blacklist "
          "and an IgG greylist (GreyListChIP on the 9 matched IgGs) removed. \"Called on\" = the days whose "
          "macs3_ge2 set overlaps the site.</p>" +
          table(["Item", "n"], [[html.escape(r["item"]), num(r["n"], 2 if "Mb" in r["item"] else 0)] for r in qc], num_cols=(1,))))

rows = []
by = {(r["method"], r["normalization"], r["sample"]): r for r in nf}
for s in samples:
    rows.append([html.escape(s["sample"])] + [num(by[(m, n, s["sample"])]["norm_factor"], 3) for m in ("diffbind", "deseq2") for n in NORMS])
S.append(("Normalization factors", "<p>DESeq2 size factors as used by each method. If background and rle differ by "
          "day, the two normalizations disagree about a global change in binding.</p>" +
          table(["Sample", "DiffBind background", "DiffBind rle", "DESeq2 background", "DESeq2 rle"], rows, num_cols=(1, 2, 3, 4))))

pca = tsv(os.path.join(D, "pca.tsv")); pv = tsv(os.path.join(D, "pca_variance.tsv"))
figs = []
for t, lab in (("vst", "PCA, vst"), ("vst_batch_removed", "PCA, batch term removed")):
    var = {r["PC"]: r["pct"] for r in pv if r["type"] == t}
    figs.append(pca_svg([r for r in pca if r["type"] == t], var, lab))
S.append(("Libraries (PCA of the top 5,000 variable sites)", "<p>Colour = day (D15 blue, D18 orange, D22 green); "
          "square = batch 3, circle = earlier library.</p>" + f'<div class="figs">{"".join(figs)}</div>'))

summ = tsv(os.path.join(D, "summary.tsv"))
rows = [[r["contrast"].replace("_", " "), r["normalization"], "DiffBind" if r["method"] == "diffbind" else "DESeq2",
         num(r["tested"]), f'<b>{num(int(r["up"]) + int(r["down"]))}</b>', num(r["up"]), num(r["down"]), num(r["up_fc1"]), num(r["down_fc1"])]
        for r in sorted(summ, key=lambda r: (CONTRASTS.index(r["contrast"]), NORMS.index(r["normalization"]), r["method"] != "diffbind"))]
S.append(("Differential sites (FDR < 0.05)", "<p>up = higher on the later day.</p>" +
          table(["Contrast", "Normalization", "Method", "Tested", "Significant", "Up", "Down", "Up, log2FC ≥ 1", "Down, log2FC ≤ −1"],
                rows, num_cols=(3, 4, 5, 6, 7, 8))))

agr = tsv(os.path.join(D, "agreement.tsv"))
rows = [[r["contrast"].replace("_", " "), r["a"], r["b"], num(r["sig_a"]), num(r["sig_b"]), num(r["both_same_direction"]), num(r["opposite"])] for r in agr]
S.append(("Agreement between methods and normalizations", "<p>Same sites, called significant by both, in the same "
          "direction.</p>" + table(["Contrast", "A", "B", "Significant in A", "Significant in B", "Both, same direction", "Opposite direction"],
                                   rows, num_cols=(3, 4, 5, 6))))

figs = []
for n in NORMS:
    for c in CONTRASTS:
        figs.append(ma_svg(tsv(os.path.join(D, "tables", f"diffbind_{n}_{c}.tsv")), f"DiffBind {n}: {c.replace('_', ' ')}"))
S.append(("MA plots (DiffBind)", "<p>Orange = up, blue = down (FDR &lt; 0.05); grey = not significant (thinned for size).</p>"
          f'<div class="figs">{"".join(figs)}</div>'))

sites = {}
with open(os.path.join(D, "sites_used.bed")) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        sites[p[3]] = p
for n in NORMS:
    parts = []
    for c in CONTRASTS:
        t = [r for r in tsv(os.path.join(D, "tables", f"diffbind_{n}_{c}.tsv")) if float(r["FDR"]) < TH]
        if not t:
            parts.append(f"<h3>{c.replace('_', ' ')}</h3><p>No significant sites.</p>")
            continue
        top = sorted(t, key=lambda r: (float(r["FDR"]), -abs(float(r["log2FC"]))))[:15]
        parts.append(f"<h3>{c.replace('_', ' ')}</h3>" + table(
            ["Site", "Location", "Called on", "Mean", "log2FC", "FDR"],
            [[r["name"], f'{r["seqnames"]}:{int(r["start"]):,}-{int(r["end"]):,}', sites[r["name"]][7], num(r["Conc"], 2),
              num(r["log2FC"], 2), f'{float(r["FDR"]):.1e}'] for r in top], num_cols=(3, 4, 5)))
    S.append((f"Top sites, DiffBind {n}", "".join(parts)))

S.append(("Files", "<p>Under <code>results/</code> (scratch; keepers copied to the folder's <code>results/</code> on "
          "/pl/active): <code>consensus/&lt;day&gt;_Ctr_ISL1.{macs3_ge2,macs3_3of3}.bed</code>, "
          "<code>tracks/</code> (bigWigs per library and per day, beds, <code>igv_session.xml</code>), and in "
          "<code>diff/</code>: <code>sites_used.bed</code>, "
          "<code>tables/{diffbind,deseq2}_&lt;norm&gt;_&lt;contrast&gt;.tsv</code> (all sites), "
          "<code>beds/diffbind_&lt;norm&gt;_&lt;contrast&gt;_{up,down}.bed</code>, <code>summary.tsv</code>, "
          "<code>agreement.tsv</code>, <code>norm_factors.tsv</code>, <code>counts_sites.tsv</code>.</p>"))

css = """
:root{--bg:#fff;--fg:#1d1d1f;--mut:#5f6368;--line:#dcdcdc;--head:#f3f4f6;--acc:#1f5fa8;--ns:#b8bcc2;--up:#e8710a;--dn:#1a73e8}
@media (prefers-color-scheme:dark){:root{--bg:#16181b;--fg:#e8e8e8;--mut:#a0a4aa;--line:#33363b;--head:#22252a;--acc:#7fb0ea;--ns:#555a61;--up:#f59e4c;--dn:#6ea8fe}}
body{background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;margin:0 auto;max-width:1200px;padding:24px 16px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 8px;color:var(--acc)}h3{font-size:14px;margin:16px 0 4px}
p{margin:6px 0 10px}.sub{color:var(--mut);margin-bottom:16px}.tw{overflow-x:auto}
table{border-collapse:collapse;margin:8px 0 4px;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:5px 10px;text-align:left;white-space:nowrap}
th{background:var(--head);font-weight:600}td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
code{font-size:12px}.figs{display:flex;flex-wrap:wrap;gap:12px}.plot{width:360px;max-width:100%;height:auto}
.plot .pt{font-size:12px;text-anchor:middle;fill:var(--fg)}.plot .al{font-size:10px;text-anchor:middle;fill:var(--mut)}
.plot .tk{font-size:9px;fill:var(--mut)}.plot .ax{stroke:var(--line)}.ns{fill:var(--ns)}.up{fill:var(--up)}.dn{fill:var(--dn)}
"""
body = "".join(f"<h2>{html.escape(t)}</h2>{c}" for t, c in S)
os.makedirs(os.path.dirname(HTML_OUT), exist_ok=True)
with open(HTML_OUT, "w") as f:
    f.write(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>ISL1 CUT&amp;RUN time course</title>
<style>{css}</style></head><body>
<h1>ISL1 CUT&amp;RUN time course</h1>
<div class="sub">Control ISL1, D15 / D18 / D22, three libraries per day (one batch 3, two earlier) · built {datetime.date.today().isoformat()}</div>
{body}</body></html>""")
print(f"wrote {HTML_OUT}")
