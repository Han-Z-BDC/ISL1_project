#!/usr/bin/env python3
"""Step 09: HTML report of the ISL1 time course: library QC (steps 01-04), consensus peaks (05),
differential binding (06-07), HOMER motifs (10-11), peak-call overlap (12). Standard library only; plots are inline SVG.
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

PO = os.path.join(OUT, "peak_overlap")
if os.path.isdir(PO):
    DAYS = ["D15", "D18", "D22"]
    cls = tsv(os.path.join(PO, "classes.tsv"))
    lab = {r["class"]: ("all three days" if r["class"] == "D15,D18,D22" else
                        r["class"].replace(",", " + ") + ("" if "," in r["class"] else " only")) for r in cls}
    nsite = sum(int(r["sites"]) for r in cls)

    def upset_svg(w=560, h=300):
        top, mid, bot, left = 30, 200, 270, 50
        mx = max(int(r["sites"]) for r in cls) * 1.12
        bw = (w - left - 10) / len(cls)
        sy = scale(0, mx, mid, top)
        out = [f'<svg viewBox="0 0 {w} {h}" class="plot wide" role="img" aria-label="Overlap classes">',
               f'<text x="{w/2}" y="14" class="pt">Sites by the days whose consensus peaks contain them ({nsite:,})</text>',
               f'<line x1="{left}" x2="{w-10}" y1="{mid}" y2="{mid}" class="ax"/>']
        for i, r in enumerate(cls):
            x = left + i * bw + bw * 0.2; n = int(r["sites"])
            out.append(f'<rect x="{x:.1f}" y="{sy(n):.1f}" width="{bw*0.6:.1f}" height="{mid-sy(n):.1f}" rx="3" class="bar">'
                       f'<title>{html.escape(lab[r["class"]])}: {n:,} sites ({r["pct_sites"]}%)</title></rect>')
            out.append(f'<text x="{x+bw*0.3:.1f}" y="{sy(n)-4:.1f}" class="tk" text-anchor="middle">{n:,}</text>')
            on = [j for j, d in enumerate(DAYS) if d in r["class"].split(",")]
            cy = lambda j: mid + 18 + j * 22
            if len(on) > 1:
                out.append(f'<line x1="{x+bw*0.3:.1f}" x2="{x+bw*0.3:.1f}" y1="{cy(on[0])}" y2="{cy(on[-1])}" class="lnk"/>')
            for j in range(3):
                out.append(f'<circle cx="{x+bw*0.3:.1f}" cy="{cy(j)}" r="6" class="{"don" if j in on else "doff"}"/>')
        for j, d in enumerate(DAYS):
            out.append(f'<text x="{left-8}" y="{mid+22+j*22}" class="tk" text-anchor="end">{d}</text>')
        out.append(f'<text x="12" y="{(top+mid)/2}" class="al" transform="rotate(-90 12 {(top+mid)/2})">sites</text>')
        out.append("</svg>")
        return "".join(out)

    libs = tsv(os.path.join(PO, "libraries.tsv"))

    def libs_svg(w=560, h=260):
        top, base, left = 30, 200, 50
        bw = (w - left - 10) / len(libs)
        sy = scale(0, 100, base, top)
        out = [f'<svg viewBox="0 0 {w} {h}" class="plot wide" role="img" aria-label="Sites called per library">',
               f'<text x="{w/2}" y="14" class="pt">Sites called by each library\'s own MACS3 peaks (% of {nsite:,})</text>',
               f'<line x1="{left}" x2="{w-10}" y1="{base}" y2="{base}" class="ax"/>']
        for v in (0, 50, 100):
            out.append(f'<text x="{left-6}" y="{sy(v)+3:.1f}" class="tk" text-anchor="end">{v}</text>')
        for i, r in enumerate(libs):
            x = left + i * bw + bw * 0.2; p = float(r["pct_sites_called"])
            b3 = r["batch"] == "batch3"
            out.append(f'<rect x="{x:.1f}" y="{sy(p):.1f}" width="{bw*0.6:.1f}" height="{base-sy(p):.1f}" rx="3" '
                       f'fill="{DAY_COL[r["day"]]}" fill-opacity="{0.45 if b3 else 1}"><title>{html.escape(r["sample"])}: '
                       f'{p:.1f}% of sites, {int(r["peaks"]):,} peaks, {float(r["pct_reads_in_sites"]):.2f}% of reads in sites'
                       f'</title></rect>')
            out.append(f'<text x="{x+bw*0.3:.1f}" y="{sy(p)-4:.1f}" class="tk" text-anchor="middle">{p:.0f}</text>')
            l = r["sample"].replace("_Ctr_", " ").replace("_ISL1", "")
            out.append(f'<text x="{x+bw*0.3:.1f}" y="{base+14}" class="tk" text-anchor="middle">{html.escape(l)}</text>')
            out.append(f'<text x="{x+bw*0.3:.1f}" y="{base+26}" class="tk" text-anchor="middle">{"batch 3" if b3 else "earlier"}</text>')
        out.append("</svg>")
        return "".join(out)

    pc = tsv(os.path.join(PO, "presence_vs_counts.tsv"))

    def presence_svg(ct, w=360, h=280):
        rows = [r for r in pc if r["contrast"] == ct]
        b, a = ct.split("_vs_")
        vals = [float(r[k]) for r in rows for k in ("p10_log2FC", "p90_log2FC")]
        lo, hi = min(vals + [0]) - 0.1, max(vals + [0]) + 0.1
        top, base, left = 26, h - 48, 44
        bw = (w - left - 10) / len(rows)
        sy = scale(lo, hi, base, top)
        out = [f'<svg viewBox="0 0 {w} {h}" class="plot" role="img" aria-label="log2FC by peak presence {b} vs {a}">',
               f'<text x="{w/2}" y="14" class="pt">{b} vs {a}: DiffBind log2FC by peak presence</text>',
               f'<line x1="{left}" x2="{w-10}" y1="{sy(0):.1f}" y2="{sy(0):.1f}" class="ax"/>']
        for v in (lo + 0.1, 0, hi - 0.1):
            out.append(f'<text x="{left-6}" y="{sy(v)+3:.1f}" class="tk" text-anchor="end">{v:.1f}</text>')
        for i, r in enumerate(rows):
            cx = left + (i + 0.5) * bw
            q1, q3, md = float(r["q1_log2FC"]), float(r["q3_log2FC"]), float(r["median_log2FC"])
            out.append(f'<line x1="{cx:.1f}" x2="{cx:.1f}" y1="{sy(float(r["p10_log2FC"])):.1f}" '
                       f'y2="{sy(float(r["p90_log2FC"])):.1f}" class="lnk"/>')
            out.append(f'<rect x="{cx-bw*0.22:.1f}" y="{sy(q3):.1f}" width="{bw*0.44:.1f}" height="{sy(q1)-sy(q3):.1f}" rx="2" '
                       f'class="bar"><title>{html.escape(r["presence"])}: {int(r["sites"]):,} sites, median {md:.2f}, '
                       f'{r["pct_log2FC_positive"]}% &gt; 0, {r["pct_up_fdr"]}% up at FDR &lt; 0.05</title></rect>')
            out.append(f'<line x1="{cx-bw*0.22:.1f}" x2="{cx+bw*0.22:.1f}" y1="{sy(md):.1f}" y2="{sy(md):.1f}" class="med"/>')
            short = r["presence"].replace("called on ", "").replace(" (third day only)", "")
            for k, part in enumerate(short.split(", ")):
                out.append(f'<text x="{cx:.1f}" y="{base+14+k*11}" class="tk" text-anchor="middle">{html.escape(part)}</text>')
        out.append(f'<text x="10" y="{(top+base)/2}" class="al" transform="rotate(-90 10 {(top+base)/2})">log2FC</text>')
        out.append("</svg>")
        return "".join(out)

    cd = tsv(os.path.join(PO, "class_days.tsv"))
    ab = [r for r in cd if r["called"].upper() == "FALSE"]
    S.append(("Overlap of the per-day peak calls (step 12)",
        "<p>Each site is marked with the days whose consensus peaks contain it (macs3_ge2: MACS3 q 0.01 peaks in at least "
        "2 of the day's 3 libraries, step 05); the 7 combinations of days are the overlap classes. Peak presence is a "
        "threshold call: it depends on each library's depth and signal-to-noise as well as on binding, and carries no "
        "test, so it is checked against the DiffBind counts (background normalization). Gained / lost = FDR &lt; 0.05 up / "
        "down in D18 vs D15 or D22 vs D15.</p>"
        f'<div class="figs">{upset_svg()}</div>' +
        table(["Called on", "Sites", "% of sites", "Median log2FC D18 vs D15", "D22 vs D15", "D22 vs D18", "% gained", "% lost"],
              [[lab[r["class"]], num(r["sites"]), num(r["pct_sites"], 1), num(r["log2FC_D18_vs_D15"], 2),
                num(r["log2FC_D22_vs_D15"], 2), num(r["log2FC_D22_vs_D18"], 2), num(r["pct_gained"], 1), num(r["pct_lost"], 1)]
               for r in cls], num_cols=(1, 2, 3, 4, 5, 6, 7)) +
        "<h3>Library call rates</h3><p>Each library's own MACS3 peaks (step 04) against the sites; lighter bars = batch 3 "
        "libraries. A site is called on a day when 2 of the day's 3 libraries call it, so the day calls follow the "
        "libraries that call most sites.</p>"
        f'<div class="figs">{libs_svg()}</div>' +
        table(["Library", "Day", "Batch", "Reads", "DiffBind size factor", "MACS3 peaks", "Sites called", "% of sites",
               "Reads in sites", "% reads in sites"],
              [[html.escape(r["sample"]), r["day"], r["batch"], num(r["lib_size"]), num(r["norm_factor"], 3), num(r["peaks"]),
                num(r["sites_called"]), num(r["pct_sites_called"], 1), num(r["reads_in_sites"]), num(r["pct_reads_in_sites"], 2)]
               for r in libs], num_cols=(3, 4, 5, 6, 7, 8, 9)) +
        "<h3>Peak presence against the counts</h3><p>Per contrast, the sites called on the later day only, on the earlier "
        "day only, on both and on neither: DiffBind log2FC (box = quartiles, line = median, whiskers = 10th-90th "
        "percentile). If presence differences reflect binding changes, the sites that appear should mostly go up.</p>"
        f'<div class="figs">{"".join(presence_svg(c) for c in ["D18_vs_D15", "D22_vs_D15", "D22_vs_D18"])}</div>' +
        table(["Contrast", "Peak presence", "Sites", "Median log2FC", "% log2FC > 0", "% up, FDR < 0.05", "% down, FDR < 0.05"],
              [[r["contrast"].replace("_vs_", " vs "), html.escape(r["presence"]), num(r["sites"]), num(r["median_log2FC"], 2),
                num(r["pct_log2FC_positive"], 1), num(r["pct_up_fdr"], 1), num(r["pct_down_fdr"], 1)] for r in pc],
              num_cols=(2, 3, 4, 5, 6)) +
        "<p>Library support on the days a site is not called: share of sites called by none or by 1 of the day's 3 "
        "libraries (2 are needed for the consensus), and share called by the batch 3 library.</p>" +
        table(["Called on", "Day not called", "Sites", "% no library calls it", "% 1 of 3 libraries", "% batch 3 library calls it"],
              [[lab[r["class"]], r["day"], num(r["sites"]), num(r["pct_libs0"], 1), num(r["pct_libs1"], 1),
                num(r["pct_batch3_lib"], 1)] for r in ab], num_cols=(2, 3, 4, 5))))

H = os.path.join(OUT, "homer")
if os.path.isdir(H):
    q = {r["item"]: r["n"] for r in tsv(os.path.join(H, "inputs", "qc.tsv"))}
    parts = ["<p>HOMER v5.1 <code>findMotifsGenome.pl</code> on each day's macs3_ge2 consensus peaks: 200 bp windows "
             "centred on the strongest summit, known motifs and de novo motifs (8, 10, 12 bp), repeats masked. "
             "<b>vs genome</b>: random genomic windows matched for GC content (HOMER default). "
             "<b>vs other days</b>: the other two days' windows, without those overlapping a consensus region of the "
             "target day; ISL1's own motif is expected to cancel out here, so this comparison shows which partner "
             "motifs differ between days.</p>",
             table(["Day", "Target windows", "Other-days windows", "Overlap the target day (dropped)",
                    "Background (other days)"],
                   [[d, num(q[f"{d} target windows"]), num(q[f"{d} other-days windows"]),
                     num(q[f"{d} other-days, overlap a {d} region (dropped)"]), num(q[f"{d} other-days background"])]
                    for d in ["D15", "D18", "D22"]], num_cols=(1, 2, 3, 4))]
    for bg, label in [("vs_genome", "vs genome"), ("vs_other_days", "vs other days")]:
        for d in ["D15", "D18", "D22"]:
            k = os.path.join(H, f"{d}_{bg}", "knownResults.txt")
            if not os.path.exists(k):
                parts.append(f"<h3>{d} {label}</h3><p>Not run.</p>")
                continue
            with open(k, newline="") as f:
                rows = list(csv.reader(f, delimiter="\t"))[1:]
            sig = [r for r in rows if float(r[4]) < TH]
            top = rows[:15]
            parts.append(f'<h3>{d} {label}</h3><p>{len(sig)} of {len(rows)} known motifs with q &lt; {TH}; top 15 by '
                         f'p-value. HOMER pages: <code>homer/{d}_{bg}/knownResults.html</code>, '
                         f'<code>homerResults.html</code> (de novo).</p>' + table(
                ["Motif", "log p", "q", "% targets", "% background"],
                [[html.escape(r[0].split("/")[0]), num(r[3], 1), f"{float(r[4]):.1e}", r[6], r[8]] for r in top],
                num_cols=(1, 2, 3, 4)))
    S.append(("HOMER motifs, consensus peaks per day", "".join(parts)))

AN = os.path.join(OUT, "homer", "annotation", "summary.tsv")
if os.path.exists(AN):
    a = tsv(AN)
    sets = ["D15", "D18", "D22", "genome_random"]
    cats = [r["category"] for r in a if r["set"] == "D15" and r["category"] != "total"]
    get = {(r["set"], r["category"]): r for r in a}
    S.append(("Genomic annotation of the consensus peaks", "<p>HOMER <code>annotatePeaks.pl</code> (GENCODE v32) on "
              "the 200 bp summit windows of all macs3_ge2 peaks of each day, and on 50,000 random genomic windows. "
              "promoter-TSS = -1 kb to +100 bp of a TSS; the 10x GTF has no separate UTR features (UTRs count as "
              "exon). Percent of windows.</p>" + table(
                  ["Category", "D15", "D18", "D22", "Random genome"],
                  [[html.escape(c)] + [num(get[(st, c)]["pct"], 1) if (st, c) in get else "–" for st in sets] for c in cats]
                  + [["windows"] + [num(get[(st, "total")]["n"]) for st in sets]], num_cols=(1, 2, 3, 4))))

MH = os.path.join(OUT, "homer", "matched")
if os.path.isdir(MH):
    q = {r["item"]: r["n"] for r in tsv(os.path.join(MH, "inputs", "qc.tsv"))}
    classes = ["promoter", "intron", "Intergenic", "exon", "TTS", "other"]
    parts = ["<p>Backgrounds matched to the peaks' genomic composition. Pools of 200 bp windows that overlap no "
             "consensus peak, the blacklist, or more than 20 N: 150,000 random genomic windows (annotated by HOMER) and "
             "60,000 windows centred within 1 kb of a random GENCODE TSS. <b>A</b>: peaks split into promoter (within "
             "1 kb of a TSS) and distal; promoter peaks vs the promoter pool, distal peaks vs genome windows more than "
             "1 kb from a TSS. <b>B</b>: all peaks vs 3 x as many random windows of each genomic class as the peaks "
             "have. HOMER normalizes GC on top in every run; settings as above.</p>",
             table(["Day", "Class", "Peaks", "Matched background (B)"],
                   [[d, c, num(q.get(f"{d} targets, {c}")), num(q.get(f"{d} matched background, {c}"))]
                    for d in ["D15", "D18", "D22"] for c in classes], num_cols=(2, 3))]
    for run, label in [("promoter_vs_promoters", "A: promoter peaks vs promoters"),
                       ("distal_vs_genome", "A: distal peaks vs genome (> 1 kb from a TSS)"),
                       ("all_vs_matched", "B: all peaks vs composition-matched windows")]:
        for d in ["D15", "D18", "D22"]:
            k = os.path.join(MH, f"{d}_{run}", "knownResults.txt")
            if not os.path.exists(k):
                parts.append(f"<h3>{d}, {label}</h3><p>Not run.</p>")
                continue
            with open(k, newline="") as f:
                rows = list(csv.reader(f, delimiter="\t"))[1:]
            sig = [r for r in rows if float(r[4]) < TH]
            parts.append(f'<h3>{d}, {label}</h3><p>{len(sig)} of {len(rows)} known motifs with q &lt; {TH}; top 15 by '
                         f'p-value. HOMER pages: <code>homer/matched/{d}_{run}/</code>.</p>' + table(
                ["Motif", "log p", "q", "% targets", "% background"],
                [[html.escape(r[0].split("/")[0]), num(r[3], 1), f"{float(r[4]):.1e}", r[6], r[8]] for r in rows[:15]],
                num_cols=(1, 2, 3, 4)))
    S.append(("HOMER motifs, composition-matched backgrounds", "".join(parts)))

S.append(("Files", "<p>Under <code>results/</code> (scratch; keepers copied to the folder's <code>results/</code> on "
          "/pl/active): <code>consensus/&lt;day&gt;_Ctr_ISL1.{macs3_ge2,macs3_3of3}.bed</code>, "
          "<code>tracks/</code> (bigWigs per library and per day, beds, <code>igv_session.xml</code>), and in "
          "<code>diff/</code>: <code>sites_used.bed</code>, "
          "<code>tables/{diffbind,deseq2}_&lt;norm&gt;_&lt;contrast&gt;.tsv</code> (all sites), "
          "<code>beds/diffbind_&lt;norm&gt;_&lt;contrast&gt;_{up,down}.bed</code>, <code>summary.tsv</code>, "
          "<code>agreement.tsv</code>, <code>norm_factors.tsv</code>, <code>counts_sites.tsv</code>, <code>counts_diffbind.tsv</code>; "
          "<code>peak_overlap/</code> (step 12: <code>sites_days.tsv</code>, <code>classes.tsv</code>, "
          "<code>class_days.tsv</code>, <code>libraries.tsv</code>, <code>presence_vs_counts.tsv</code>); "
          "<code>homer/inputs/</code> (target and background beds) and <code>homer/&lt;day&gt;_{vs_genome,vs_other_days}/</code>, <code>homer/annotation/</code>, "
          "<code>homer/matched/</code>.</p>"))

css = """
:root{--bg:#fff;--fg:#1d1d1f;--mut:#5f6368;--line:#dcdcdc;--head:#f3f4f6;--acc:#1f5fa8;--ns:#b8bcc2;--up:#e8710a;--dn:#1a73e8}
@media (prefers-color-scheme:dark){:root{--bg:#16181b;--fg:#e8e8e8;--mut:#a0a4aa;--line:#33363b;--head:#22252a;--acc:#7fb0ea;--ns:#555a61;--up:#f59e4c;--dn:#6ea8fe}}
body{background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;margin:0 auto;max-width:1200px;padding:24px 16px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 8px;color:var(--acc)}h3{font-size:14px;margin:16px 0 4px}
p{margin:6px 0 10px}.sub{color:var(--mut);margin-bottom:16px}.tw{overflow-x:auto}
table{border-collapse:collapse;margin:8px 0 4px;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:5px 10px;text-align:left;white-space:nowrap}
th{background:var(--head);font-weight:600}td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
code{font-size:12px}.figs{display:flex;flex-wrap:wrap;gap:12px}.plot{width:360px;max-width:100%;height:auto}.plot.wide{width:560px}
.plot .pt{font-size:12px;text-anchor:middle;fill:var(--fg)}.plot .al{font-size:10px;text-anchor:middle;fill:var(--mut)}
.plot .tk{font-size:9px;fill:var(--mut)}.plot .ax{stroke:var(--line)}.ns{fill:var(--ns)}.up{fill:var(--up)}.dn{fill:var(--dn)}
.bar{fill:var(--acc)}.med{stroke:var(--bg);stroke-width:2}.lnk{stroke:var(--fg);stroke-width:2}.don{fill:var(--fg)}.doff{fill:var(--line)}
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
