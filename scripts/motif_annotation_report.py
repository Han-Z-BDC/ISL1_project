#!/usr/bin/env python3
"""HTML report of HOMER motifs (07) and ChIPseeker annotation / target genes (08) for both final
peak sets (06) of every group.

usage: motif_annotation_report.py <final_peaks_summary.tsv> <homer_dir> <annotation_dir> <out.html>
"""
import csv
import html
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

SETS = ["conservative", "optimal"]
GROUP_ORDER = ["D12_Ctr", "D15_Ctr", "D22_Ctr", "D12_KO", "D15_KO"]
CATS = ["Promoter (<=1 kb)", "Promoter (1-3 kb)", "UTR / exon", "Intron", "Downstream (<=300 bp)",
        "Distal intergenic", "Other"]
BINS = ["< -100 kb", "-100 to -10 kb", "-10 to -1 kb", "-1 to 0 kb", "0 to 1 kb", "1 to 10 kb",
        "10 to 100 kb", "> 100 kb"]
TOP_KNOWN, TOP_DENOVO, TOP_GENES = 10, 5, 25
PANEL_NEAR = 50000   # panel cells beyond this distance show as "–"


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def esc(x):
    return html.escape(str(x))


def label(g):
    day, geno = g.split("_")
    return f"{day} {'Ctr' if geno == 'Ctr' else 'KO+rescue'}"


def short_motif(name):
    return name.split("/")[0]


def known(homer_dir, s, g):
    p = Path(homer_dir) / s / g / "knownResults.txt"
    if not p.exists():
        return []
    rows = []
    with open(p) as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        for x in r:
            rows.append({"name": x[0], "consensus": x[1], "neglog10p": -float(x[3]) / math.log(10),
                         "q": x[4], "pct_t": x[6], "pct_b": x[8]})
    return rows


def denovo(homer_dir, s, g):
    d = Path(homer_dir) / s / g / "homerResults"
    out = []
    for i in range(1, TOP_DENOVO + 1):
        m = d / f"motif{i}.motif"
        if not m.exists():
            break
        hdr = m.read_text().splitlines()[0].lstrip(">").split("\t")
        best = re.sub(r"^.*BestGuess:", "", hdr[1]) if len(hdr) > 1 else ""
        stats = dict(kv.split(":", 1) for kv in hdr[5].split(",") if ":" in kv) if len(hdr) > 5 else {}
        svg = d / f"motif{i}.logo.svg"
        logo = re.sub(r"<\?xml[^>]*\?>|<!DOCTYPE[^>]*>", "", svg.read_text()) if svg.exists() else ""
        out.append({"consensus": hdr[0], "best": best, "p": stats.get("P", ""), "t": stats.get("T", ""),
                    "b": stats.get("B", ""), "logo": logo})
    return out


def stacked_bars(rows_by_group, groups, title):
    """Horizontal 100% stacked bars, one per group, categories in fixed slot order."""
    left, bw, bh, gap = 120, 560, 22, 10
    h = len(groups) * (bh + gap) + 6
    p = [f'<svg viewBox="0 0 {left + bw + 10} {h}" style="max-width:{left + bw + 10}px" role="img" aria-label="{esc(title)}">']
    for i, g in enumerate(groups):
        y = i * (bh + gap) + 2
        counts = rows_by_group.get(g, {})
        tot = sum(counts.values()) or 1
        p.append(f'<text x="{left - 10}" y="{y + bh / 2 + 4}" class="rl">{esc(label(g))}</text>')
        x = left
        for k, c in enumerate(CATS):
            n = counts.get(c, 0)
            if not n:
                continue
            w = bw * n / tot
            tip = f"{label(g)} · {c}: {n:,} peaks ({100 * n / tot:.1f}%)"
            p.append(f'<rect x="{x:.1f}" y="{y}" width="{max(0.5, w - 2):.1f}" height="{bh}" rx="3" class="c{k}" data-tip="{esc(tip)}"/>')
            x += w
    p.append("</svg>")
    legend = "".join(f'<span class="key"><i class="sw c{k}"></i>{esc(c)}</span>' for k, c in enumerate(CATS))
    return f'<div class="legend">{legend}</div>{"".join(p)}'


def heat_class(v, vmax):
    return f"h{min(6, int(7 * v / vmax))}" if vmax > 0 else "h0"


def render_set(s, groups, fps, homer_dir, ann_dir, ann_sum, tss_sum, panel):
    out = []
    # Overview
    rows = []
    kn = {g: known(homer_dir, s, g) for g in groups}
    tg = {g: read_tsv(Path(ann_dir) / s / f"{g}.target_genes.tsv") for g in groups}
    for g in groups:
        c = ann_sum[s].get(g, {})
        tot = sum(c.values()) or 1
        prom = c.get(CATS[0], 0) + c.get(CATS[1], 0)
        top = kn[g][0] if kn[g] else None
        top_p = f"{top['neglog10p']:.0f}" if top else "–"
        rows.append(f"<tr><td class='l'>{esc(label(g))}</td><td>{int(fps[(g, s)]['final_peaks']):,}</td>"
                    f"<td class='l'>{esc(fps[(g, s)]['source'].replace('_', ' '))}</td>"
                    f"<td>{100 * prom / tot:.0f}%</td><td>{100 * c.get('Intron', 0) / tot:.0f}%</td>"
                    f"<td>{100 * c.get('Distal intergenic', 0) / tot:.0f}%</td>"
                    f"<td>{len(tg[g]):,}</td><td>{sum(r['promoter_bound'] == 'TRUE' for r in tg[g]):,}</td>"
                    f"<td class='l'>{esc(short_motif(top['name'])) if top else '–'}</td>"
                    f"<td>{top_p}</td></tr>")
    out.append("<h2>Overview</h2><div class='scroll'><table><thead><tr><th class='l'>Group</th><th>Peaks</th>"
               "<th class='l'>Source</th><th>Promoter<br>(±3 kb)</th><th>Intron</th><th>Distal<br>intergenic</th>"
               "<th>Nearest<br>genes</th><th>Promoter-bound<br>genes</th><th class='l'>Top known motif</th>"
               "<th>−log10 p</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"
               "<p class='small mute'>KO+rescue sets are small (under 250 peaks), so their motif and annotation "
               "results rest on few sites.</p>")

    # Genomic annotation
    out.append("<h2>Where the peaks are</h2><p>Annotation of each summit against UCSC knownGene (promoter = TSS ± 3 kb). "
               "Hover a segment for counts.</p>")
    out.append(stacked_bars(ann_sum[s], groups, "Genomic annotation"))
    out.append("<details><summary class='small'>Table view</summary><div class='scroll'><table><thead><tr><th class='l'>Group</th>"
               + "".join(f"<th>{esc(c)}</th>" for c in CATS) + "</tr></thead><tbody>"
               + "".join("<tr><td class='l'>" + esc(label(g)) + "</td>" + "".join(
                   f"<td>{ann_sum[s].get(g, {}).get(c, 0):,}</td>" for c in CATS) + "</tr>" for g in groups)
               + "</tbody></table></div></details>")
    out.append("<h3>Distance from summit to nearest TSS</h3><div class='scroll'><table><thead><tr><th class='l'>Group</th>"
               + "".join(f"<th>{esc(b)}</th>" for b in BINS) + "</tr></thead><tbody>")
    for g in groups:
        c = tss_sum[s].get(g, {}); tot = sum(c.values()) or 1
        out.append("<tr><td class='l'>" + esc(label(g)) + "</td>"
                   + "".join(f"<td>{100 * c.get(b, 0) / tot:.0f}%</td>" for b in BINS) + "</tr>")
    out.append("</tbody></table></div><p class='small mute'>Negative = upstream of the TSS (strand-aware).</p>")

    # Known motifs: comparison heat table, then per-group detail
    union = []
    for g in groups:
        for r in kn[g][:5]:
            if r["name"] not in union:
                union.append(r["name"])
    vals = {(g, r["name"]): r["neglog10p"] for g in groups for r in kn[g]}
    vmax = max([vals.get((g, m), 0) for g in groups for m in union] or [1])
    out.append("<h2>Known motifs (HOMER)</h2><p>Top 5 known motifs of each group, compared across groups as −log10 p "
               "(darker = stronger). 200 bp around summits vs GC-matched genomic background.</p>")
    out.append("<div class='scroll'><table class='heatt'><thead><tr><th class='l'>Motif</th>"
               + "".join(f"<th>{esc(label(g))}</th>" for g in groups) + "</tr></thead><tbody>")
    for m in union:
        cells = []
        for g in groups:
            v = vals.get((g, m))
            cells.append("<td>–</td>" if v is None else
                         f"<td class='{heat_class(v, vmax)}' data-tip='{esc(short_motif(m))} in {esc(label(g))}: −log10 p {v:.0f}'>{v:.0f}</td>")
        out.append(f"<tr><td class='l' title='{esc(m)}'>{esc(short_motif(m))}</td>{''.join(cells)}</tr>")
    out.append("</tbody></table></div>")
    for g in groups:
        out.append(f"<details><summary>{esc(label(g))}: top {TOP_KNOWN} known motifs</summary><div class='scroll'><table><thead><tr>"
                   "<th class='l'>Motif</th><th class='l'>Consensus</th><th>−log10 p</th><th>q</th><th>% targets</th>"
                   "<th>% background</th></tr></thead><tbody>"
                   + "".join(f"<tr><td class='l' title='{esc(r['name'])}'>{esc(short_motif(r['name']))}</td>"
                             f"<td class='l mono'>{esc(r['consensus'])}</td><td>{r['neglog10p']:.0f}</td><td>{esc(r['q'])}</td>"
                             f"<td>{esc(r['pct_t'])}</td><td>{esc(r['pct_b'])}</td></tr>" for r in kn[g][:TOP_KNOWN])
                   + "</tbody></table></div></details>")

    # De novo
    out.append("<h2>De novo motifs (HOMER)</h2><p>Top de novo motifs per group with HOMER's best match to known motifs.</p>")
    for g in groups:
        dn = denovo(homer_dir, s, g)
        items = "".join(
            f"<div class='motif'><div class='logo'>{m['logo'] or '<span class=mono>' + esc(m['consensus']) + '</span>'}</div>"
            f"<div><div class='mono'>{esc(m['consensus'])}</div><div class='small'>{esc(short_motif(m['best']))}</div>"
            f"<div class='small mute'>p {esc(m['p'])} · targets {esc(m['t'])} · background {esc(m['b'])}</div></div></div>"
            for m in dn) or "<p class='small mute'>No de novo motifs.</p>"
        out.append(f"<details><summary>{esc(label(g))}</summary>{items}</details>")

    # Target genes
    out.append("<h2>Target genes</h2><p>Nearest gene of each summit (ChIPseeker). Promoter-bound = a summit within ±3 kb "
               f"of the TSS. Top {TOP_GENES} by peak strength; full lists in "
               f"<span class='mono'>results_group/annotation/{s}/&lt;group&gt;.target_genes.tsv</span> on Alpine.</p>")
    ctr = [g for g in groups if g.endswith("_Ctr")]
    if len(ctr) == 3:
        sets_all = {g: {r["SYMBOL"] for r in tg[g]} for g in ctr}
        sets_pro = {g: {r["SYMBOL"] for r in tg[g] if r["promoter_bound"] == "TRUE"} for g in ctr}
        def combos(S):
            a, b, c = (S[g] for g in ctr)
            return [("all three", len(a & b & c)), (f"{label(ctr[0])} only", len(a - b - c)),
                    (f"{label(ctr[1])} only", len(b - a - c)), (f"{label(ctr[2])} only", len(c - a - b)),
                    (f"{label(ctr[0])} + {label(ctr[1])} only", len((a & b) - c)),
                    (f"{label(ctr[0])} + {label(ctr[2])} only", len((a & c) - b)),
                    (f"{label(ctr[1])} + {label(ctr[2])} only", len((b & c) - a))]
        ca, cp = combos(sets_all), combos(sets_pro)
        out.append("<h3>Shared across Ctr stages</h3><div class='scroll'><table><thead><tr><th class='l'>Genes bound in</th>"
                   "<th>nearest genes</th><th>promoter-bound</th></tr></thead><tbody>"
                   + "".join(f"<tr><td class='l'>{esc(k)}</td><td>{v:,}</td><td>{cp[i][1]:,}</td></tr>" for i, (k, v) in enumerate(ca))
                   + "</tbody></table></div>")
    for g in groups:
        out.append(f"<details><summary>{esc(label(g))}: top {TOP_GENES} of {len(tg[g]):,} genes</summary><div class='scroll'><table><thead><tr>"
                   "<th class='l'>Gene</th><th>Peaks</th><th>Best −log10 p</th><th>Closest summit<br>to TSS (kb)</th>"
                   "<th class='l'>Promoter</th></tr></thead><tbody>"
                   + "".join(f"<tr><td class='l'>{esc(r['SYMBOL'])}</td><td>{esc(r['n_peaks'])}</td>"
                             f"<td>{float(r['best_neglog10p']):.1f}</td><td>{float(r['min_abs_dist_to_tss']) / 1000:.1f}</td>"
                             f"<td class='l'>{'yes' if r['promoter_bound'] == 'TRUE' else ''}</td></tr>" for r in tg[g][:TOP_GENES])
                   + "</tbody></table></div></details>")

    # Islet panel
    genes = []
    for r in panel:
        if r["set"] == s and r["gene"] not in genes:
            genes.append(r["gene"])
    pd = {(r["group"], r["gene"]): r["dist_bp"] for r in panel if r["set"] == s}
    out.append("<h2>Islet gene panel</h2><p>Distance from each gene's TSS to the nearest ISL1 summit, in kb. "
               "<b>P</b> = within ±3 kb (promoter); – = none within 50 kb.</p><div class='scroll'><table><thead><tr><th class='l'>Gene</th>"
               + "".join(f"<th>{esc(label(g))}</th>" for g in groups) + "</tr></thead><tbody>")
    for gene in genes:
        cells = []
        for g in groups:
            d = pd.get((g, gene), "NA")
            if d in ("NA", ""):
                cells.append("<td>–</td>")
                continue
            d = float(d)
            cells.append("<td>–</td>" if d > PANEL_NEAR else
                         f"<td class='{'pb' if d <= 3000 else ''}'>{d / 1000:.1f}{' P' if d <= 3000 else ''}</td>")
        out.append(f"<tr><td class='l'>{esc(gene)}</td>{''.join(cells)}</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


CSS = """
:root{color-scheme:light;--surface:#fcfcfb;--raise:#f4f3f0;--ink:#0b0b0b;--ink2:#52514e;--mute:#7a7975;--line:#e4e3df;
--c0:#2a78d6;--c1:#eb6834;--c2:#1baf7a;--c3:#eda100;--c4:#e87ba4;--c5:#008300;--c6:#a8a7a1;--pb:#e3eefc;
--h0:#cde2fb;--h1:#9ec5f4;--h2:#6da7ec;--h3:#3987e5;--h4:#256abf;--h5:#184f95;--h6:#0d366b;
--t0:#0b0b0b;--t1:#0b0b0b;--t2:#0b0b0b;--t3:#fff;--t4:#fff;--t5:#fff;--t6:#fff}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--surface:#1a1a19;--raise:#232322;--ink:#f0efec;--ink2:#c3c2b7;--mute:#97968f;--line:#353533;
--c0:#3987e5;--c1:#d95926;--c2:#199e70;--c3:#c98500;--c4:#d55181;--c5:#008300;--c6:#6b6a65;--pb:#1d3352;
--h0:#0d366b;--h1:#104281;--h2:#184f95;--h3:#256abf;--h4:#3987e5;--h5:#6da7ec;--h6:#b7d3f6;
--t0:#f0efec;--t1:#f0efec;--t2:#f0efec;--t3:#f0efec;--t4:#0b0b0b;--t5:#0b0b0b;--t6:#0b0b0b}}
:root[data-theme="dark"]{color-scheme:dark;--surface:#1a1a19;--raise:#232322;--ink:#f0efec;--ink2:#c3c2b7;--mute:#97968f;--line:#353533;
--c0:#3987e5;--c1:#d95926;--c2:#199e70;--c3:#c98500;--c4:#d55181;--c5:#008300;--c6:#6b6a65;--pb:#1d3352;
--h0:#0d366b;--h1:#104281;--h2:#184f95;--h3:#256abf;--h4:#3987e5;--h5:#6da7ec;--h6:#b7d3f6;
--t0:#f0efec;--t1:#f0efec;--t2:#f0efec;--t3:#f0efec;--t4:#0b0b0b;--t5:#0b0b0b;--t6:#0b0b0b}
*{box-sizing:border-box}
body{margin:0;background:var(--surface);color:var(--ink);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,sans-serif}
main{max-width:1180px;margin:0 auto;padding:28px 16px 60px}
h1{font-size:24px;margin:0 0 6px}h2{font-size:18px;margin:34px 0 8px;padding-top:12px;border-top:1px solid var(--line)}h3{font-size:15px;margin:18px 0 6px}
p{max-width:860px;color:var(--ink2)}.mute{color:var(--mute)}.small{font-size:12px}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
.scroll{overflow-x:auto}
table{border-collapse:collapse;font-variant-numeric:tabular-nums;margin:6px 0 10px}
th,td{padding:5px 9px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th{color:var(--ink2);font-weight:600;font-size:12px;vertical-align:bottom}td.l,th.l{text-align:left}
details{margin:6px 0;padding:6px 10px;background:var(--raise);border-radius:8px}summary{cursor:pointer;font-weight:600;font-size:13px}
.tabs{display:flex;gap:6px;margin:14px 0 4px;position:sticky;top:0;background:var(--surface);padding:8px 0;z-index:5}
.tabs button{font:inherit;padding:6px 14px;border-radius:8px;border:1px solid var(--line);background:var(--raise);color:var(--ink);cursor:pointer}
.tabs button[aria-pressed="true"]{background:var(--ink);color:var(--surface);border-color:var(--ink)}
section[hidden]{display:none}
svg{display:block;width:100%;height:auto}
.rl{text-anchor:end;font-size:12px;fill:var(--ink2)}
.legend{display:flex;flex-wrap:wrap;gap:12px;font-size:12px;color:var(--ink2);margin:6px 0}
.key{display:inline-flex;align-items:center;gap:6px}.sw{width:12px;height:12px;border-radius:3px;display:inline-block}
.c0{fill:var(--c0);background:var(--c0)}.c1{fill:var(--c1);background:var(--c1)}.c2{fill:var(--c2);background:var(--c2)}
.c3{fill:var(--c3);background:var(--c3)}.c4{fill:var(--c4);background:var(--c4)}.c5{fill:var(--c5);background:var(--c5)}.c6{fill:var(--c6);background:var(--c6)}
.heatt td{min-width:70px}
.h0{background:var(--h0);color:var(--t0)}.h1{background:var(--h1);color:var(--t1)}.h2{background:var(--h2);color:var(--t2)}
.h3{background:var(--h3);color:var(--t3)}.h4{background:var(--h4);color:var(--t4)}.h5{background:var(--h5);color:var(--t5)}.h6{background:var(--h6);color:var(--t6)}
td.pb{background:var(--pb);font-weight:600}
.motif{display:flex;gap:14px;align-items:center;padding:8px 0;border-bottom:1px solid var(--line)}
.logo{width:240px;flex:none;background:#fff;border-radius:6px;padding:2px}.logo svg{width:100%;height:auto}
#tip{position:fixed;pointer-events:none;background:var(--ink);color:var(--surface);font-size:12px;padding:5px 8px;border-radius:6px;display:none;max-width:320px;z-index:9}
"""

JS = """
const tip=document.getElementById('tip');
document.addEventListener('mousemove',e=>{const t=e.target.closest('[data-tip]');
if(!t){tip.style.display='none';return}tip.textContent=t.dataset.tip;tip.style.display='block';
tip.style.left=Math.min(e.clientX+14,innerWidth-tip.offsetWidth-8)+'px';tip.style.top=(e.clientY+14)+'px'});
const show=s=>{document.querySelectorAll('section[data-set]').forEach(x=>x.hidden=x.dataset.set!==s);
document.querySelectorAll('.tabs button').forEach(b=>b.setAttribute('aria-pressed',b.dataset.set===s));
try{localStorage.setItem('isl1-set',s)}catch(e){}};
document.querySelectorAll('.tabs button').forEach(b=>b.addEventListener('click',()=>show(b.dataset.set)));
let s0='conservative';try{s0=localStorage.getItem('isl1-set')||s0}catch(e){}show(s0);
"""


def main(fps_tsv, homer_dir, ann_dir, out_html):
    fps = {(r["group"], r["set"]): r for r in read_tsv(fps_tsv)}
    groups = [g for g in GROUP_ORDER if (g, "conservative") in fps]
    ann_sum = {s: defaultdict(dict) for s in SETS}
    for r in read_tsv(Path(ann_dir) / "annotation_summary.tsv"):
        ann_sum[r["set"]][r["group"]][r["category"]] = int(r["n"])
    tss_sum = {s: defaultdict(dict) for s in SETS}
    for r in read_tsv(Path(ann_dir) / "tss_distance_summary.tsv"):
        tss_sum[r["set"]][r["group"]][r["bin"]] = int(r["n"])
    panel = read_tsv(Path(ann_dir) / "islet_panel.tsv")

    body = "".join(f'<section data-set="{s}">{render_set(s, groups, fps, homer_dir, ann_dir, ann_sum, tss_sum, panel)}</section>'
                   for s in SETS)
    Path(out_html).parent.mkdir(parents=True, exist_ok=True)
    Path(out_html).write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>ISL1 motifs and targets</title>
<style>{CSS}</style></head><body><main>
<h1>ISL1 CUT&amp;RUN: motifs, genomic annotation and target genes</h1>
<p>Final ISL1 peak sets per group: <b>conservative</b> (true-replicate IDR) and <b>optimal</b> (ENCODE: the larger of
true-replicate and pooled-pseudoreplicate IDR). D15 Ctr has one replicate (self-pseudoreplicate IDR, unreplicated).
"KO+rescue" = ISL1 KO with ISL1 rescue. Switch sets with the buttons; hover cells and bars for detail.</p>
<div class="tabs" role="group" aria-label="Peak set"><button data-set="conservative">Conservative</button><button data-set="optimal">Optimal</button></div>
{body}
<p class="small mute">HOMER findMotifsGenome.pl (-size 200, summit-centred, hg38) via scripts/07_homer_motifs.sb; ChIPseeker
annotatePeak (TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, TSS ± 3 kb) via scripts/08_annotate_peaks.sb;
report by scripts/motif_annotation_report.py.</p>
</main><div id="tip"></div><script>{JS}</script></body></html>
""")
    print(f"report: {out_html}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
