#!/usr/bin/env python3
"""HTML report of the consensus regions (10) and DESeq2 results (11) for every strategy x peak set.

usage: differential_report.py <consensus_dir> <out.html>
"""
import csv
import html
import sys
from collections import Counter
from pathlib import Path

SETS = [("a_conservative", "a · Ctr only · conservative"), ("a_optimal", "a · Ctr only · optimal"),
        ("b_conservative", "b · Ctr + KO+rescue · conservative"), ("b_optimal", "b · Ctr + KO+rescue · optimal")]
DAYS = ["D12", "D15", "D22"]   # PCA colour (3 categorical slots, validated all-pairs)
TOP = 20


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def esc(x):
    return html.escape(str(x))


def label(g):
    day, geno = g.split("_")
    return f"{day} {'Ctr' if geno == 'Ctr' else 'KO+rescue'}"


def fnum(x, fmt):
    try:
        return format(float(x), fmt)
    except ValueError:
        return "–"


def pca_svg(pca, var, samples):
    info = {r["sample"]: r for r in samples}
    xs = [float(r["PC1"]) for r in pca]; ys = [float(r["PC2"]) for r in pca]
    w, h, m = 620, 400, 56
    def sx(v): return m + (v - min(xs)) / ((max(xs) - min(xs)) or 1) * (w - 2 * m)
    def sy(v): return h - m + 10 - (v - min(ys)) / ((max(ys) - min(ys)) or 1) * (h - 2 * m)
    p = [f'<svg viewBox="0 0 {w} {h}" style="max-width:{w}px" role="img" aria-label="PCA">',
         f'<line x1="{m - 20}" y1="{h - 30}" x2="{w - 10}" y2="{h - 30}" class="axis"/>',
         f'<line x1="{m - 20}" y1="10" x2="{m - 20}" y2="{h - 30}" class="axis"/>',
         f'<text x="{(w + m) / 2}" y="{h - 8}" class="al">PC1 ({float(var[0]["percent"]):.0f}%)</text>',
         f'<text transform="translate(14,{h / 2}) rotate(-90)" class="al">PC2 ({float(var[1]["percent"]):.0f}%)</text>']
    for r in pca:
        s = info[r["sample"]]
        x, y = sx(float(r["PC1"])), sy(float(r["PC2"]))
        k = DAYS.index(s["day"]) + 1
        tip = f"{r['sample']} · batch {s['batch']} · PC1 {float(r['PC1']):.1f}, PC2 {float(r['PC2']):.1f}"
        mark = (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" class="pt s{k}"/>' if s["genotype"] == "Ctr" else
                f'<path d="M{x:.1f},{y - 8:.1f} L{x + 7.5:.1f},{y + 6:.1f} L{x - 7.5:.1f},{y + 6:.1f} Z" class="pt s{k}"/>')
        p.append(f'<g data-tip="{esc(tip)}"><circle cx="{x:.1f}" cy="{y:.1f}" r="14" class="hit"/>{mark}'
                 f'<text x="{x + 11:.1f}" y="{y + 4:.1f}" class="pl">{esc(s["day"])} {"Ctr" if s["genotype"] == "Ctr" else "KO"} r{esc(s["replicate"])}·b{esc(s["batch"])}</text></g>')
    p.append("</svg>")
    legend = ("".join(f'<span class="key"><i class="sw s{i + 1}"></i>{d}</span>' for i, d in enumerate(DAYS))
              + '<span class="key"><svg width="14" height="14"><circle cx="7" cy="7" r="5" class="ko-key"/></svg>Ctr</span>'
              + '<span class="key"><svg width="14" height="14"><path d="M7,1 L13,12 L1,12 Z" class="ko-key"/></svg>KO+rescue</span>')
    return f'<div class="legend">{legend}</div>{"".join(p)}'


def render(key, cdir):
    d = Path(cdir) / key
    regions = [l.rstrip("\n").split("\t") for l in open(d / "regions.bed")]
    samples = read_tsv(d / "deseq2" / "samples.tsv")
    out = []
    by_n = Counter(int(r[5]) for r in regions)
    by_groups = Counter(r[4] for r in regions)
    out.append(f"<h2>Consensus regions</h2><p>{len(regions):,} regions: summits ± 250 bp of every group in the strategy, "
               "merged where they overlap.</p><div class='grid2'><div class='scroll'><table><thead><tr><th>Groups supporting</th>"
               "<th>Regions</th></tr></thead><tbody>"
               + "".join(f"<tr><td>{n}</td><td>{by_n[n]:,}</td></tr>" for n in sorted(by_n))
               + "</tbody></table></div><div class='scroll'><table><thead><tr><th class='l'>Exact combination</th><th>Regions</th></tr></thead><tbody>"
               + "".join(f"<tr><td class='l'>{esc(' + '.join(label(g) for g in k.split(',')))}</td><td>{v:,}</td></tr>"
                         for k, v in by_groups.most_common(10))
               + "</tbody></table></div></div>")
    out.append("<h2>Libraries</h2><div class='scroll'><table><thead><tr><th class='l'>Sample</th><th class='l'>Group</th><th>Batch</th>"
               "<th>Fragments in regions</th><th>Size factor</th></tr></thead><tbody>"
               + "".join(f"<tr><td class='l'>{esc(s['sample'])}</td><td class='l'>{esc(label(s['group']))}</td><td>{esc(s['batch'])}</td>"
                         f"<td>{int(float(s['fragments_counted'])):,}</td><td>{float(s['size_factor']):.2f}</td></tr>" for s in samples)
               + "</tbody></table></div><p class='small mute'>Size factors: DESeq2 median-of-ratios on these regions. "
               "Fragments: read pairs counted by featureCounts, duplicates kept (as for peak calling).</p>")
    pca = read_tsv(d / "deseq2" / "pca.tsv"); var = read_tsv(d / "deseq2" / "pca_variance.tsv")
    out.append("<h2>PCA</h2><p>Variance-stabilized counts, 500 most variable regions, no batch correction "
               "(the DESeq2 model includes batch; this plot shows the raw structure). Labels: day, genotype, replicate, batch.</p>"
               + pca_svg(pca, var, samples))
    summ = read_tsv(d / "deseq2" / "de_summary.tsv")
    out.append("<h2>Differential binding</h2><p>DESeq2, <span class='mono'>~ batch + group</span>, Wald test. Significant: "
               "padj &lt; 0.05 and |log2FC| ≥ 1. Up = more ISL1 signal in the first group.</p><div class='scroll'><table><thead><tr>"
               "<th class='l'>Contrast</th><th>Regions tested</th><th>Up</th><th>Down</th></tr></thead><tbody>"
               + "".join(f"<tr><td class='l'>{esc(label(r['a']))} vs {esc(label(r['b']))}</td><td>{int(r['tested']):,}</td>"
                         f"<td>{int(r['up']):,}</td><td>{int(r['down']):,}</td></tr>" for r in summ)
               + "</tbody></table></div>")
    for r in summ:
        rows = read_tsv(d / "deseq2" / f"de_{r['contrast']}.tsv")[:TOP]
        out.append(f"<details><summary>{esc(label(r['a']))} vs {esc(label(r['b']))}: top {TOP} by padj</summary><div class='scroll'><table><thead><tr>"
                   "<th class='l'>Region</th><th class='l'>Location</th><th class='l'>Nearest gene</th><th>Distance to TSS</th>"
                   "<th>Mean count</th><th>log2FC</th><th>padj</th><th class='l'>Called in</th></tr></thead><tbody>"
                   + "".join(f"<tr><td class='l mono'>{esc(x['region'])}</td><td class='l mono'>{esc(x['chr'])}:{int(x['start']):,}-{int(x['end']):,}</td>"
                             f"<td class='l'>{esc(x['SYMBOL'] or '–')}</td><td>{fnum(x['distanceToTSS'], ',.0f')}</td>"
                             f"<td>{fnum(x['baseMean'], ',.0f')}</td><td>{fnum(x['log2FC'], '+.2f')}</td><td>{fnum(x['padj'], '.1e')}</td>"
                             f"<td class='l small'>{esc(', '.join(label(g) for g in x['groups'].split(',')))}</td></tr>" for x in rows)
                   + f"</tbody></table></div><p class='small mute'>Full table: <span class='mono'>results_group/consensus/{key}/deseq2/de_{esc(r['contrast'])}.tsv</span> on Alpine.</p></details>")
    return "".join(out)


CSS = """
:root{color-scheme:light;--surface:#fcfcfb;--raise:#f4f3f0;--ink:#0b0b0b;--ink2:#52514e;--mute:#7a7975;--line:#e4e3df;
--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--surface:#1a1a19;--raise:#232322;--ink:#f0efec;--ink2:#c3c2b7;--mute:#97968f;--line:#353533;
--s1:#3987e5;--s2:#d95926;--s3:#199e70}}
:root[data-theme="dark"]{color-scheme:dark;--surface:#1a1a19;--raise:#232322;--ink:#f0efec;--ink2:#c3c2b7;--mute:#97968f;--line:#353533;
--s1:#3987e5;--s2:#d95926;--s3:#199e70}
*{box-sizing:border-box}
body{margin:0;background:var(--surface);color:var(--ink);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,sans-serif}
main{max-width:1180px;margin:0 auto;padding:28px 16px 60px}
h1{font-size:24px;margin:0 0 6px}h2{font-size:18px;margin:30px 0 8px;padding-top:12px;border-top:1px solid var(--line)}
p{max-width:860px;color:var(--ink2)}.mute{color:var(--mute)}.small{font-size:12px}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
.scroll{overflow-x:auto}.grid2{display:flex;flex-wrap:wrap;gap:24px}
table{border-collapse:collapse;font-variant-numeric:tabular-nums;margin:6px 0 10px}
th,td{padding:5px 9px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th{color:var(--ink2);font-weight:600;font-size:12px;vertical-align:bottom}td.l,th.l{text-align:left}
details{margin:6px 0;padding:6px 10px;background:var(--raise);border-radius:8px}summary{cursor:pointer;font-weight:600;font-size:13px}
.tabs{display:flex;flex-wrap:wrap;gap:6px;margin:14px 0 4px;position:sticky;top:0;background:var(--surface);padding:8px 0;z-index:5}
.tabs button{font:inherit;padding:6px 12px;border-radius:8px;border:1px solid var(--line);background:var(--raise);color:var(--ink);cursor:pointer}
.tabs button[aria-pressed="true"]{background:var(--ink);color:var(--surface);border-color:var(--ink)}
section[hidden]{display:none}
svg{display:block;width:100%;height:auto}
.axis{stroke:var(--line);stroke-width:1}.al{text-anchor:middle;font-size:12px;fill:var(--ink2)}
.pt{stroke:var(--surface);stroke-width:2}.pt.s1{fill:var(--s1)}.pt.s2{fill:var(--s2)}.pt.s3{fill:var(--s3)}
.pl{font-size:11px;fill:var(--ink2)}.hit{fill:transparent}
.legend{display:flex;flex-wrap:wrap;gap:14px;font-size:12px;color:var(--ink2);margin:6px 0}
.key{display:inline-flex;align-items:center;gap:6px}.sw{width:12px;height:12px;border-radius:3px;display:inline-block}
.sw.s1{background:var(--s1)}.sw.s2{background:var(--s2)}.sw.s3{background:var(--s3)}.ko-key{fill:var(--ink2)}
#tip{position:fixed;pointer-events:none;background:var(--ink);color:var(--surface);font-size:12px;padding:5px 8px;border-radius:6px;display:none;max-width:340px;z-index:9}
"""

JS = """
const tip=document.getElementById('tip');
document.addEventListener('mousemove',e=>{const t=e.target.closest('[data-tip]');
if(!t){tip.style.display='none';return}tip.textContent=t.dataset.tip;tip.style.display='block';
tip.style.left=Math.min(e.clientX+14,innerWidth-tip.offsetWidth-8)+'px';tip.style.top=(e.clientY+14)+'px'});
const show=k=>{document.querySelectorAll('section[data-key]').forEach(x=>x.hidden=x.dataset.key!==k);
document.querySelectorAll('.tabs button').forEach(b=>b.setAttribute('aria-pressed',b.dataset.key===k));
try{localStorage.setItem('isl1-diff',k)}catch(e){}};
document.querySelectorAll('.tabs button').forEach(b=>b.addEventListener('click',()=>show(b.dataset.key)));
let k0='a_conservative';try{k0=localStorage.getItem('isl1-diff')||k0}catch(e){}show(k0);
"""


def main(cdir, out_html):
    body = "".join(f'<section data-key="{k}">{render(k, cdir)}</section>' for k, _ in SETS)
    tabs = "".join(f'<button data-key="{k}">{esc(n)}</button>' for k, n in SETS)
    Path(out_html).parent.mkdir(parents=True, exist_ok=True)
    Path(out_html).write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>ISL1 differential binding</title>
<style>{CSS}</style></head><body><main>
<h1>ISL1 CUT&amp;RUN: consensus regions and differential binding</h1>
<p>Strategy <b>a</b>: Ctr groups only (D12, D15, D22); KO+rescue left out. Strategy <b>b</b>: Ctr and KO+rescue groups.
Each strategy is built from the conservative and from the optimal final peak sets. D15 Ctr has one library, so D15 Ctr
contrasts rest on a single sample for that group.</p>
<div class="tabs" role="group" aria-label="Consensus set">{tabs}</div>
{body}
<p class="small mute">Regions and counts: scripts/10_consensus_counts.sb (bedtools merge, featureCounts -p --countReadPairs -B -C).
DESeq2: scripts/deseq2_consensus.R via scripts/11_deseq2.sb. Report: scripts/differential_report.py.</p>
</main><div id="tip"></div><script>{JS}</script></body></html>
""")
    print(f"report: {out_html}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
