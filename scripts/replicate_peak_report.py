#!/usr/bin/env python3
"""Per-replicate ISL1 peak statistics and replicate overlap for the final group-IgG calls.

Writes replicate_peak_report.html plus replicate_stats.tsv, replicate_overlap.tsv and
peak_jaccard.tsv to <out_dir>.

usage: replicate_peak_report.py <samples.tsv> <sample_metadata.tsv> <qc_alignment_dir>
           <group_macs3_dir> <group_relaxed_dir> <idr_dir> <frip.tsv> <bedtools> <out_dir>
"""
import csv
import html
import re
import statistics
import subprocess
import sys
import tempfile
from collections import Counter
from itertools import combinations
from pathlib import Path

TOP_N = 1000            # strongest peaks (by -log10 q) per replicate for the top-peak overlap
RELAXED_CAP = 300000    # 03c keeps the top 300k relaxed peaks
WEAK_TOP = 0.70         # note a replicate whose top peaks are recovered below this in a partner
GROUP_ORDER = ["D12_Ctr", "D12_KO", "D15_Ctr", "D15_KO", "D22_Ctr"]


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def unique_peaks(narrowpeak):
    """--call-summits repeats an interval once per summit; keep each interval once, best q."""
    best = {}
    for line in open(narrowpeak):
        r = line.rstrip("\n").split("\t")
        if len(r) < 9:
            continue
        key = (r[0], int(r[1]), int(r[2]))
        best[key] = max(best.get(key, float("-inf")), float(r[8]))
    return best


def write_bed(keys, path, index=False):
    with open(path, "w") as f:
        for i, (c, s, e) in enumerate(sorted(keys)):
            f.write(f"{c}\t{s}\t{e}" + (f"\t{i}" if index else "") + "\n")


def bedtools(tool, *args):
    return subprocess.run([tool, *args], capture_output=True, text=True, check=True).stdout


def n_hit(tool, a, b):
    return sum(1 for l in bedtools(tool, "intersect", "-u", "-a", a, "-b", b).splitlines() if l)


def flagstat(path):
    txt = open(path).read()
    get = lambda pat: int(re.search(pat, txt, re.M).group(1))
    total = get(r"^(\d+) \+ \d+ in total")
    dup = get(r"^(\d+) \+ \d+ duplicates$")
    read1 = get(r"^(\d+) \+ \d+ read1$")
    return {"fragments": read1, "dup_rate": dup / total, "unique_fragments": (total - dup) // 2}


def membership(tool, tmp, tag, sets):
    """Venn counts over merged regions: union of all replicates' peaks, merged; each merged
    region is assigned the set of replicates with a peak overlapping it."""
    names = list(sets)
    for n in names:
        write_bed(sets[n], f"{tmp}/{tag}.{n}.bed")
    write_bed(set().union(*sets.values()), f"{tmp}/{tag}.all.bed")
    merged = bedtools(tool, "merge", "-i", f"{tmp}/{tag}.all.bed")
    with open(f"{tmp}/{tag}.union.bed", "w") as f:
        for i, l in enumerate(merged.splitlines()):
            f.write(f"{l}\t{i}\n")
    n_union = len(merged.splitlines())
    member = [0] * n_union
    for bit, n in enumerate(names):
        for l in bedtools(tool, "intersect", "-u", "-a", f"{tmp}/{tag}.union.bed",
                          "-b", f"{tmp}/{tag}.{n}.bed").splitlines():
            member[int(l.split("\t")[3])] |= 1 << bit
    return names, Counter(member), n_union


def main(samples_tsv, metadata_tsv, qc_dir, macs3_dir, relaxed_dir, idr_dir, frip_tsv, tool, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    meta = {r["sample_name"]: r for r in read_tsv(metadata_tsv)}
    frip = {r["sample"]: r for r in read_tsv(frip_tsv)}
    tf = [r["sample_name"] for r in read_tsv(samples_tsv) if r["role"] == "TF"]
    group_of = {s: f"{meta[s]['day']}_{meta[s]['genotype']}" for s in tf}
    groups = {g: sorted((s for s in tf if group_of[s] == g), key=lambda s: int(meta[s]["replicate"]))
              for g in GROUP_ORDER if any(group_of[s] == g for s in tf)}

    peaks = {s: unique_peaks(Path(macs3_dir) / f"{s}_peaks.filtered.narrowPeak") for s in tf}
    top = {s: set(sorted(peaks[s], key=peaks[s].get, reverse=True)[:TOP_N]) for s in tf}

    stats = {}
    for s in tf:
        fs = flagstat(Path(qc_dir) / f"{s}.flagstat.txt")
        widths = [e - b for _, b, e in peaks[s]]
        relaxed = sum(1 for _ in open(Path(relaxed_dir) / f"{s}_peaks.relaxed.narrowPeak"))
        stats[s] = {
            "sample": s, "group": group_of[s], "replicate": int(meta[s]["replicate"]),
            "batch": meta[s]["batch"], **fs, "peaks": len(peaks[s]),
            "median_width": int(statistics.median(widths)) if widths else 0,
            "frip": int(frip[s]["in_peaks"]) / int(frip[s]["fragments"]),
            "relaxed_peaks": relaxed,
        }

    pair_rows, venns = [], {}
    with tempfile.TemporaryDirectory() as tmp:
        for s in tf:
            write_bed(peaks[s], f"{tmp}/{s}.all.bed")
            write_bed(top[s], f"{tmp}/{s}.top.bed")

        # Pairwise recovery within each group, both directions.
        for g, reps in groups.items():
            for a, b in ((x, y) for x in reps for y in reps if x != y):
                pair_rows.append({
                    "group": g, "a": a, "b": b,
                    "a_in_b": n_hit(tool, f"{tmp}/{a}.all.bed", f"{tmp}/{b}.all.bed") / max(1, len(peaks[a])),
                    "a_top_in_b": n_hit(tool, f"{tmp}/{a}.top.bed", f"{tmp}/{b}.all.bed") / max(1, len(top[a])),
                })
            if len(reps) >= 2:
                venns[g] = {
                    "all": membership(tool, tmp, f"{g}.all", {r: set(peaks[r]) for r in reps}),
                    "top": membership(tool, tmp, f"{g}.top", {r: top[r] for r in reps}),
                }

        # KO peaks at same-day Ctr ISL1 sites (union of that day's Ctr replicate peaks).
        for g, reps in groups.items():
            day, geno = g.split("_")
            ctr = groups.get(f"{day}_Ctr")
            if geno != "KO" or not ctr:
                continue
            write_bed(set().union(*(peaks[c] for c in ctr)), f"{tmp}/{day}.ctr_union.bed")
            for s in reps:
                stats[s]["at_ctr_sites"] = n_hit(tool, f"{tmp}/{s}.all.bed", f"{tmp}/{day}.ctr_union.bed") / max(1, len(peaks[s]))

        # Jaccard (bp) between every pair of replicates.
        jac = {}
        for a, b in combinations(tf, 2):
            hdr, val = bedtools(tool, "jaccard", "-a", f"{tmp}/{a}.all.bed", "-b", f"{tmp}/{b}.all.bed").splitlines()[:2]
            jac[(a, b)] = jac[(b, a)] = float(val.split()[hdr.split().index("jaccard")])

    for s in tf:
        rec = [p["a_top_in_b"] for p in pair_rows if p["a"] == s]
        got = [p["a_top_in_b"] for p in pair_rows if p["b"] == s]
        stats[s]["top_recovered_min"] = min(rec) if rec else None
        stats[s]["partner_top_found_min"] = min(got) if got else None
        notes = []
        if stats[s]["relaxed_peaks"] >= 0.99 * RELAXED_CAP:
            notes.append("relaxed calls at the 300k cap")
        # Top-peak notes for Ctr only: KO+rescue groups have a few hundred weak peaks per
        # replicate, so top-1000 recovery there is dominated by small numbers.
        is_ctr = group_of[s].endswith("_Ctr")
        if is_ctr and got and min(got) < WEAK_TOP:
            notes.append(f"misses partner top peaks (as low as {100 * min(got):.0f}%)")
        if is_ctr and rec and min(rec) < WEAK_TOP:
            notes.append(f"own top peaks weakly recovered (as low as {100 * min(rec):.0f}%)")
        if stats[s].get("at_ctr_sites") is not None and stats[s]["peaks"] > 5000:
            notes.append(f"{stats[s]['peaks']:,} KO peaks, {100 * stats[s]['at_ctr_sites']:.0f}% at Ctr sites")
        stats[s]["notes"] = "; ".join(notes)

    idr_pairs = read_tsv(Path(idr_dir) / "idr_summary.tsv")
    self_idr = {r["sample"]: int(r["N_self"]) for r in read_tsv(Path(idr_dir) / "self_pseudorep_idr.tsv")}
    qc = {r["group"]: r for r in read_tsv(Path(idr_dir) / "encode_idr_qc.tsv")}
    for s in tf:
        stats[s]["self_pseudorep_idr"] = self_idr.get(s)

    cols = ["sample", "group", "replicate", "batch", "fragments", "dup_rate", "unique_fragments", "peaks",
            "median_width", "frip", "relaxed_peaks", "top_recovered_min", "partner_top_found_min",
            "at_ctr_sites", "self_pseudorep_idr", "notes"]
    with open(out / "replicate_stats.tsv", "w") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(stats[s] for s in tf)
    with open(out / "replicate_overlap.tsv", "w") as f:
        w = csv.DictWriter(f, ["group", "a", "b", "a_in_b", "a_top_in_b"], delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)
    order = [s for g in groups.values() for s in g]
    with open(out / "peak_jaccard.tsv", "w") as f:
        f.write("\t".join(["sample"] + order) + "\n")
        for a in order:
            f.write("\t".join([a] + [("1" if a == b else f"{jac[(a, b)]:.4f}") for b in order]) + "\n")

    (out / "replicate_peak_report.html").write_text(
        render(groups, stats, pair_rows, venns, jac, order, idr_pairs, qc, meta))
    print(f"report: {out / 'replicate_peak_report.html'}")


# ---------------------------------------------------------------- rendering

def esc(x):
    return html.escape(str(x))


def short(s, meta):
    return f"rep{meta[s]['replicate']}"


def pct(v):
    return "–" if v is None else f"{100 * v:.0f}%"


VENN2 = {"circles": [(140, 125), (240, 125)], "r": 95, "w": 380, "h": 250,
         "labels": {0b01: (92, 125), 0b10: (288, 125), 0b11: (190, 125)}}
VENN3 = {"circles": [(150, 118), (250, 118), (200, 205)], "r": 90, "w": 400, "h": 315,
         "labels": {0b001: (112, 92), 0b010: (288, 92), 0b100: (200, 262), 0b011: (200, 88),
                    0b101: (150, 182), 0b110: (250, 182), 0b111: (200, 148)}}


def svg_venn(names, counts, n_union, meta, title):
    lay = VENN2 if len(names) == 2 else VENN3
    parts = [f'<svg viewBox="0 0 {lay["w"]} {lay["h"]}" class="venn" role="img" aria-label="{esc(title)}">']
    for i, (cx, cy) in enumerate(lay["circles"]):
        parts.append(f'<circle cx="{cx}" cy="{cy}" r="{lay["r"]}" class="vc s{i + 1}"/>')
    for mask, (x, y) in lay["labels"].items():
        n = counts.get(mask, 0)
        who = " ∩ ".join(short(names[b], meta) for b in range(len(names)) if mask >> b & 1)
        tip = f"{who} only: {n:,} regions ({100 * n / max(1, n_union):.1f}% of {n_union:,})"
        parts.append(f'<g data-tip="{esc(tip)}"><circle cx="{x}" cy="{y}" r="24" class="hit"/>'
                     f'<text x="{x}" y="{y + 5}" class="vn">{n:,}</text></g>')
    parts.append("</svg>")
    legend = "".join(
        f'<span class="key"><i class="sw s{i + 1}"></i>{esc(short(n, meta))} · batch {esc(meta[n]["batch"])}</span>'
        for i, n in enumerate(names))
    return (f'<figure class="panel"><figcaption>{esc(title)}</figcaption>{"".join(parts)}'
            f'<div class="legend">{legend}</div>'
            f'<div class="mute small">{n_union:,} merged regions; areas not to scale</div></figure>')


def svg_upset(names, counts, n_union, meta, title, set_sizes):
    masks = sorted((m for m in counts if m), key=lambda m: -counts[m])
    colw, left, barh, rowh = 34, 150, 150, 24
    w = left + colw * len(masks) + 10
    h = barh + 30 + rowh * len(names) + 10
    vmax = max(counts[m] for m in masks)
    p = [f'<svg viewBox="0 0 {w} {h}" class="upset" style="max-width:{w}px" role="img" aria-label="{esc(title)}">']
    p.append(f'<line x1="{left - 4}" y1="{barh + 20}" x2="{w - 6}" y2="{barh + 20}" class="axis"/>')
    for j, m in enumerate(masks):
        x = left + j * colw + colw / 2
        bh = max(6, barh * counts[m] / vmax)
        who = " ∩ ".join(short(names[b], meta) for b in range(len(names)) if m >> b & 1)
        tip = f"{who} only: {counts[m]:,} regions ({100 * counts[m] / max(1, n_union):.1f}%)"
        p.append(f'<g data-tip="{esc(tip)}"><rect x="{x - colw / 2}" y="0" width="{colw}" height="{h}" class="hit"/>'
                 f'<path d="M{x - 9},{barh + 20} v{-(bh - 4):.1f} q0,-4 4,-4 h10 q4,0 4,4 v{bh - 4:.1f} z" class="bar"/>'
                 f'<text x="{x}" y="{barh + 14 - bh:.1f}" class="bl">{counts[m]:,}</text>')
        ys = [barh + 30 + i * rowh + rowh / 2 for i in range(len(names)) if m >> i & 1]
        if len(ys) > 1:
            p.append(f'<line x1="{x}" y1="{min(ys)}" x2="{x}" y2="{max(ys)}" class="link"/>')
        for i in range(len(names)):
            y = barh + 30 + i * rowh + rowh / 2
            p.append(f'<circle cx="{x}" cy="{y}" r="6" class="{"dot on" if m >> i & 1 else "dot"}"/>')
        p.append("</g>")
    for i, n in enumerate(names):
        y = barh + 30 + i * rowh + rowh / 2 + 4
        p.append(f'<text x="{left - 12}" y="{y}" class="rl">{esc(short(n, meta))} · b{esc(meta[n]["batch"])} ({set_sizes[n]:,})</text>')
    p.append("</svg>")
    return (f'<figure class="panel wide"><figcaption>{esc(title)}</figcaption>{"".join(p)}'
            f'<div class="mute small">{n_union:,} merged regions; bars = regions in exactly that combination; '
            f'row label = replicate · batch (its own peak count)</div></figure>')


def heatmap(order, jac, stats, meta):
    cell, left, topm = 34, 170, 70
    n = len(order)
    vals = [jac[(a, b)] for a in order for b in order if a != b]
    vmax = max(vals) if vals else 1
    w, h = left + cell * n + 60, topm + cell * n + 4
    p = [f'<svg viewBox="0 0 {w} {h}" class="heat" style="max-width:{w}px" role="img" aria-label="Peak Jaccard between replicates">']
    for j, b in enumerate(order):
        x = left + j * cell + cell / 2
        p.append(f'<text transform="translate({x + 4},{topm - 8}) rotate(-55)" class="cl">'
                 f'{esc(stats[b]["group"].replace("_", " "))} {esc(short(b, meta))}</text>')
    for i, a in enumerate(order):
        y = topm + i * cell
        p.append(f'<text x="{left - 8}" y="{y + cell / 2 + 4}" class="rl">'
                 f'{esc(stats[a]["group"].replace("_", " "))} {esc(short(a, meta))} · b{esc(meta[a]["batch"])}</text>')
        for j, b in enumerate(order):
            x = left + j * cell
            if a == b:
                p.append(f'<rect x="{x + 1}" y="{y + 1}" width="{cell - 2}" height="{cell - 2}" rx="3" class="diag"/>')
                continue
            v = jac[(a, b)]
            k = min(6, int(7 * v / vmax))
            tip = f"{a} vs {b}: Jaccard {v:.3f}"
            p.append(f'<g data-tip="{esc(tip)}"><rect x="{x + 1}" y="{y + 1}" width="{cell - 2}" height="{cell - 2}" rx="3" class="h{k}"/>'
                     f'<text x="{x + cell / 2}" y="{y + cell / 2 + 4}" class="hv t{k}">{v:.2f}</text></g>')
    p.append("</svg>")
    return "".join(p)


CSS = """
:root{color-scheme:light;--surface:#fcfcfb;--raise:#f4f3f0;--ink:#0b0b0b;--ink2:#52514e;--mute:#7a7975;--line:#e4e3df;
--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--bar:#2a78d6;--flag:#fbe9e2;
--h0:#cde2fb;--h1:#9ec5f4;--h2:#6da7ec;--h3:#3987e5;--h4:#256abf;--h5:#184f95;--h6:#0d366b;
--t0:#0b0b0b;--t1:#0b0b0b;--t2:#0b0b0b;--t3:#ffffff;--t4:#ffffff;--t5:#ffffff;--t6:#ffffff}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--surface:#1a1a19;--raise:#232322;--ink:#f0efec;--ink2:#c3c2b7;--mute:#97968f;--line:#353533;
--s1:#3987e5;--s2:#d95926;--s3:#199e70;--bar:#3987e5;--flag:#46261d;
--h0:#0d366b;--h1:#104281;--h2:#184f95;--h3:#256abf;--h4:#3987e5;--h5:#6da7ec;--h6:#b7d3f6;
--t0:#f0efec;--t1:#f0efec;--t2:#f0efec;--t3:#f0efec;--t4:#0b0b0b;--t5:#0b0b0b;--t6:#0b0b0b}}
:root[data-theme="dark"]{color-scheme:dark;--surface:#1a1a19;--raise:#232322;--ink:#f0efec;--ink2:#c3c2b7;--mute:#97968f;--line:#353533;
--s1:#3987e5;--s2:#d95926;--s3:#199e70;--bar:#3987e5;--flag:#46261d;
--h0:#0d366b;--h1:#104281;--h2:#184f95;--h3:#256abf;--h4:#3987e5;--h5:#6da7ec;--h6:#b7d3f6;
--t0:#f0efec;--t1:#f0efec;--t2:#f0efec;--t3:#f0efec;--t4:#0b0b0b;--t5:#0b0b0b;--t6:#0b0b0b}
*{box-sizing:border-box}
body{margin:0;background:var(--surface);color:var(--ink);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,sans-serif}
main{max-width:1180px;margin:0 auto;padding:28px 16px 60px}
h1{font-size:24px;margin:0 0 6px}h2{font-size:18px;margin:36px 0 8px;padding-top:12px;border-top:1px solid var(--line)}
h3{font-size:15px;margin:18px 0 6px}
p{max-width:860px;color:var(--ink2)}.mute{color:var(--mute)}.small{font-size:12px}
.scroll{overflow-x:auto}
table{border-collapse:collapse;font-variant-numeric:tabular-nums;margin:6px 0 10px}
th,td{padding:5px 9px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}
th{color:var(--ink2);font-weight:600;font-size:12px;vertical-align:bottom}
td.l,th.l{text-align:left}td.notes{text-align:left;white-space:normal;color:var(--ink2);font-size:12px;padding-top:0}
tr.noted:has(+ tr.noterow) td{border-bottom:0}
tr.noted td{background:var(--flag)}
.grid{display:flex;flex-wrap:wrap;gap:18px}
figure.panel{margin:0;padding:12px 14px;background:var(--raise);border-radius:10px;flex:1 1 340px;max-width:440px}
figure.panel.wide{max-width:none;flex:1 1 100%}
figcaption{font-weight:600;font-size:13px;margin-bottom:4px}
svg{display:block;width:100%;height:auto;max-width:100%}
.vc{fill-opacity:.16;stroke-width:2}.vc.s1{fill:var(--s1);stroke:var(--s1)}.vc.s2{fill:var(--s2);stroke:var(--s2)}.vc.s3{fill:var(--s3);stroke:var(--s3)}
.vn{text-anchor:middle;font-size:15px;font-weight:600;fill:var(--ink)}
.hit{fill:transparent}
.legend{display:flex;flex-wrap:wrap;gap:12px;font-size:12px;color:var(--ink2);margin-top:4px}
.key{display:inline-flex;align-items:center;gap:6px}.sw{width:12px;height:12px;border-radius:3px;display:inline-block}
.sw.s1{background:var(--s1)}.sw.s2{background:var(--s2)}.sw.s3{background:var(--s3)}
.upset .bar{fill:var(--bar)}.upset .bl{text-anchor:middle;font-size:11px;fill:var(--ink2)}
.axis{stroke:var(--line);stroke-width:1}.link{stroke:var(--ink);stroke-width:2}
.dot{fill:var(--line)}.dot.on{fill:var(--ink)}
.rl{text-anchor:end;font-size:12px;fill:var(--ink2)}.cl{font-size:11px;fill:var(--ink2)}
.hv{text-anchor:middle;font-size:10px}
.h0{fill:var(--h0)}.h1{fill:var(--h1)}.h2{fill:var(--h2)}.h3{fill:var(--h3)}.h4{fill:var(--h4)}.h5{fill:var(--h5)}.h6{fill:var(--h6)}
.t0{fill:var(--t0)}.t1{fill:var(--t1)}.t2{fill:var(--t2)}.t3{fill:var(--t3)}.t4{fill:var(--t4)}.t5{fill:var(--t5)}.t6{fill:var(--t6)}
.diag{fill:var(--raise);stroke:var(--line)}
.status{font-weight:600}
#tip{position:fixed;pointer-events:none;background:var(--ink);color:var(--surface);font-size:12px;padding:5px 8px;border-radius:6px;display:none;max-width:320px;z-index:9}
"""

JS = """
const tip=document.getElementById('tip');
document.addEventListener('mousemove',e=>{const t=e.target.closest('[data-tip]');
if(!t){tip.style.display='none';return}tip.textContent=t.dataset.tip;tip.style.display='block';
const x=Math.min(e.clientX+14,innerWidth-tip.offsetWidth-8);tip.style.left=x+'px';tip.style.top=(e.clientY+14)+'px'});
"""


def render(groups, stats, pair_rows, venns, jac, order, idr_pairs, qc, meta):
    tf_rows = []
    for g, reps in groups.items():
        for s in reps:
            r = stats[s]
            cls = ' class="noted"' if r["notes"] else ""
            tf_rows.append(
                f"<tr{cls}><td class='l'>{esc(g.replace('_', ' '))}</td><td class='l'>{esc(short(s, meta))}</td>"
                f"<td>{esc(r['batch'])}</td><td>{r['fragments'] / 1e6:.1f}</td><td>{100 * r['dup_rate']:.1f}%</td>"
                f"<td>{r['unique_fragments'] / 1e6:.1f}</td><td>{r['peaks']:,}</td><td>{r['median_width']}</td>"
                f"<td>{100 * r['frip']:.1f}%</td><td>{r['relaxed_peaks']:,}</td>"
                f"<td>{pct(r['top_recovered_min'])}</td><td>{pct(r['partner_top_found_min'])}</td>"
                f"<td>{pct(r.get('at_ctr_sites'))}</td>"
                f"<td>{'–' if r['self_pseudorep_idr'] is None else format(r['self_pseudorep_idr'], ',')}</td></tr>"
                + (f"<tr class='noted noterow'><td></td><td colspan='13' class='notes'>{esc(r['notes'])}</td></tr>"
                   if r["notes"] else ""))

    qc_rows = []
    for g in groups:
        q = qc.get(g)
        if not q:
            continue
        qc_rows.append(
            f"<tr><td class='l'>{esc(g.replace('_', ' '))}</td><td>{esc(q['n_reps'])}</td>"
            f"<td>{int(q['N_self_min']):,} – {int(q['N_self_max']):,}</td><td>{esc(q['self_consistency_ratio'])}</td>"
            f"<td>{'NA' if q['Nt'] == 'NA' else format(int(q['Nt']), ',')}</td>"
            f"<td>{'NA' if q['Np'] == 'NA' else format(int(q['Np']), ',')}</td>"
            f"<td>{esc(q['rescue_ratio'])}</td><td class='l status'>{esc(q['status'])}"
            # KO = ISL1 KO + ISL1 rescue, so these ratios describe real (rescued) binding.
            + (" <span class='mute small'>(KO + ISL1 rescue)</span>" if g.endswith("_KO") else "")
            + "</td></tr>")

    sections = []
    for g, reps in groups.items():
        body = [f"<h2>{esc(g.replace('_', ' '))}</h2>"]
        if len(reps) < 2:
            s = reps[0]
            body.append(f"<p>One replicate ({esc(short(s, meta))}, batch {esc(meta[s]['batch'])}): "
                        f"{stats[s]['peaks']:,} peaks, no replicate overlap possible. Its two halves "
                        f"(pseudoreplicates) share {stats[s]['self_pseudorep_idr'] or 0:,} peaks at IDR &lt; 0.05.</p>")
            sections.append("".join(body))
            continue
        v = venns[g]
        sizes = {r: stats[r]["peaks"] for r in reps}
        tops = {r: min(TOP_N, stats[r]["peaks"]) for r in reps}
        if len(reps) <= 3:
            body.append('<div class="grid">'
                        + svg_venn(*v["all"], meta, "All peaks (q 0.01)")
                        + svg_venn(*v["top"], meta, f"Top {TOP_N} peaks per replicate")
                        + "</div>")
        else:
            body.append(f"<p class='small mute'>Four replicates: shown as UpSet plots (a four-way Venn is hard to read). "
                        f"Each bar is the regions shared by exactly the replicates with filled dots.</p>")
            body.append(svg_upset(*v["all"], meta, "All peaks (q 0.01)", sizes))
            body.append(svg_upset(*v["top"], meta, f"Top {TOP_N} peaks per replicate", tops))
        prs = [p for p in pair_rows if p["group"] == g]
        body.append("<h3>Pairwise recovery</h3><div class='scroll'><table><thead><tr><th class='l'>Peaks of</th>"
                    "<th class='l'>found in</th><th>all peaks</th><th>top 1000</th></tr></thead><tbody>"
                    + "".join(f"<tr><td class='l'>{esc(short(p['a'], meta))} · b{esc(meta[p['a']]['batch'])}</td>"
                              f"<td class='l'>{esc(short(p['b'], meta))} · b{esc(meta[p['b']]['batch'])}</td>"
                              f"<td>{pct(p['a_in_b'])}</td><td>{pct(p['a_top_in_b'])}</td></tr>" for p in prs)
                    + "</tbody></table></div>")
        ip = [r for r in idr_pairs if r["group"] == g]
        if ip:
            body.append("<h3>IDR between replicates (relaxed calls, IDR &lt; 0.05)</h3><div class='scroll'><table><thead><tr>"
                        "<th class='l'>Pair</th><th>relaxed peaks</th><th>passing</th><th class='l'>used</th></tr></thead><tbody>"
                        + "".join(
                            f"<tr><td class='l'>{esc(' vs '.join(short(x, meta) for x in r['pair'].split('_vs_')))}</td>"
                            f"<td>{int(r['peaks_a']):,} / {int(r['peaks_b']):,}</td><td>{int(r['passing']):,}</td>"
                            f"<td class='l'>{'reproducible set' if r['chosen'] == 'yes' else ''}</td></tr>" for r in ip)
                        + "</tbody></table></div>")
        sections.append("".join(body))

    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>ISL1 replicate peak QC</title>
<style>{CSS}</style></head><body><main>
<h1>ISL1 CUT&amp;RUN: replicate peak QC</h1>
<p>Final calls: every ISL1 replicate against its group's merged IgG (D12_KO_rep1_IgG excluded), MACS3 BAMPE,
q 0.01, blacklist filtered, one row per interval (summits collapsed). Use this to judge which replicates agree
before choosing the reproducible and consensus peak sets. Hover a number for detail.</p>

<h2>Per replicate</h2>
<div class="scroll"><table><thead><tr><th class="l">Group</th><th class="l">Rep</th><th>Batch</th>
<th>Fragments<br>(M)</th><th>Dup<br>rate</th><th>Unique<br>(M)</th><th>Peaks<br>q 0.01</th><th>Median<br>width</th>
<th>FRiP</th><th>Relaxed<br>p 0.01</th><th>Own top 1000<br>found in partners<br>(min)</th>
<th>Partners' top 1000<br>found in it<br>(min)</th><th>KO peaks at<br>Ctr ISL1 sites</th><th>Self-pseudorep<br>IDR peaks</th></tr></thead><tbody>{''.join(tf_rows)}</tbody></table></div>
<p class="small mute">Fragments: read pairs after MAPQ 20 / proper-pair filtering (duplicates kept, as MACS3 sees them).
FRiP: share of those fragments inside the replicate's own peaks. Relaxed: the p 0.01 calls used for IDR, capped at
300,000; a library at the cap has much more weak signal. "Partners' top 1000 found in it" is the lowest share of
another replicate's strongest 1000 peaks that this replicate also calls: low values mean it misses strong sites the
others see. Replicates with notes are shaded, the note on the line below (top-peak notes for Ctr groups only).</p>

<h2>ENCODE IDR reproducibility</h2>
<div class="scroll"><table><thead><tr><th class="l">Group</th><th>Replicates</th><th>Self-pseudorep IDR peaks<br>(min – max over replicates)</th>
<th>Self-consistency<br>ratio</th><th>Nt<br>(true replicates)</th><th>Np<br>(pooled pseudoreps)</th><th>Rescue<br>ratio</th>
<th class="l">ENCODE</th></tr></thead><tbody>{''.join(qc_rows)}</tbody></table></div>
<p class="small mute">Self-consistency = largest / smallest self-pseudoreplicate IDR count across a group's replicates (do the
replicates have similar internal quality?). Rescue = max(Nt, Np) / min(Nt, Np) (does pooling recover much more than the
replicates agree on?). Both &le; 2 pass, one &gt; 2 borderline, both &gt; 2 fail. Nt is the best replicate pair, as used for the
reproducible set. All IDR on relaxed calls, ranked by p value, global IDR &lt; 0.05.</p>

{''.join(sections)}

<h2>All replicates: peak overlap (Jaccard, bp)</h2>
<p>Base-pair Jaccard of every pair of replicates' peak sets, across groups. Batch structure shows up here as blocks.</p>
<div class="scroll heatwrap">{heatmap(order, jac, stats, meta)}</div>
<details><summary class="small">Table view</summary><div class="scroll"><table><thead><tr><th class="l"></th>
{''.join(f'<th>{esc(stats[b]["group"].replace("_", " "))} {esc(short(b, meta))}</th>' for b in order)}</tr></thead><tbody>
{''.join('<tr><th class="l">' + esc(stats[a]['group'].replace('_', ' ')) + ' ' + esc(short(a, meta)) + '</th>' + ''.join('<td>' + ('–' if a == b else f'{jac[(a, b)]:.3f}') + '</td>' for b in order) + '</tr>' for a in order)}
</tbody></table></div></details>
<p class="small mute">Generated by scripts/replicate_peak_report.py (scripts/05_replicate_peak_report.sb).</p>
</main><div id="tip"></div><script>{JS}</script></body></html>
"""


if __name__ == "__main__":
    if len(sys.argv) != 10:
        sys.exit(__doc__)
    main(*sys.argv[1:])
