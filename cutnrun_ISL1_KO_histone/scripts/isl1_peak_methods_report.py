#!/usr/bin/env python3
"""ISL1 CUT&RUN peak-calling methods report (control ISL1 only), HTML. Standard library only.

Reads the outputs of steps 03/13 (per-library calls), 05 (batch 3 rep1/rep2 sets), 15 (ctr_trio
three-replicate sets) and 18 (kept sets); uses bedtools for overlaps.
  python isl1_peak_methods_report.py <OUT_DIR> <STUDY_DIR> <bedtools> <out.html>
"""
import csv
import datetime
import html
import os
import statistics
import subprocess
import sys

OUT, STUDY, BT, HTML_OUT = sys.argv[1:5]
TRIO = {"D15": [1, 3, 4], "D18": [1, 3, 4], "D22": [2, 4, 5]}
REF = os.path.join(OUT, "isl1_consensus", "reference.bed")
FINAL = os.path.join(OUT, "isl1_ctr_final")


def tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def nlines(path):
    with open(path) as f:
        return sum(1 for line in f if line.strip() and not line.startswith(("#", "track")))


def widths(path):
    w = []
    with open(path) as f:
        for line in f:
            p = line.split("\t")
            if len(p) >= 3 and not line.startswith(("#", "track")):
                w.append(int(p[2]) - int(p[1]))
    return w


def n_overlapping(a, b):
    r = subprocess.run([BT, "intersect", "-u", "-a", a, "-b", b], capture_output=True, text=True, check=True)
    return sum(1 for line in r.stdout.splitlines() if line.strip())


def num(x, dec=0):
    if x in (None, "", "NA"):
        return "–"
    try:
        v = float(x)
    except ValueError:
        return html.escape(str(x))
    if v != v:
        return "–"
    return f"{v:,.{dec}f}"


def pct(a, b, dec=0):
    return "–" if not b else f"{100 * a / b:.{dec}f}%"


def table(head, rows, num_cols=()):
    h = "".join(f'<th class="{"n" if i in num_cols else ""}">{html.escape(c)}</th>' for i, c in enumerate(head))
    body = []
    for r in rows:
        body.append("<tr>" + "".join(f'<td class="{"n" if i in num_cols else ""}">{c}</td>' for i, c in enumerate(r)) + "</tr>")
    return f'<div class="tw"><table><thead><tr>{h}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


sections = []

# ---- 1. methods overview
sections.append(("Methods in this run", "<p>All calls use each library's own deduplicated block IgG "
    "(matched IgG), MACS3, blacklist filtered. Only control ISL1 libraries are used for the consensus "
    "sets; KO lines have no ISL1 to call.</p>" + table(
    ["Method", "Rule", "Replicates", "Chance estimate", "Step / output"],
    [["Per-library MACS3", "narrow, q 0.01, summits; regions = distinct peaks", "each library", "–", "03, 13 · <code>macs3/</code>, <code>ctr_matched/</code>"],
     ["ENCODE overlap (Nt / Np)", "pooled relaxed call (p 0.01); keep peaks overlapping both replicate (Nt) or both pseudo-pool (Np) lists by ≥ 50% of either peak", "batch 3 rep1 + rep2", "–", "05 · <code>reproducibility/overlap/</code>"],
     ["IDR (Nt / Np)", "IDR &lt; 0.05 on relaxed calls, pooled peaks as <code>--peak-list</code>", "batch 3 rep1 + rep2", "–", "05 · <code>reproducibility/idr/</code>"],
     ["A: MACS3 replicate threshold", "merge per-library q 0.01 regions; keep regions with peaks from ≥ k replicates (nf-core <code>replicate_threshold</code>)", "chosen 3 per day", "each replicate shuffled independently, counted alike", "15 · <code>ctr_trio/&lt;day&gt;/A/</code>"],
     ["B: ENCODE overlap, ≥ k", "pooled relaxed peaks overlapping the relaxed peaks of ≥ k replicates by ≥ 50%", "chosen 3 per day", "–", "15 · <code>ctr_trio/&lt;day&gt;/B/</code>"],
     ["IDR, ≥ k", "IDR on every replicate pair; k = 2: passes in ≥ 1 pair; k = 3: passes in pairs covering all 3 replicates", "chosen 3 per day", "–", "15 · <code>ctr_trio/&lt;day&gt;/idr/</code>"],
     ["IDR, best pair", "the single pair with most passing peaks (ENCODE \"optimal\")", "2 of the chosen 3", "–", "15"],
     ["Kept: macs3_ge2", "= A, k = 2", "chosen 3 per day", "as A", "18 · <code>isl1_ctr_final/macs3_ge2/</code>"],
     ["Kept: idr_any_pair", "= IDR, k = 2", "chosen 3 per day", "–", "18 · <code>isl1_ctr_final/idr_any_pair/</code>"]])
    + "<p class=\"note\">Not part of this run (exploratory, not repeated on the 10x-filtered data): SEACR, "
    "the MACS3 rep1/rep2 pair consensus (step 07), the all-replicate sets (<code>ctr_multi</code>) and "
    "the group-IgG test.</p>"
    "<p class=\"note\"><b>Reference</b> = the earlier batches' D15 Ctr conservative ISL1 set (9,329 peaks), "
    "QC only. It was built from the library that is now D18 rep3, so it is not independent of the earlier "
    "libraries.</p>"))

# ---- 2. per-library calls
reps = tsv(os.path.join(STUDY, "study_config", "ctr_replicates.tsv"))
rows = []
for r in sorted(reps, key=lambda r: (int(r["day"][1:]), int(r["replicate"]))):
    day, rep = r["day"], int(r["replicate"])
    st_path = (os.path.join(OUT, "qc", "peaks", f"{r['original_sample']}.peak_stats.tsv") if r["batch_set"] == "batch3"
               else os.path.join(OUT, "ctr_matched", "qc", f"{r['sample']}.peak_stats.tsv"))
    st = tsv(st_path)[0] if os.path.exists(st_path) else {}
    src = "this study (10x-filtered)" if "ISL1_KO_histone_cutNrun/results/bam" in r["isl1_bam"] else "ISL1_cutNrun"
    used = "yes" if rep in TRIO.get(day, []) else ""
    frip = st.get("frip")
    rows.append([day, f"rep{rep}", r["batch_set"], html.escape(r["original_sample"]), src,
                 num(int(st["fragments"]) / 1e6, 1) + " M" if st.get("fragments") else "–",
                 num(st.get("filtered_regions", st.get("regions"))), f"{100 * float(frip):.1f}%" if frip else "–", used])
sections.append(("Control ISL1 libraries: per-library calls",
    table(["Day", "Rep", "Set", "Original name", "BAM from", "Fragments", "Peak regions", "FRiP", "In trio"],
          rows, num_cols=(5, 6, 7))))

# ---- 3. batch 3 rep1/rep2 (step 05)
rows = []
for day in ["D9", "D15", "D18", "D22"]:
    q = tsv(os.path.join(OUT, "reproducibility", "qc", f"{day}_H3CTR_ISL1.repro_qc.tsv"))
    for m in q:
        rows.append([day, m["method"], num(m["N1"]), num(m["N2"]), f"<b>{num(m['Nt'])}</b>", num(m["Np"]),
                     m["self_consistency"], m["rescue"], m["encode_flag"], num(m["Nt_Mb"], 1),
                     f"{100 * float(m['Nt_frip_rep1']):.1f} / {100 * float(m['Nt_frip_rep2']):.1f}%"])
sections.append(("Batch 3 rep1 + rep2 only (ENCODE design, step 05)",
    "<p>N1 / N2: each replicate against its own pseudoreplicate halves. Self = N1/N2, rescue = Np/Nt; "
    "both ≤ 2 = pass.</p>" +
    table(["Day", "Method", "N1", "N2", "Nt", "Np", "Self", "Rescue", "Flag", "Nt Mb", "Nt FRiP rep1 / rep2"],
          rows, num_cols=(2, 3, 4, 5, 6, 7, 9, 10))))

# ---- 4. trio sets (step 15)
label = {"A": "A: MACS3 ≥ k", "B": "B: overlap ≥ k", "IDR": "IDR ≥ k"}
rows = []
for day, keep in TRIO.items():
    d = os.path.join(OUT, "ctr_trio", day)
    for s in tsv(os.path.join(d, "summary.tsv")):
        m, k = s["method"], s["k"]
        bed = (os.path.join(d, {"A": "A", "B": "B", "IDR": "idr"}[m], f"k{k}.bed") if m in label
               else os.path.join(d, "idr", "optimal.bed"))
        w = widths(bed) if os.path.exists(bed) else []
        name = label.get(m, "IDR best pair (" + m.replace("IDR_optimal_", "").replace("_", "–") + ")")
        rows.append([day, ",".join(map(str, keep)), name, k, f"<b>{num(s['regions'])}</b>", num(s["chance"]),
                     num(s["over_chance"], 1), num(s["Mb"], 2), num(statistics.median(w)) if w else "–",
                     num(s["cross_batch"]), num(s["earlier_only"]),
                     f"{num(s['top1000_in_reference_pct'], 1)} / {num(s['all_in_reference_pct'], 1)}%"])
sections.append(("Three chosen control replicates per day (ctr_trio, step 15)",
    "<p>Chosen replicates: the better batch 3 library plus the earlier libraries (D15 and D18: rep1, 3, 4; "
    "D22: rep2, 4, 5). Cross-batch = region supported by the batch 3 library and an earlier library; "
    "no region in any set is supported by the batch 3 library alone.</p>" +
    table(["Day", "Reps", "Method", "k", "Regions", "Chance", "× chance", "Mb", "Median bp",
           "Cross-batch", "Earlier only", "In reference: top 1000 / all"],
          rows, num_cols=(3, 4, 5, 6, 7, 8, 9, 10, 11))))

rows = []
for day in TRIO:
    for p in tsv(os.path.join(OUT, "ctr_trio", day, "idr", "pair_counts.tsv")):
        rows.append([day, f"rep{p['rep_a']} – rep{p['rep_b']}", f"{p['batch_a']} × {p['batch_b']}", num(p["passing_peaks"])])
sections.append(("IDR by replicate pair", table(["Day", "Pair", "Batches", "Peaks passing IDR < 0.05"], rows, num_cols=(3,))))

# ---- 5. kept sets (step 18) and how they overlap
rows = []
for day in TRIO:
    a = os.path.join(FINAL, "macs3_ge2", f"{day}_Ctr_ISL1.bed")
    b = os.path.join(FINAL, "idr_any_pair", f"{day}_Ctr_ISL1.bed")
    na, nb = nlines(a), nlines(b)
    a_in_b, b_in_a = n_overlapping(a, b), n_overlapping(b, a)
    wa, wb = widths(a), widths(b)
    ra, rb = n_overlapping(a, REF), n_overlapping(b, REF)
    rows.append([day, f"<b>{num(na)}</b>", num(statistics.median(wa)), num(sum(wa) / 1e6, 2), pct(ra, na),
                 f"<b>{num(nb)}</b>", num(statistics.median(wb)), num(sum(wb) / 1e6, 2), pct(rb, nb),
                 f"{num(a_in_b)} ({pct(a_in_b, na)})", f"{num(b_in_a)} ({pct(b_in_a, nb)})"])
sections.append(("Kept sets (step 18) and their overlap",
    "<p>Files: <code>results/isl1_ctr_final/{macs3_ge2,idr_any_pair}/&lt;day&gt;_Ctr_ISL1.bed</code>. "
    "Columns: chr, start, end, name, n_reps, reps, max −log10 p, n_batch3, n_earlier.</p>" +
    table(["Day", "macs3_ge2", "Median bp", "Mb", "In reference", "idr_any_pair", "Median bp", "Mb", "In reference",
           "macs3_ge2 in IDR", "IDR in macs3_ge2"], rows, num_cols=tuple(range(1, 11)))))

# ---- write
css = """
:root{--bg:#fff;--fg:#1d1d1f;--mut:#5f6368;--line:#dcdcdc;--head:#f3f4f6;--acc:#1f5fa8}
@media (prefers-color-scheme:dark){:root{--bg:#16181b;--fg:#e8e8e8;--mut:#a0a4aa;--line:#33363b;--head:#22252a;--acc:#7fb0ea}}
body{background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;margin:0 auto;max-width:1200px;padding:24px 16px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 8px;color:var(--acc)}
p{margin:6px 0 10px}.note{color:var(--mut);font-size:13px}.sub{color:var(--mut);margin-bottom:16px}
.tw{overflow-x:auto}table{border-collapse:collapse;margin:8px 0 4px;font-size:13px}
th,td{border-bottom:1px solid var(--line);padding:5px 10px;text-align:left;white-space:nowrap}
th{background:var(--head);font-weight:600}td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
code{font-size:12px}
"""
stamp = datetime.date.today().isoformat()
body = "".join(f"<h2>{html.escape(t)}</h2>{c}" for t, c in sections)
os.makedirs(os.path.dirname(HTML_OUT), exist_ok=True)
with open(HTML_OUT, "w") as f:
    f.write(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>ISL1 peak methods</title>
<style>{css}</style></head><body>
<h1>ISL1 CUT&amp;RUN: peak-calling methods</h1>
<div class="sub">Control ISL1, D9-D22 · batch 3 libraries on the 10x-filtered reads · built {stamp}</div>
{body}</body></html>""")
print(f"wrote {HTML_OUT}")
