#!/usr/bin/env python3
"""Per ISL1 sample: compare MACS3 peaks called against its own IgG (matched) vs the pooled IgG.

A sample whose two peak sets disagree much more than the rest of the cohort points at its own
IgG. Writes igg_method_comparison.tsv and igg_method_comparison.html to <out_dir>.

usage: compare_igg_methods.py <samples.tsv> <sample_metadata.tsv> <matched_macs3_dir>
                              <pooled_macs3_dir> <pool_manifest.tsv> <bedtools> <out_dir>
"""
import csv
import html
import math
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

TOP_N = 1000          # strongest peaks (by -log10 q) checked for recovery in the other set
# Flag a sample when either
#   min_overlap < MIN_OVERLAP: under this share of one peak set is found in the other, or
#   log2(n_pooled / n_matched) is more than OUTLIER_MADS scaled MADs from the cohort median.
# A median/MAD cutoff on min_overlap itself does not work here: Ctr samples (thousands of
# peaks) and KO samples (hundreds) spread the cohort so widely that the cutoff goes negative.
MIN_OVERLAP = 0.10
OUTLIER_MADS = 3


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


def write_bed(keys, path):
    with open(path, "w") as f:
        for c, s, e in sorted(keys):
            f.write(f"{c}\t{s}\t{e}\n")


def bedtools(tool, *args):
    return subprocess.run([tool, *args], capture_output=True, text=True, check=True).stdout


def n_hit(tool, a, b):
    return sum(1 for l in bedtools(tool, "intersect", "-u", "-a", a, "-b", b).splitlines() if l)


def jaccard(tool, a, b):
    hdr, val = bedtools(tool, "jaccard", "-a", a, "-b", b).splitlines()[:2]
    return float(val.split()[hdr.split().index("jaccard")])


def frac(n, d):
    return n / d if d else float("nan")


def main(samples_tsv, metadata_tsv, matched_dir, pooled_dir, manifest_tsv, tool, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    meta = {r["sample_name"]: r for r in read_tsv(metadata_tsv)}
    manifest = {r["control"]: r for r in read_tsv(manifest_tsv) if r["control"] != "pooled_total"}

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for s in read_tsv(samples_tsv):
            if s["role"] != "TF":
                continue
            name, igg = s["sample_name"], s["control_sample"]
            m = unique_peaks(Path(matched_dir) / f"{name}_peaks.filtered.narrowPeak")
            p = unique_peaks(Path(pooled_dir) / f"{name}_peaks.filtered.narrowPeak")
            beds = {}
            for tag, peaks in (("m", m), ("p", p), ("mtop", sorted(m, key=m.get, reverse=True)[:TOP_N]),
                               ("ptop", sorted(p, key=p.get, reverse=True)[:TOP_N])):
                beds[tag] = f"{tmp}/{name}.{tag}.bed"
                write_bed(peaks, beds[tag])
            ig = manifest.get(igg, {})
            n_in, n_dd = int(ig.get("fragments_input", 0)), int(ig.get("fragments_deduplicated", 0))
            row = {
                "sample": name, "day": meta[name]["day"], "genotype": meta[name]["genotype"],
                "batch": meta[name]["batch"], "igg": igg,
                "igg_fragments": n_in, "igg_unique_fragments": n_dd,
                "igg_dup_rate": 1 - frac(n_dd, n_in),
                "n_matched": len(m), "n_pooled": len(p),
                "matched_in_pooled": frac(n_hit(tool, beds["m"], beds["p"]), len(m)),
                "pooled_in_matched": frac(n_hit(tool, beds["p"], beds["m"]), len(p)),
                "top_matched_in_pooled": frac(n_hit(tool, beds["mtop"], beds["p"]), min(TOP_N, len(m))),
                "top_pooled_in_matched": frac(n_hit(tool, beds["ptop"], beds["m"]), min(TOP_N, len(p))),
                "jaccard_bp": jaccard(tool, beds["m"], beds["p"]) if m and p else float("nan"),
            }
            row["min_overlap"] = min(row["matched_in_pooled"], row["pooled_in_matched"])
            row["log2_peak_ratio"] = (math.log2(len(p) / len(m)) if m and p else float("nan"))
            rows.append(row)

    vals = [r["log2_peak_ratio"] for r in rows if r["log2_peak_ratio"] == r["log2_peak_ratio"]]
    med = statistics.median(vals)
    mad = 1.4826 * statistics.median(abs(v - med) for v in vals)
    for r in rows:
        why = []
        if r["min_overlap"] < MIN_OVERLAP:
            why.append("low overlap")
        if mad and abs(r["log2_peak_ratio"] - med) > OUTLIER_MADS * mad:
            why.append("peak-count ratio")
        r["flag"] = "; ".join(why)
    rows.sort(key=lambda r: r["min_overlap"])

    cols = list(rows[0])
    with open(out / "igg_method_comparison.tsv", "w") as f:
        w = csv.DictWriter(f, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (out / "igg_method_comparison.html").write_text(render_html(rows, med, mad))
    print(f"{len(rows)} samples; log2(pooled/matched peaks) median {med:.2f}, scaled MAD {mad:.2f}; "
          f"flagged: {[r['sample'] for r in rows if r['flag']] or 'none'}")


def render_html(rows, med, mad):
    pct = lambda v: "n/a" if v != v else f"{100 * v:.1f}%"
    head = ["Sample", "Day", "Geno", "Batch", "Own IgG unique frags", "IgG dup", "Peaks matched",
            "Peaks pooled", "log2 pooled/matched", "Matched found in pooled", "Pooled found in matched",
            f"Top {TOP_N} matched in pooled", f"Top {TOP_N} pooled in matched", "Jaccard (bp)", ""]
    body = []
    for r in rows:
        cells = [r["sample"], r["day"], r["genotype"], r["batch"], f"{r['igg_unique_fragments'] / 1e6:.1f} M",
                 pct(r["igg_dup_rate"]), f"{r['n_matched']:,}", f"{r['n_pooled']:,}",
                 f"{r['log2_peak_ratio']:.2f}",
                 pct(r["matched_in_pooled"]), pct(r["pooled_in_matched"]),
                 pct(r["top_matched_in_pooled"]), pct(r["top_pooled_in_matched"]),
                 f"{r['jaccard_bp']:.3f}", r["flag"]]
        cls = ' class="flag"' if r["flag"] else ""
        body.append(f"<tr{cls}>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in cells) + "</tr>")
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>IgG method comparison</title>
<style>
:root{{--bg:#fff;--fg:#1d1d1f;--mute:#6b6b70;--line:#e3e3e6;--flag:#fde8e4}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1b1b1d;--fg:#ececee;--mute:#9a9aa0;--line:#333336;--flag:#4a2622}}}}
body{{background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,system-ui,sans-serif;margin:24px;max-width:1200px}}
table{{border-collapse:collapse;font-variant-numeric:tabular-nums}}
th,td{{padding:6px 10px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}}
th:nth-child(-n+4),td:nth-child(-n+4){{text-align:left}}
th{{color:var(--mute);font-weight:600;font-size:12px}} tr.flag td{{background:var(--flag)}}
p{{max-width:820px}} .mute{{color:var(--mute)}}
</style></head><body>
<h1>ISL1 CUT&amp;RUN: own IgG vs pooled IgG</h1>
<p>Each ISL1 sample was called with MACS3 twice: against its own deduplicated IgG (matched) and against
one pool of all 13 IgGs, each deduplicated and downsampled to the same depth (pooled). Rows are sorted by
the smaller of the two overlap percentages. A row is flagged when that value is under
{100 * MIN_OVERLAP:.0f}% (low overlap), or when log2(pooled / matched peak count) is more than
{OUTLIER_MADS} scaled MADs from the cohort median ({med:.2f}; MAD {mad:.2f}; peak-count ratio).</p>
<p><b>Reading the direction.</b> Matched peaks missing from pooled: the own IgG lacks background the pool has
(too shallow or too flat), so extra peaks appear at sticky or open regions. Pooled peaks missing from matched:
the own IgG has signal at real sites (carry-over), or it captures sample-specific background that the pool
dilutes. KO samples have few true peaks, so their percentages rest on small counts.</p>
<table><thead><tr>{''.join(f'<th>{h}</th>' for h in head)}</tr></thead><tbody>
{''.join(body)}
</tbody></table>
<p class="mute">Peaks: blacklist-filtered narrowPeak, one row per interval (summits collapsed). Overlap = any
1 bp. Top {TOP_N} by -log10 q. Unique fragments and duplicate rate from the pool manifest
(02_build_igg_background.sb, pooled).</p>
</body></html>
"""


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    main(*sys.argv[1:])
