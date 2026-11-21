"""Report for step 7.3i (alpha vs beta bias of the progenitors), built on the Mac.

Reads the 7.3i outputs copied from Alpine into results/trajectory/alpha_bias/{D15,D15_D22}/:
  scp -r -i ~/.ssh/id_rsa \\
    "hzhu10@xsede.org@login.rc.colorado.edu:/scratch/alpine/hzhu10@xsede.org/ISL1_scRNAseq/out_all/day_analysis/D15/trajectory/alpha_bias/D15" \\
    "hzhu10@xsede.org@login.rc.colorado.edu:/scratch/alpine/hzhu10@xsede.org/ISL1_scRNAseq/out_all/day_analysis/D15/trajectory/alpha_bias/D15_D22" \\
    results/trajectory/alpha_bias/
Writes results/trajectory/alpha_bias/index.html (both datasets; pseudotime / velocity kernel switch).
Run from scRNAseq/: ~/anaconda3/bin/python3 analysis/scripts/build_alpha_bias_report.py
"""
import html
import sys
from pathlib import Path

import pandas as pd

R = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2] / "results" / "trajectory" / "alpha_bias"
SETS = [("D15", "D15 only",
         "7.3d CellRank fates (run on all cells). α = pre_alpha, β = pre_beta2, EC = pre_EC. pre_EC is an end state here, "
         "so pre_EC cells have no bias and are not included."),
        ("D15_D22", "D15 + D22 (alpha_beta_EC subset)",
         "7.4d CellRank fates (run on all cells). α = SC_alpha, β = SC_beta, EC = SC_EC1 + SC_EC2 (D22 end states), so the "
         "D15 pre_EC cells are not end states and get a bias.")]


def table(df, pct=()):
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        tds = []
        for c, v in r.items():
            if isinstance(v, float):
                s = f"{v:+.2f}" if c in pct else f"{v:.3f}"
                tds.append(f"<td class=n>{s}</td>")
            elif isinstance(v, (int,)) or str(v).isdigit():
                tds.append(f"<td class=n>{v}</td>")
            else:
                tds.append(f"<td>{html.escape(str(v))}</td>")
        rows.append("<tr>" + "".join(tds) + "</tr>")
    return f"<div class=tw><table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"


sec = []
for key, title, note in SETS:
    ch = pd.read_csv(R / key / "check1.csv")
    co = pd.read_csv(R / key / "correlations.csv")
    co["y"] = co.y.map({"P_EC": "P(EC)", "pseudotime": "pseudotime"})
    blocks = []
    for k in ["pseudotime", "velocity"]:
        c1 = ch[ch.kernel == k].drop(columns="kernel")
        cc = co[co.kernel == k].drop(columns="kernel")
        figs = "".join(f"<figure><img src=\"{key}/{f}_{k}.png\" alt=\"{html.escape(cap)}\" loading=lazy>"
                       f"<figcaption>{html.escape(cap)}</figcaption></figure>"
                       for f, cap in [("umap", "α bias on the 2D UMAP, control | KO (× = P(α) + P(β) < 0.05)"),
                                      ("hist", "Distribution of the α bias per subcluster, control (filled) vs KO (line)"),
                                      ("vs_ec", "α bias vs P(EC)"),
                                      ("vs_pt", "α bias vs single-value pseudotime")])
        blocks.append(f"<div class=\"kern\" data-k=\"{k}\"{'' if k == 'pseudotime' else ' hidden'}>"
                      f"<h3>Check 1: is there a ratio for every cell?</h3>"
                      f"<p class=mut>S = P(α) + P(β). S = 0: no ratio; S &lt; 0.05: a ratio exists but rests on small "
                      f"probabilities. Values are medians unless noted.</p>{table(c1)}"
                      f"{figs}<h3>Spearman correlations</h3>{table(cc, pct=set(cc.columns) - {'y', 'genotype'})}</div>")
    sec.append(f"<section><h2>{html.escape(title)}</h2><p class=sub>{html.escape(note)}</p>{''.join(blocks)}</section>")

page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Progenitor Alpha Bias</title>
<style>
:root{{color-scheme:light;--bg:#fcfcfb;--card:#fff;--fg:#0b0b0b;--fg2:#52514e;--mut:#7a7975;--line:#e4e3de}}
@media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{color-scheme:dark;--bg:#1a1a19;--card:#232322;--fg:#f0efec;--fg2:#c3c2b7;--mut:#94938c;--line:#383835}}}}
:root[data-theme="dark"]{{color-scheme:dark;--bg:#1a1a19;--card:#232322;--fg:#f0efec;--fg2:#c3c2b7;--mut:#94938c;--line:#383835}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}}
main{{max-width:1500px;margin:0 auto;padding:24px 16px 56px}}
h1{{font-size:24px;margin:0 0 6px}} h2{{font-size:20px;margin:36px 0 4px}} h3{{font-size:15px;margin:22px 0 6px}}
p{{margin:4px 0}} .sub{{color:var(--fg2)}} .mut{{color:var(--mut);font-size:13px}}
.ctl{{position:sticky;top:0;background:var(--bg);padding:10px 0;z-index:1;display:flex;gap:10px;align-items:center;border-bottom:1px solid var(--line)}}
.seg{{display:inline-flex;border:1px solid var(--line);border-radius:7px;overflow:hidden}}
.seg button{{font:inherit;font-size:13px;padding:5px 12px;border:0;background:var(--card);color:var(--fg2);cursor:pointer}}
.seg button[aria-pressed="true"]{{background:var(--fg);color:var(--bg)}}
figure{{margin:14px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px}}
figure img{{width:100%;height:auto;display:block;background:#fff;border-radius:6px}}
figcaption{{color:var(--fg2);font-size:13px;margin-top:6px}}
.tw{{overflow-x:auto}} table{{border-collapse:collapse;font-size:12.5px;width:100%;background:var(--card)}}
th,td{{padding:5px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}
th{{color:var(--fg2);font-weight:600}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
</style></head><body><main>
<h1>α vs β bias of the endocrine progenitors</h1>
<p class=sub>Per cell, α bias = P(α) / (P(α) + P(β)) from the existing CellRank fate probabilities (nothing rerun): 0 = β-leaning,
1 = α-leaning. Cells: ENP1, ENP2, ENP_FEV_high, pre_alpha_KO (and pre_EC for D15 + D22). Single-value pseudotime =
Slingshot pseudotimes averaged with the curve weights, scaled 0-1 (as used by the CellRank pseudotime kernel). Velocity
runs backwards on these data, so the pseudotime kernel is the main one and the velocity kernel a check.</p>
<div class=ctl><span class=mut>Kernel</span><span class=seg role=group aria-label="Kernel">
<button data-k=pseudotime aria-pressed=true>Pseudotime</button><button data-k=velocity aria-pressed=false>Velocity</button></span></div>
{''.join(sec)}
</main><script>
document.querySelectorAll('.seg button').forEach(b => b.onclick = () => {{
  document.querySelectorAll('.seg button').forEach(x => x.setAttribute('aria-pressed', x === b));
  document.querySelectorAll('.kern').forEach(d => d.hidden = d.dataset.k !== b.dataset.k);
}});
</script></body></html>"""
(R / "index.html").write_text(page)
print("wrote", R / "index.html")
