"""Step 7.1e: one-page report for a 7.1d compartment, built around an interactive 3D UMAP.

Request (user, 2026-11-18): build the 3D UMAP and a report for the D15 endocrine subclusters.

Inputs (7.1d, ${OUT_DIR}/day_analysis/{day}/subcluster_{name}/):
  {day}_{name}.h5ad (obs: parent, genotype, sample, sub_r*; obsm X_umap_3d; X = unregressed log1p),
  composition_r{res}.csv, markers_r{res}.csv, dotplot_genes_r{res}.png; obsm X_umap_2d for the 2D panels
  + labels: study_config/annotation/{day}_{name}_r{res}.tsv (subcluster, label, group, status) -- optional;
    without it the page shows subcluster numbers only.
Writes: {day}_{name}_r{res}_report.html in the same folder (self-contained except plotly.js from cdnjs).
The 3D and 2D coordinates are stored once and recoloured in the browser (7.1d's umap3d.html repeats them per gene).
2D UMAP (user, 2026-11-18): interactive, larger points, panels all / control / KO with the other genotype in grey.
Run: sbatch analysis/scripts/7.1e_subcluster_report.sb D15 endocrine 0.6
"""
import base64
import html
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = Path(os.environ["OUT_DIR"])
GENES = ["CHGA", "NEUROG3", "NKX2-2", "NEUROD1", "INSM1", "PAX4", "ARX", "IRX2", "PAX6", "ISL1", "MAFB", "MAFA",
         "INS", "GCG", "SST", "GHRL", "PPY", "PDX1", "NKX6-1", "FEV", "LMX1A", "TPH1", "SLC18A1", "DDC",
         "PHOX2A", "SCGN", "FOXP2", "TOX", "DCC", "MKI67"]
SAMPLE_ORDER = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def main():
    day, name, res = sys.argv[1], sys.argv[2], float(sys.argv[3])
    out = OUT_DIR / "day_analysis" / day / f"subcluster_{name}"
    rk = f"sub_r{str(res).replace('.', 'p')}"
    a = sc.read_h5ad(out / f"{day}_{name}.h5ad")
    comp = pd.read_csv(out / f"composition_r{res}.csv", index_col=0)
    mk = pd.read_csv(out / f"markers_r{res}.csv")
    par = json.loads((out / "parent.json").read_text())

    lab_f = REPO / "study_config" / "annotation" / f"{day}_{name}_r{res}.tsv"
    cl_ids = sorted(comp.index.astype(int))
    if lab_f.exists():
        lab = pd.read_csv(lab_f, sep="\t", index_col=0)
        lab.index = lab.index.astype(int)
        assert set(lab.index) == set(cl_ids), f"{lab_f} does not list every subcluster at r{res}"
        draft = (lab["status"] == "draft").any()
    else:
        lab = pd.DataFrame({"label": [f"subcluster {c}" for c in cl_ids], "group": "unlabelled",
                            "status": "none"}, index=cl_ids)
        draft = False
    groups = list(dict.fromkeys(lab.loc[[c for c in lab.index], "group"]))   # file order = colour order

    # ---- per-cell data, stored once
    u = a.obsm["X_umap_3d"]; u2 = a.obsm["X_umap_2d"]
    sub = a.obs[rk].astype(int).to_numpy()
    genes = [g for g in GENES if g in a.var_names]
    X = a[:, genes].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    data = {
        "x": np.round(u[:, 0], 2).tolist(), "y": np.round(u[:, 1], 2).tolist(), "z": np.round(u[:, 2], 2).tolist(),
        "u": np.round(u2[:, 0], 2).tolist(), "v": np.round(u2[:, 1], 2).tolist(),
        "sub": sub.tolist(),
        "geno": (a.obs["genotype"].astype(str) == "KO").astype(int).tolist(),
        "sample": [SAMPLE_ORDER.index(s) for s in a.obs["sample"].astype(str)],
        "parent": a.obs["parent"].astype(int).tolist(),
        "genes": genes,
        "expr": {g: np.round(X[:, i], 1).tolist() for i, g in enumerate(genes)},
        "clusters": [{"id": int(c), "label": str(lab.loc[c, "label"]), "group": str(lab.loc[c, "group"])}
                     for c in cl_ids],
        "groups": groups,
        "samples": SAMPLE_ORDER,
        "parents": sorted(int(p) for p in a.obs["parent"].astype(int).unique()),
    }

    # ---- summary table
    top = (mk[mk["logfoldchanges"] > 0].groupby("group").head(8).groupby("group")["names"]
           .apply(lambda v: ", ".join(v)))
    top.index = top.index.astype(int)
    rows = []
    for c in cl_ids:
        r = comp.loc[c]
        rows.append({"id": c, "label": lab.loc[c, "label"], "group": lab.loc[c, "group"], "cells": int(r["cells"]),
                     "ctr": float(r["pct_of_ctr_subset"]), "ko": float(r["pct_of_KO_subset"]),
                     "markers": top.get(c, "")})
    order = sorted(rows, key=lambda d: (groups.index(d["group"]), d["id"]))
    data["comp"] = order
    n_ctr = int((a.obs["genotype"] == "ctr").sum()); n_ko = int((a.obs["genotype"] == "KO").sum())
    hormone = [d for d in rows if d["group"].startswith("hormone")]
    koonly = [d for d in rows if d["group"] == "KO-only"]
    kpi = {
        "cells": f"{a.n_obs:,}", "ctr": f"{n_ctr:,}", "ko": f"{n_ko:,}",
        "hormone_ctr": round(sum(d["ctr"] for d in hormone), 1), "hormone_ko": round(sum(d["ko"] for d in hormone), 1),
        "koonly_ctr": round(sum(d["ctr"] for d in koonly), 1), "koonly_ko": round(sum(d["ko"] for d in koonly), 1),
    }

    def esc(s):
        return html.escape(str(s))
    trs = "".join(
        f"<tr><td class=n>{d['id']}</td><td><span class=sw data-g=\"{groups.index(d['group'])}\"></span>{esc(d['label'])}</td>"
        f"<td>{esc(d['group'])}</td><td class=n>{d['cells']:,}</td><td class=n>{d['ctr']:.1f}</td>"
        f"<td class=n>{d['ko']:.1f}</td><td class=n>{(d['ko'] / d['ctr']) if d['ctr'] > 0 else float('inf'):.2f}</td>"
        f"<td class=mk>{esc(d['markers'])}</td></tr>" for d in order)
    pc = ", ".join(str(c) for c in par["clusters"])
    draft_note = ('<span class=badge>Draft labels</span>' if draft else "")
    title = f"{day} {name.capitalize()} Subclusters"

    page = TEMPLATE
    for k, v in {
        "@@TITLE@@": esc(title), "@@DAY@@": day, "@@NAME@@": esc(name), "@@RES@@": str(res), "@@PRES@@": str(par["resolution"]),
        "@@PCL@@": pc, "@@DRAFT@@": draft_note, "@@TRS@@": trs,
        "@@CELLS@@": kpi["cells"], "@@NCTR@@": kpi["ctr"], "@@NKO@@": kpi["ko"],
        "@@HCTR@@": str(kpi["hormone_ctr"]), "@@HKO@@": str(kpi["hormone_ko"]),
        "@@KCTR@@": str(kpi["koonly_ctr"]), "@@KKO@@": str(kpi["koonly_ko"]),
        "@@DOT@@": b64(out / f"dotplot_genes_r{res}.png"),
        "@@DATA@@": json.dumps(data, separators=(",", ":")),
    }.items():
        page = page.replace(k, v)
    dest = out / f"{day}_{name}_r{res}_report.html"
    dest.write_text(page)
    print(f"wrote {dest} ({dest.stat().st_size / 1e6:.1f} MB)")


TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>@@TITLE@@</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js"></script>
<style>
:root{color-scheme:light;--bg:#fcfcfb;--card:#ffffff;--fg:#0b0b0b;--fg2:#52514e;--mut:#7a7975;--line:#e4e3de;
--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--s4:#eda100;--s5:#e87ba4;--s6:#008300;--s7:#6250d6;--s8:#e34948;
--ctr:#2a78d6;--ko:#eb6834;--ctr2:#86b6ef;--ko2:#f3a682;--seq0:#e9eef5;--seq1:#184f95}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--bg:#1a1a19;--card:#232322;--fg:#f0efec;
--fg2:#c3c2b7;--mut:#94938c;--line:#383835;--s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--s5:#d55181;--s6:#0ca30c;
--s7:#9085e9;--s8:#e66767;--ctr:#3987e5;--ko:#d95926;--ctr2:#1c5cab;--ko2:#9c3f1b;--seq0:#2c2c2a;--seq1:#86b6ef}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#1a1a19;--card:#232322;--fg:#f0efec;--fg2:#c3c2b7;--mut:#94938c;--line:#383835;
--s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--s5:#d55181;--s6:#0ca30c;--s7:#9085e9;--s8:#e66767;--ctr:#3987e5;
--ko:#d95926;--ctr2:#1c5cab;--ko2:#9c3f1b;--seq0:#2c2c2a;--seq1:#86b6ef}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1180px;margin:0 auto;padding:28px 16px 64px}
h1{font-size:26px;margin:0 0 6px;letter-spacing:-.01em} h2{font-size:18px;margin:40px 0 6px}
p{margin:6px 0} .sub{color:var(--fg2)} .mut{color:var(--mut);font-size:13px}
.badge{display:inline-block;font-size:12px;padding:2px 8px;border-radius:999px;border:1px solid var(--line);color:var(--fg2);margin-left:8px;vertical-align:middle}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px;margin:20px 0 8px}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.kpi .v{font-size:26px;font-weight:600;font-variant-numeric:tabular-nums} .kpi .l{color:var(--fg2);font-size:13px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px}
.ctl{display:flex;flex-wrap:wrap;gap:10px 16px;align-items:center;margin:4px 0 10px;font-size:14px}
.ctl label{color:var(--fg2)} select{font:inherit;padding:4px 8px;border-radius:6px;border:1px solid var(--line);background:var(--card);color:var(--fg)}
.seg{display:inline-flex;border:1px solid var(--line);border-radius:7px;overflow:hidden}
.seg button{font:inherit;font-size:13px;padding:4px 10px;border:0;background:transparent;color:var(--fg2);cursor:pointer}
.seg button[aria-pressed="true"]{background:var(--fg);color:var(--bg)}
#u3{width:100%;height:640px}
#comp{width:100%;height:440px}
#u2{width:100%;height:480px}
.tw{overflow-x:auto} table{border-collapse:collapse;font-size:13px;width:100%}
th,td{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{color:var(--fg2);font-weight:600} td.n,th.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
td.mk{color:var(--fg2);font-size:12px}
.sw{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:7px;vertical-align:baseline}
img{max-width:100%;height:auto;background:#fff;border-radius:6px;display:block}
@media (max-width:640px){#u3{height:440px}}
</style></head><body><main>
<h1>@@DAY@@ endocrine subclusters @@DRAFT@@</h1>
<p class=sub>@@CELLS@@ cells (@@NCTR@@ control, @@NKO@@ ISL1 KO) from the @@DAY@@ overview clusters @@PCL@@ at resolution
@@PRES@@, re-clustered on their own. Subclusters shown at resolution @@RES@@.</p>
<div class=kpis>
 <div class=kpi><div class=v>@@HCTR@@% → @@HKO@@%</div><div class=l>hormone+ (β-, α-, δ-like) share of the compartment, control → KO</div></div>
 <div class=kpi><div class=v>@@KCTR@@% → @@KKO@@%</div><div class=l>KO-only states (2 subclusters), control → KO</div></div>
 <div class=kpi><div class=v>@@CELLS@@</div><div class=l>endocrine-compartment cells, 4 samples (2 control, 2 KO)</div></div>
</div>

<h2>3D UMAP</h2>
<p class=mut>Drag to rotate, scroll to zoom, hover a cell for its subcluster. Labels sit at each subcluster's centre.</p>
<div class=card>
<div class=ctl>
 <label for=col>Colour by</label>
 <select id=col></select>
 <span class=seg role=group aria-label="Cells shown">
  <button data-f=all aria-pressed=true>All cells</button><button data-f=0 aria-pressed=false>Control</button><button data-f=1 aria-pressed=false>KO</button>
 </span>
 <label><input type=checkbox id=lbl checked> subcluster labels</label>
</div>
<div id=u3></div>
</div>

<h2>Composition: control vs KO</h2>
<p class=mut>Share of each genotype's compartment cells in each subcluster (each genotype sums to 100%). Grouped by label group.</p>
<div class=card><div id=comp></div></div>

<h2>Subclusters</h2>
<div class="card tw"><table>
<thead><tr><th class=n>Sub</th><th>Label</th><th>Group</th><th class=n>Cells</th><th class=n>% ctr</th><th class=n>% KO</th><th class=n>KO / ctr</th><th>Top markers (vs rest of compartment)</th></tr></thead>
<tbody>@@TRS@@</tbody></table></div>

<h2>2D UMAP, split by genotype</h2>
<p class=mut>Same graph as the 3D UMAP, in two dimensions. Colours follow &ldquo;Colour by&rdquo; above; in the control and KO
panels the other genotype is drawn in grey underneath. Hover a cell for its subcluster, drag to zoom, double-click to reset.</p>
<div class=card><div id=u2></div></div>

<h2>Endocrine genes by subcluster</h2>
<div class=card><img src="data:image/png;base64,@@DOT@@" alt="Dot plot of endocrine genes across subclusters"></div>

<h2>Methods</h2>
<p class=sub>Per-day merge of the 4 samples, no batch correction (step 7.1). Overview graph from SCT B (SCTransform v2,
S−G2M regressed, 3,000 variable genes without MT- genes) and Leiden (7.1b). The compartment's cells were taken out and SCT B
refit on them alone, then 30 PCs, 15 neighbours, UMAP in 2D and 3D from the same graph, Leiden (igraph, 2 iterations,
seed 1234) (7.1d). Markers: Wilcoxon, each subcluster vs the rest of the compartment, on unregressed log1p values.
Gene colours on the 3D UMAP are log1p values capped at each gene's 99th percentile. Labels are from marker genes and
are a draft until confirmed.</p>
</main>
<script>
const D = @@DATA@@;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const N = D.x.length;
const gcol = i => css('--s' + (i + 1));
const subGroup = {}; D.clusters.forEach(c => subGroup[c.id] = D.groups.indexOf(c.group));
const subLabel = {}; D.clusters.forEach(c => subLabel[c.id] = c.label);
let filt = 'all', mode = 'group';

const sel = document.getElementById('col');
[['group', 'Label group'], ['sub', 'Subcluster'], ['geno', 'Genotype'], ['sample', 'Sample'], ['parent', 'Overview cluster']]
  .forEach(([v, t]) => sel.add(new Option(t, v)));
const og = document.createElement('optgroup'); og.label = 'Gene (log1p)';
D.genes.forEach(g => og.appendChild(new Option(g, 'gene:' + g))); sel.appendChild(og);

function idx() { const r = []; for (let i = 0; i < N; i++) if (filt === 'all' || D.geno[i] == filt) r.push(i); return r; }
function pick(a, ix) { return ix.map(i => a[i]); }
function pct(arr, p) { const s = arr.filter(v => v > 0).sort((a, b) => a - b); return s.length ? s[Math.floor(p * (s.length - 1))] : 1; }

function catTraces(ix, keyFn, names, colours) {
  const by = {}; ix.forEach(i => { const k = keyFn(i); (by[k] = by[k] || []).push(i); });
  return Object.keys(names).filter(k => by[k]).map(k => ({
    type: 'scatter3d', mode: 'markers', name: names[k] + ' (' + by[k].length + ')',
    x: pick(D.x, by[k]), y: pick(D.y, by[k]), z: pick(D.z, by[k]),
    text: by[k].map(i => 'subcluster ' + D.sub[i] + ': ' + subLabel[D.sub[i]] + '<br>' + D.samples[D.sample[i]]),
    hovertemplate: '%{text}<extra></extra>', marker: { size: 2, color: colours[k], opacity: .85 }
  }));
}

function labelTrace(ix) {
  const acc = {}; ix.forEach(i => { const s = D.sub[i]; (acc[s] = acc[s] || [[], [], []]); acc[s][0].push(D.x[i]); acc[s][1].push(D.y[i]); acc[s][2].push(D.z[i]); });
  const med = a => { const s = [...a].sort((p, q) => p - q); return s[Math.floor(s.length / 2)]; };
  const ks = Object.keys(acc);
  return { type: 'scatter3d', mode: 'text', showlegend: false, hoverinfo: 'skip',
    x: ks.map(k => med(acc[k][0])), y: ks.map(k => med(acc[k][1])), z: ks.map(k => med(acc[k][2])),
    text: ks.map(k => '<b>' + k + '</b>'), textfont: { size: 14, color: css('--fg') } };
}

function draw() {
  const ix = idx(); let tr;
  if (mode === 'group') {
    const nm = {}, co = {}; D.groups.forEach((g, i) => { nm[i] = g; co[i] = gcol(i); });
    tr = catTraces(ix, i => subGroup[D.sub[i]], nm, co);
  } else if (mode === 'sub') {
    const nm = {}, co = {}; D.clusters.forEach(c => { nm[c.id] = c.id + ' ' + c.label; co[c.id] = gcol(subGroup[c.id]); });
    tr = catTraces(ix, i => D.sub[i], nm, co);
  } else if (mode === 'geno') {
    tr = catTraces(ix, i => D.geno[i], { 0: 'control', 1: 'ISL1 KO' }, { 0: css('--ctr'), 1: css('--ko') });
  } else if (mode === 'sample') {
    const co = { 0: css('--ctr'), 1: css('--ctr2'), 2: css('--ko'), 3: css('--ko2') };
    tr = catTraces(ix, i => D.sample[i], Object.fromEntries(D.samples.map((s, i) => [i, s])), co);
  } else if (mode === 'parent') {
    const nm = {}, co = {}; D.parents.forEach((p, i) => { nm[p] = 'overview cluster ' + p; co[p] = gcol(i); });
    tr = catTraces(ix, i => D.parent[i], nm, co);
  } else {
    const g = mode.slice(5), e = D.expr[g], vmax = pct(e, .99) || 1;
    // draw low values first so expressing cells sit on top
    const ord = [...ix].sort((a, b) => e[a] - e[b]);
    tr = [{ type: 'scatter3d', mode: 'markers', name: g, showlegend: false,
      x: pick(D.x, ord), y: pick(D.y, ord), z: pick(D.z, ord),
      text: ord.map(i => g + ' ' + e[i] + '<br>subcluster ' + D.sub[i] + ': ' + subLabel[D.sub[i]]),
      hovertemplate: '%{text}<extra></extra>',
      marker: { size: 2, opacity: .9, color: pick(e, ord), cmin: 0, cmax: vmax,
        colorscale: [[0, css('--seq0')], [1, css('--seq1')]],
        colorbar: { title: { text: g, font: { color: css('--fg2') } }, thickness: 12, len: .6, tickfont: { color: css('--fg2') } } } }];
  }
  if (document.getElementById('lbl').checked) tr.push(labelTrace(ix));
  const ax = t => ({ title: { text: t }, showgrid: true, gridcolor: css('--line'), zeroline: false, showticklabels: false,
    showbackground: false, color: css('--fg2') });
  const old = document.getElementById('u3').layout;
  Plotly.react('u3', tr, {
    paper_bgcolor: 'rgba(0,0,0,0)', font: { family: 'system-ui, sans-serif', color: css('--fg') },
    margin: { l: 0, r: 0, t: 0, b: 0 }, legend: { itemsizing: 'constant', font: { size: 12 }, bgcolor: 'rgba(0,0,0,0)' },
    scene: { xaxis: ax('UMAP 1'), yaxis: ax('UMAP 2'), zaxis: ax('UMAP 3'), camera: old && old.scene ? old.scene.camera : undefined },
    uirevision: 'keep'
  }, { responsive: true, displaylogo: false });
}

function catSpec() {
  const nm = {}, co = {};
  if (mode === 'group') { D.groups.forEach((g, i) => { nm[i] = g; co[i] = gcol(i); }); return { key: i => subGroup[D.sub[i]], nm, co }; }
  if (mode === 'sub') { D.clusters.forEach(c => { nm[c.id] = c.id + ' ' + c.label; co[c.id] = gcol(subGroup[c.id]); }); return { key: i => D.sub[i], nm, co }; }
  if (mode === 'geno') return { key: i => D.geno[i], nm: { 0: 'control', 1: 'ISL1 KO' }, co: { 0: css('--ctr'), 1: css('--ko') } };
  if (mode === 'sample') return { key: i => D.sample[i], nm: Object.fromEntries(D.samples.map((s, i) => [i, s])),
    co: { 0: css('--ctr'), 1: css('--ctr2'), 2: css('--ko'), 3: css('--ko2') } };
  if (mode === 'parent') { D.parents.forEach((p, i) => { nm[p] = 'overview cluster ' + p; co[p] = gcol(i); }); return { key: i => D.parent[i], nm, co }; }
  return null;
}

function draw2d() {
  const el = document.getElementById('u2'), stack = el.clientWidth < 760, MS = 5;
  const panels = [['all', 'All cells'], [0, 'Control (KO in grey)'], [1, 'ISL1 KO (control in grey)']];
  const lo = a => Math.min(...a), hi = a => Math.max(...a), pad = (a, b) => [a - .03 * (b - a), b + .03 * (b - a)];
  const xr = pad(lo(D.u), hi(D.u)), yr = pad(lo(D.v), hi(D.v));
  const spec = catSpec(), tr = [], ann = [], lay = {};
  const hov = i => 'subcluster ' + D.sub[i] + ': ' + subLabel[D.sub[i]] + '<br>' + D.samples[D.sample[i]];
  panels.forEach(([f, t], p) => {
    const n = p === 0 ? '' : String(p + 1), xa = 'x' + n, ya = 'y' + n;
    const fg = [], bg = [];
    for (let i = 0; i < N; i++) (f === 'all' || D.geno[i] == f ? fg : bg).push(i);
    if (bg.length) tr.push({ type: 'scattergl', mode: 'markers', x: pick(D.u, bg), y: pick(D.v, bg), xaxis: xa, yaxis: ya,
      showlegend: false, hoverinfo: 'skip', marker: { size: MS - 1, color: css('--line') } });
    if (spec) {
      const by = {}; fg.forEach(i => { const k = spec.key(i); (by[k] = by[k] || []).push(i); });
      Object.keys(spec.nm).filter(k => by[k]).forEach(k => tr.push({ type: 'scattergl', mode: 'markers',
        name: spec.nm[k], legendgroup: 'g' + k, showlegend: p === 0, x: pick(D.u, by[k]), y: pick(D.v, by[k]),
        xaxis: xa, yaxis: ya, text: by[k].map(hov), hovertemplate: '%{text}<extra></extra>',
        marker: { size: MS, color: spec.co[k], opacity: .85 } }));
    } else {
      const g = mode.slice(5), e = D.expr[g], vmax = pct(e, .99) || 1, ord = [...fg].sort((a, b) => e[a] - e[b]);
      tr.push({ type: 'scattergl', mode: 'markers', showlegend: false, x: pick(D.u, ord), y: pick(D.v, ord),
        xaxis: xa, yaxis: ya, text: ord.map(i => g + ' ' + e[i] + '<br>' + hov(i)), hovertemplate: '%{text}<extra></extra>',
        marker: { size: MS, opacity: .9, color: pick(e, ord), cmin: 0, cmax: vmax, showscale: p === 2,
          colorscale: [[0, css('--seq0')], [1, css('--seq1')]],
          colorbar: { title: { text: g, font: { color: css('--fg2') } }, thickness: 12, len: stack ? .25 : .7, tickfont: { color: css('--fg2') } } } });
    }
    if (document.getElementById('lbl').checked) {
      const acc = {}; fg.forEach(i => { (acc[D.sub[i]] = acc[D.sub[i]] || [[], []]); acc[D.sub[i]][0].push(D.u[i]); acc[D.sub[i]][1].push(D.v[i]); });
      const med = a => { const s = [...a].sort((q, r) => q - r); return s[Math.floor(s.length / 2)]; };
      const ks = Object.keys(acc).filter(k => acc[k][0].length >= 10);
      tr.push({ type: 'scatter', mode: 'text', showlegend: false, hoverinfo: 'skip', xaxis: xa, yaxis: ya,
        x: ks.map(k => med(acc[k][0])), y: ks.map(k => med(acc[k][1])), text: ks.map(k => '<b>' + k + '</b>'),
        textfont: { size: 13, color: css('--fg') } });
    }
    const axis = { range: null, showgrid: false, zeroline: false, showticklabels: false, ticks: '', showline: false };
    lay['xaxis' + n] = { ...axis, range: xr }; lay['yaxis' + n] = { ...axis, range: yr };
    ann.push({ text: t + ' (' + fg.length.toLocaleString() + ')', xref: xa + ' domain', yref: ya + ' domain', x: .5, y: 1.0,
      yanchor: 'bottom', showarrow: false, font: { size: 13, color: css('--fg2') } });
  });
  Plotly.react('u2', tr, Object.assign(lay, {
    grid: { rows: stack ? 3 : 1, columns: stack ? 1 : 3, pattern: 'independent', xgap: .04, ygap: .08 },
    height: stack ? 1150 : 480, annotations: ann, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { family: 'system-ui, sans-serif', color: css('--fg') }, margin: { l: 4, r: 4, t: 24, b: 4 },
    legend: { orientation: 'h', x: 0, y: -.02, yanchor: 'top', itemsizing: 'constant', font: { size: 12 }, bgcolor: 'rgba(0,0,0,0)' },
    uirevision: 'keep2d'
  }), { responsive: true, displaylogo: false });
}

function comp() {
  const c = D.comp, y = c.map(d => d.id + '  ' + d.label);
  const mk = (k, nm, col) => ({ type: 'bar', orientation: 'h', name: nm, y: y, x: c.map(d => d[k]),
    marker: { color: css(col), line: { color: css('--card'), width: 2 } },
    hovertemplate: '%{y}<br>' + nm + ': %{x:.1f}% of compartment cells<extra></extra>' });
  Plotly.react('comp', [mk('ctr', 'control', '--ctr'), mk('ko', 'ISL1 KO', '--ko')], {
    barmode: 'group', bargap: .35, bargroupgap: .08, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    font: { family: 'system-ui, sans-serif', color: css('--fg2'), size: 12 },
    margin: { l: 10, r: 10, t: 10, b: 40 },
    xaxis: { title: { text: '% of genotype\'s compartment cells' }, gridcolor: css('--line'), zeroline: false, ticksuffix: '%' },
    yaxis: { autorange: 'reversed', automargin: true, ticks: '' },
    legend: { orientation: 'h', x: 1, xanchor: 'right', y: 1.06 }
  }, { responsive: true, displaylogo: false });
}

document.querySelectorAll('.sw').forEach(s => s.style.background = gcol(+s.dataset.g));
sel.onchange = () => { mode = sel.value; draw(); draw2d(); };
document.getElementById('lbl').onchange = () => { draw(); draw2d(); };
document.querySelectorAll('.seg button').forEach(b => b.onclick = () => {
  filt = b.dataset.f; document.querySelectorAll('.seg button').forEach(x => x.setAttribute('aria-pressed', x === b)); draw();
});
draw(); draw2d(); comp();
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { draw(); draw2d(); comp(); });
let rz; addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(draw2d, 200); });
</script></body></html>"""


if __name__ == "__main__":
    main()
