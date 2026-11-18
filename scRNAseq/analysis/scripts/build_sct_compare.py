"""Build the log1p vs SCTransform report (step 7.1c) on the Mac from outputs fetched from Alpine.

Input:  <fetched>/{D15,D22,D40}/  = copies of ${OUT_DIR}/day_analysis/{day}/sct/
Output: <out.html>
Usage:  ~/anaconda3/bin/python3 build_sct_compare.py <fetched_dir> <out.html>
"""
import base64
import io
import sys
from pathlib import Path

import pandas as pd
from PIL import Image

SRC, OUT = Path(sys.argv[1]), Path(sys.argv[2])
DAYS = [d for d in ("D15", "D22", "D40") if (SRC / d / "method_overlap.csv").exists()]
M = ["log1p", "SCT", "SCT_A", "SCT_B"]
LABEL = {"log1p": "log1p", "SCT": "SCT (S−G2M)", "SCT_A": "SCT A (S−G2M + %MT)", "SCT_B": "SCT B (S−G2M, no MT− HVGs)"}


def img(p, alt):
    """PNG -> JPEG (quality 80) to keep the page small."""
    buf = io.BytesIO()
    Image.open(p).convert("RGB").save(buf, "JPEG", quality=80, optimize=True)
    return f'<img src="data:image/jpeg;base64,{base64.b64encode(buf.getvalue()).decode()}" alt="{alt}">'


def table(df, index=False):
    return f'<div class=tw>{df.to_html(border=0, classes="t", index=index)}</div>'


summary, sections = [], []
for d in DAYS:
    s = SRC / d
    ov = pd.read_csv(s / "method_overlap.csv")
    pcc = pd.read_csv(s / "pc_covariates.csv")
    mix = pd.read_csv(s / "knn_mixing.csv")
    pcs = [f"PC{k}" for k in range(1, 11)]
    pcc["max |rho| PC1-10"] = pcc[pcs].abs().max(axis=1)
    pcc["PC of max"] = pcc[pcs].abs().idxmax(axis=1)
    r2 = pcc.pivot(index="covariate", columns="method", values="R2_30PCs")[M]
    rho = pcc.pivot(index="covariate", columns="method", values="max |rho| PC1-10")[M]
    cov_t = pd.concat({"R² from 30 PCs": r2, "max |Spearman ρ|, PC1–10": rho}, axis=1).loc[r2.index]
    mx = mix.pivot(index="sample", columns="method", values="replicate_mixing_obs_over_exp")[M]
    gs = mix.pivot(index="sample", columns="method", values="genotype_separation_obs_over_exp")[M]
    mix_t = pd.concat({"replicate mixing (obs/exp)": mx, "genotype separation (obs/exp)": gs}, axis=1)
    mix_t = mix_t.reindex(["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"])
    var = pd.read_csv(s / "pca_variance.csv").head(10).set_index("PC").T.mul(100).round(1)
    hvo = pd.read_csv(s / "hvg_only.csv").fillna("")
    mk = pd.read_csv(s / "marker_genes_in_hvg.csv")
    mk = mk.set_index("gene"); mk.columns = [c.replace("_HVG", "") for c in mk.columns]
    mk_diff = mk[mk.nunique(axis=1) > 1].replace({True: "yes", False: "–"})
    ov_t = ov.set_index("method")[["hvg_shared_with_log1p", "hvg_jaccard_log1p", "hvg_shared_with_SCT",
                                   "knn_jaccard_vs_log1p", "knn_jaccard_vs_SCT"]]
    ov_t.index = [LABEL[i] for i in ov_t.index]
    for t in (cov_t, mix_t):
        t.columns = pd.MultiIndex.from_tuples([(a, LABEL[b]) for a, b in t.columns])
    n_cells = int(ov.cells.iloc[0])
    for m in M:
        summary.append({"day": d, "method": LABEL[m], "cells": n_cells,
                        "R² log10 counts": r2.loc["log10 counts", m], "R² % MT": r2.loc["% MT", m],
                        "R² nuclear fraction": r2.loc["nuclear fraction", m], "R² genotype": r2.loc["genotype (KO=1)", m],
                        "replicate mixing": round(mx[m].mean(), 2), "genotype separation": round(gs[m].mean(), 2)})
    sections.append(f"""
<h2 id="{d}">{d} <span class=sub>({n_cells:,} cells)</span></h2>
<p class=sub>Variable genes shared (of 3,000) and mean per-cell overlap of the 15-neighbour sets (Jaccard; 1 = identical graphs).</p>
{table(ov_t, index=True)}
<h3>Samples and genotype</h3><p class=sub>Rows top to bottom: {", ".join(LABEL[m] for m in M)}.</p>
{img(s / 'cmp_samples.png', f'{d} UMAPs by sample and genotype, one row per method')}
<h3>Neighbour mixing</h3>
<p class=sub>Replicate mixing: of a cell's neighbours with its own genotype, the share from the other replicate,
divided by the share expected if replicates mixed at random (1 = fully mixed, lower = replicates apart).
Genotype separation: share of neighbours with the same genotype divided by the share expected at random
(1 = genotypes mixed, higher = genotypes apart).</p>{table(mix_t, index=True)}
<h3>QC metrics</h3>{img(s / 'cmp_qc.png', f'{d} UMAPs coloured by QC metrics')}
<p class=sub>How much each covariate is carried by the embedding: R² of a linear fit of the covariate on the 30 PCs
used for the graph, and the strongest Spearman correlation with any of PC1–10. Lower for technical covariates
(counts, genes, %MT, nuclear fraction) means less technical structure in the embedding; counts and genes
also differ between real cell types, so they will not be zero.</p>{table(cov_t, index=True)}
<p class=sub>Variance explained by PC1–10 (% of the first 50 PCs).</p>{table(var, index=True)}
<h3>Cell cycle</h3>{img(s / 'cmp_cellcycle.png', f'{d} UMAPs coloured by cell cycle')}
<h3>Marker genes (log1p values on both embeddings)</h3>
{"".join(f"<h4>{LABEL[m]} UMAP</h4>" + img(s / f"genes_{m}.png", f"{d} {LABEL[m]} UMAP coloured by marker genes") for m in M)}
<h3>Variable genes in one method only</h3>
<p class=sub>Base SCT vs log1p, top 40 by each method's own ranking (SCT: residual variance; log1p: normalised dispersion).</p>
{table(hvo)}
<p class=sub>Marker genes that are variable in some methods but not others:</p>{table(mk_diff, index=True) if len(mk_diff) else "<p>none</p>"}""")

summ = pd.DataFrame(summary).set_index(["day", "method"])
html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>log1p vs SCTransform</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1240px;margin:0 auto;padding:24px 16px 60px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:21px;margin:40px 0 8px;border-top:1px solid var(--line);padding-top:18px}}
h3{{font-size:17px;margin:24px 0 6px}} h4{{font-size:14px;margin:12px 0 2px;color:var(--mut)}}
.sub{{color:var(--mut);font-weight:normal}} .tw{{overflow-x:auto}} nav a{{margin-right:14px}}
table.t{{border-collapse:collapse;font-size:13px;margin:6px 0}}
.t th,.t td{{padding:3px 9px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>log1p vs SCTransform variants, per day</h1>
<p class=sub>Each day: its 4 samples merged, no batch correction, damaged cells (nuclear fraction &gt; 0.6) removed.
Both routes use the same cells and counts (SoupX), regress the same S&minus;G2M score, keep 3,000 variable genes,
and build the graph from 30 PCs and 15 neighbours, then UMAP (seed 1234). Only the normalization differs.
<b>log1p</b>: normalize_total (10,000) &rarr; log1p &rarr; HVGs ranked within each sample &rarr; S&minus;G2M regressed &rarr; scaled.
<b>SCTransform v2</b> (Seurat 5.3, glmGamPoi), one model on the merged object, Pearson residuals of 3,000 variable features, three variants:
<b>SCT</b> S&minus;G2M regressed; <b>SCT A</b> S&minus;G2M and % MT regressed; <b>SCT B</b> S&minus;G2M regressed and mitochondrial (MT-) genes
removed from the variable features (replaced by the next-ranked genes).
No clustering. Gene panels show the same log1p values on both UMAPs.</p>
<nav>{" ".join(f'<a href="#{d}">{d}</a>' for d in DAYS)}</nav>
<h2>Summary</h2>
<p class=sub>R²: share of the covariate explained by the 30 PCs (lower = less technical signal for counts, %MT and
nuclear fraction). Replicate mixing: mean over the 4 samples, 1 = replicates fully mixed. Genotype separation:
mean over the 4 samples, 1 = genotypes mixed at random.</p>
{table(summ, index=True)}
{''.join(sections)}
</main></body></html>"""
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(html)
print("wrote", OUT, f"{OUT.stat().st_size / 1e6:.1f} MB")
