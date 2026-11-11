"""Step 1.4: per-sample cell accounting across the QC steps (1.1 -> 1.3 -> 5.2 merge).

Each removed cell is counted once, at the first step that removes it:
  Cell Ranger cells -> 1.1 low counts/genes (4 MAD, either) -> 1.1 high %MT (> cutoff, among cells
  that passed counts/genes) -> 1.2 SoupX (removes ambient UMIs, not cells) -> 1.3 doublets
  (scDblFinder or Scrublet; split into scDblFinder only / Scrublet only / both) -> 5.2 merged.
The 1.1 flags are recomputed with the pipeline's own is_mad_outlier() and qc_RNA.yaml and checked
against the saved passed_qc_step1_2 column, so the split matches what was actually removed.

Writes ${OUT_DIR}/qc_accounting/qc_cell_accounting.tsv and qc_cell_accounting.html.
Run: sbatch analysis/scripts/1.4_qc_cell_accounting.sb
"""
import base64
import importlib.util
import io
import json
import os
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = Path(os.environ["OUT_DIR"])
RES = OUT_DIR / "qc_accounting"


def load_qc_module():
    spec = importlib.util.spec_from_file_location("rna_qc", REPO / "pipeline" / "scripts" / "1.1_rna_qc.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    RES.mkdir(parents=True, exist_ok=True)
    q = load_qc_module()
    cfg_all = yaml.safe_load((REPO / "pipeline" / "configs" / "qc_RNA.yaml").read_text())
    samples = pd.read_csv(os.environ["SAMPLE_SHEET"], sep="\t")["sample_id"].tolist()
    merged = ad.read_h5ad(OUT_DIR / "integration" / "merged_rna_raw.h5ad", backed="r")
    rename = pd.read_csv(REPO / "study_config" / "rename_manifest.csv").set_index("old_name")["new_name"]
    n_merged = merged.obs["sample_id"].astype(str).value_counts()

    rows = []
    for s in samples:
        a = ad.read_h5ad(OUT_DIR / "rna_qc" / "h5ad" / f"{s}.qc_annotated.h5ad", backed="r")
        obs = a.obs
        n_mads = cfg_all["n_mads"]
        low = q.is_mad_outlier(obs["log1p_total_counts"].to_numpy(), n_mads) | \
            q.is_mad_outlier(obs["log1p_n_genes_by_counts"].to_numpy(), n_mads)
        mt_cut = cfg_all["mt_pct_max_by_sample"][s]
        mt = obs["pct_counts_mt"].to_numpy() > mt_cut
        keep = ~(low | mt)
        assert (keep == obs["passed_qc_step1_2"].to_numpy()).all(), f"{s}: recomputed 1.1 flags differ"
        umis_in = float(obs.loc[keep, "total_counts"].sum())

        soup = json.loads((OUT_DIR / "remove_background" / "summary" / f"{s}.soupx_summary.json").read_text())
        d = ad.read_h5ad(OUT_DIR / "detect_doublets" / "h5ad" / f"{s}.doublets_annotated.h5ad", backed="r").obs
        sd, sr = d["scDblFinder_predicted_doublet"].astype(bool), d["scrublet_predicted_doublet"].astype(bool)
        n_after_dbl = int((~(sd | sr)).sum())
        rows.append({
            "sample": rename[s],
            "cellranger_cells": a.n_obs,
            "qc_low_counts_genes": int(low.sum()),
            "qc_high_mt": int((mt & ~low).sum()),
            "after_qc": int(keep.sum()),
            "soupx_cells_removed": soup.get("n_cells_excluded_by_precluster_override", 0),
            "soupx_pct_umis_removed": round(100 * soup["total_counts_removed"] / umis_in, 1),
            "soupx_tfidfMin": soup["tfidf_min_used"],
            "doublet_scDblFinder_only": int((sd & ~sr).sum()),
            "doublet_Scrublet_only": int((~sd & sr).sum()),
            "doublet_both": int((sd & sr).sum()),
            "doublets_removed": int((sd | sr).sum()),
            "after_doublets": n_after_dbl,
            "in_merged_object": int(n_merged.get(rename[s], 0)),
            "mt_cutoff_pct": mt_cut,
        })
    t = pd.DataFrame(rows)
    t["pct_retained"] = (100 * t["in_merged_object"] / t["cellranger_cells"]).round(1)
    tot = t.select_dtypes("number").sum()
    t.to_csv(RES / "qc_cell_accounting.tsv", sep="\t", index=False)
    print(t.to_string(index=False))

    # stacked bars: where each Cell Ranger cell ended up
    parts = [("in_merged_object", "kept", "#4c9a6a"), ("qc_low_counts_genes", "low counts/genes (4 MAD)", "#9aa3ad"),
             ("qc_high_mt", "%MT > 10", "#e3a72f"), ("doublet_scDblFinder_only", "doublet: scDblFinder only", "#c97b84"),
             ("doublet_Scrublet_only", "doublet: Scrublet only", "#8e6bbf"), ("doublet_both", "doublet: both tools", "#b2364a")]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, frac in zip(axs, (False, True)):
        bottom = np.zeros(len(t))
        for col, lab, c in parts:
            v = t[col].to_numpy(float)
            if frac:
                v = 100 * v / t["cellranger_cells"].to_numpy()
            ax.bar(t["sample"], v, bottom=bottom, color=c, label=lab, width=0.75)
            bottom += v
        ax.set_ylabel("% of Cell Ranger cells" if frac else "cells")
        ax.tick_params(axis="x", rotation=60, labelsize=8)
        for lab in ax.get_xticklabels():
            lab.set_ha("right")
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].legend(fontsize=8, frameon=False, loc="upper right")
    axs[1].set_ylim(0, 100)
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=130); plt.close(fig)
    img = base64.b64encode(buf.getvalue()).decode()

    def fmt(v):
        return f"{v:,}" if isinstance(v, (int, np.integer)) else str(v)

    cols = [("sample", "Sample"), ("cellranger_cells", "Cell Ranger cells"), ("qc_low_counts_genes", "Low counts/genes"),
            ("qc_high_mt", "%MT &gt; 10"), ("after_qc", "After QC 1.1"), ("soupx_cells_removed", "SoupX cells"),
            ("soupx_pct_umis_removed", "SoupX % UMIs removed"), ("doublet_scDblFinder_only", "scDblFinder only"),
            ("doublet_Scrublet_only", "Scrublet only"), ("doublet_both", "Both"), ("doublets_removed", "Doublets removed"),
            ("after_doublets", "After doublets"), ("in_merged_object", "In merged object"), ("pct_retained", "% retained")]
    head = "".join(f"<th>{h}</th>" for _, h in cols)
    body = "".join("<tr>" + "".join(f"<td>{fmt(r[c])}</td>" for c, _ in cols) + "</tr>" for r in t.to_dict("records"))
    total = {c: (fmt(int(tot[c])) if c in tot and c not in ("pct_retained", "soupx_pct_umis_removed", "mt_cutoff_pct") else "")
             for c, _ in cols}
    total["sample"] = "Total"
    total["pct_retained"] = f"{100 * tot['in_merged_object'] / tot['cellranger_cells']:.1f}"
    body += "<tr class=tot>" + "".join(f"<td>{total[c]}</td>" for c, _ in cols) + "</tr>"
    tf = ", ".join(f"{r['sample']} ({r['soupx_tfidfMin']})" for r in rows if r["soupx_tfidfMin"] != "default(1)")
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>QC Cell Accounting</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc;--card:#fff}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36;--card:#232321}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1180px;margin:0 auto;padding:24px 16px 60px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:18px;margin:28px 0 8px}}
.sub{{color:var(--mut);margin:0 0 16px}}
.tw{{overflow-x:auto}} table{{border-collapse:collapse;font-size:13px}}
th,td{{padding:4px 8px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap;font-variant-numeric:tabular-nums}}
th{{vertical-align:bottom;white-space:normal;max-width:90px}} th:first-child,td:first-child{{text-align:left}}
tr.tot td{{font-weight:600;border-top:2px solid var(--fg)}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px}}
ul{{padding-left:20px}}
</style></head><body><main>
<h1>QC cell accounting</h1>
<p class=sub>ISL1 scRNA-seq, 12 samples (D15, D22 from Diff37; D40 combined sequencing from Diff34) · out_all/</p>
<div class=tw><table><tr>{head}</tr>{body}</table></div>
<h2>Where each Cell Ranger cell went</h2>
<img src="data:image/png;base64,{img}" alt="stacked bar of cells kept and removed per sample">
<h2>How to read this</h2>
<ul>
<li>Each removed cell is counted once, at the first step that removed it. Cells failing both the
counts/genes and %MT filters are counted under low counts/genes.</li>
<li><b>QC 1.1</b>: low counts or low genes = more than 4 MADs below the sample's own median (log scale);
%MT cutoff fixed at 10% for every sample. No ribosomal filter.</li>
<li><b>SoupX 1.2</b> removes ambient UMIs from each cell, not cells. The percent is UMIs removed out of all
UMIs in the post-QC cells. Fallback tfidfMin used: {tf or "none"}; all other samples used SoupX defaults.</li>
<li><b>Doublets 1.3</b>: a cell is removed if scDblFinder or Scrublet calls it a doublet. scDblFinder prior
1.5% per 1,000 cells (10x OCM); Scrublet's expected rate = scDblFinder's observed rate in that sample.</li>
<li><b>Merge 5.2</b>: cells in the merged object, which should equal "after doublets".</li>
</ul>
</main></body></html>"""
    (RES / "qc_cell_accounting.html").write_text(html)
    print("wrote", RES)


if __name__ == "__main__":
    main()
