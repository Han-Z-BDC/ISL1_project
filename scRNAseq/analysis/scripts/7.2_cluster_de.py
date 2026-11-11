"""Step 7.2: KO vs ctr differential expression inside one cluster of one day.

Request (user, 2026-11-09): D15, Harmony (per sample) clustering at r=0.8 (7.1 --cluster), cluster 2.

Design: each sample is the replicate (2 ctr vs 2 KO; rep1/rep2 not paired), so the main test is a
pseudobulk: raw Cell Ranger counts (layer "counts") summed over the cluster's cells per sample ->
edgeR quasi-likelihood (filterByExpr, TMM, ~ genotype, robust QL fit) via 7.2_pseudobulk_edger.R.
A cell-level Wilcoxon (scanpy, log1p of SoupX counts normalised to 1e4) is reported alongside for
the % of cells expressing each gene; its p-values treat cells as independent and are not used for calls.

Writes ${OUT_DIR}/day_analysis/{day}/de_{method}_r{res}_c{cluster}/:
  pb_counts.csv, pb_meta.csv, edger_KO_vs_ctr.csv, edger_logCPM.csv, wilcoxon_KO_vs_ctr.csv,
  de_table.csv (merged), volcano.png, heatmap_top.png, umap_cluster.png, de_report.html
Run: sbatch analysis/scripts/7.2_cluster_de.sb D15 harmony 0.8 2
"""
import argparse
import base64
import os
import subprocess
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = Path(os.environ["OUT_DIR"])
NAMES = {"none": "no integration", "harmony": "Harmony (per sample)", "scvi": "scVI (per sample)"}
SAMPLE_ORDER = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
GENO_PAL = {"ctr": "#3b6fb6", "KO": "#d0573f"}
MARKERS = ["PDX1", "PTF1A", "NKX6-1", "ARX", "GCG", "PAX4", "PAX6", "INS", "SLC18A1", "FEV", "SST",
           "NKX2-2", "RFX6", "ISL1"]
FDR = 0.05


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("day"); ap.add_argument("method"); ap.add_argument("resolution", type=float)
    ap.add_argument("cluster", type=int)
    args = ap.parse_args()
    day, meth, res, clu = args.day, args.method, args.resolution, args.cluster
    dd = OUT_DIR / "day_analysis" / day
    out = dd / f"de_{meth}_r{res}_c{clu}"
    out.mkdir(parents=True, exist_ok=True)

    d = sc.read_h5ad(OUT_DIR / "integration_per_day" / f"{day}.h5ad")
    col = f"{meth}_r{str(res).replace('.', 'p')}"
    lab = pd.read_csv(dd / f"clusters_{day}.csv", index_col=0)[col].reindex(d.obs_names).to_numpy()
    in_c = lab == clu
    cells = d.obs_names[in_c]
    obs = d.obs.loc[cells].copy()
    print(f"{day} {meth} r{res} cluster {clu}: {len(cells)} cells", obs["sample_short"].value_counts().to_dict(), flush=True)

    raw = ad.read_h5ad(OUT_DIR / "integration" / "merged_rna_raw.h5ad", backed="r")
    idx = np.sort(raw.obs_names.get_indexer(cells))
    sub = raw[idx].to_memory()[cells]
    del raw

    # ---- pseudobulk + edgeR ----
    C = sub.layers["counts"]
    C = C.tocsr() if sp.issparse(C) else sp.csr_matrix(C)
    samp = obs["sample_short"].astype(str).to_numpy()
    pb = pd.DataFrame({s: np.asarray(C[samp == s].sum(0)).ravel() for s in SAMPLE_ORDER}, index=sub.var_names)
    pb = pb.astype(int)
    pb.to_csv(out / "pb_counts.csv")
    n_cells = pd.Series(samp).value_counts().reindex(SAMPLE_ORDER)
    meta = pd.DataFrame({"genotype": [s.split("_")[0] for s in SAMPLE_ORDER], "n_cells": n_cells.to_numpy(),
                         "lib_size": pb.sum(0).to_numpy()}, index=SAMPLE_ORDER)
    meta.to_csv(out / "pb_meta.csv")
    r = subprocess.run([os.environ["TOOLS"] + "/Rscript", str(REPO / "analysis/scripts/7.2_pseudobulk_edger.R"), str(out)],
                       capture_output=True, text=True)
    print(r.stdout, r.stderr[-3000:], flush=True)
    r.check_returncode()
    rlog = r.stdout
    eg = pd.read_csv(out / "edger_KO_vs_ctr.csv").set_index("gene")
    lcpm = pd.read_csv(out / "edger_logCPM.csv", index_col=0)[SAMPLE_ORDER]

    # ---- cell-level: % expressing + Wilcoxon (descriptive) ----
    a = ad.AnnData(X=sub.layers["soupx_counts"].copy(), obs=obs[["genotype", "sample_short"]].copy(),
                   var=pd.DataFrame(index=sub.var_names))
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, "genotype", groups=["KO"], reference="ctr", method="wilcoxon", pts=True, use_raw=False)
    w = sc.get.rank_genes_groups_df(a, group="KO").set_index("names")
    w = w.rename(columns={"logfoldchanges": "wilcox_log2FC", "pvals_adj": "wilcox_padj"})
    X = a.X.tocsr() if sp.issparse(a.X) else a.X
    g_ = obs["genotype"].astype(str).to_numpy()
    for g in ("ctr", "KO"):     # fraction of cells with any count (pct_nz_reference is not returned here)
        w[f"pct_{g}"] = pd.Series(np.asarray((X[g_ == g] > 0).mean(0)).ravel(), index=a.var_names)
    w[["wilcox_log2FC", "wilcox_padj", "pct_ctr", "pct_KO"]].to_csv(out / "wilcoxon_KO_vs_ctr.csv")
    pct_s = pd.DataFrame({s: np.asarray((X[samp == s] > 0).mean(0)).ravel() * 100 for s in SAMPLE_ORDER},
                         index=a.var_names).round(1)

    tab = eg[["logFC", "logCPM", "F", "PValue", "FDR"]].copy()
    tab = tab.join(lcpm.add_prefix("logCPM_")).join((w[["pct_ctr", "pct_KO"]] * 100).round(1))
    tab = tab.join(pct_s.add_prefix("pct_")).join(w[["wilcox_log2FC", "wilcox_padj"]])
    # same direction in both replicates of each genotype (min KO > max ctr, or the reverse)
    ko_min, ko_max = tab[["logCPM_KO_rep1", "logCPM_KO_rep2"]].min(1), tab[["logCPM_KO_rep1", "logCPM_KO_rep2"]].max(1)
    c_min, c_max = tab[["logCPM_ctr_rep1", "logCPM_ctr_rep2"]].min(1), tab[["logCPM_ctr_rep1", "logCPM_ctr_rep2"]].max(1)
    tab["reps_separate"] = ((ko_min > c_max) | (ko_max < c_min))
    tab = tab.sort_values("PValue")
    tab.to_csv(out / "de_table.csv")
    sig = tab[tab.FDR < FDR]
    up, down = sig[sig.logFC > 0], sig[sig.logFC < 0]
    print(f"FDR<{FDR}: {len(sig)} (up {len(up)}, down {len(down)})", flush=True)

    # ---- figures ----
    u = np.asarray(d.obsm[f"X_umap_{meth}"]); geno = d.obs["genotype"].astype(str).to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True, sharey=True)
    for ax, g in zip(axes, ("ctr", "KO")):
        ax.scatter(u[:, 0], u[:, 1], s=1, c="#e6e6e6", linewidths=0)
        m = in_c & (geno == g)
        ax.scatter(u[m, 0], u[m, 1], s=4, c=GENO_PAL[g], linewidths=0)
        ax.set_title(f"cluster {clu}: {g} cells (n={m.sum():,})", fontsize=12)
        ax.set_xticks([]); ax.set_yticks([]); [s_.set_visible(False) for s_ in ax.spines.values()]
    fig.suptitle(f"{day}, {NAMES[meth]} UMAP, r={res}", fontsize=13)
    fig.tight_layout(); fig.savefig(out / "umap_cluster.png", dpi=90); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 7.5))
    y = -np.log10(tab.PValue.clip(lower=1e-300))
    col_ = np.where(tab.FDR < FDR, np.where(tab.logFC > 0, GENO_PAL["KO"], GENO_PAL["ctr"]), "#bdbdbd")
    ax.scatter(tab.logFC, y, s=7, c=col_, linewidths=0)
    lab_g = list(up.head(15).index) + list(down.head(15).index) + [g for g in MARKERS if g in tab.index and tab.loc[g, "FDR"] < FDR]
    for g in dict.fromkeys(lab_g):
        ax.annotate(g, (tab.loc[g, "logFC"], -np.log10(max(tab.loc[g, "PValue"], 1e-300))), fontsize=8,
                    xytext=(3, 2), textcoords="offset points",
                    weight="bold" if g in MARKERS else "normal")
    ax.axvline(0, c="#888", lw=0.6)
    ax.set_xlabel("log2 fold change, KO vs ctr (edgeR, pseudobulk)"); ax.set_ylabel("-log10 p")
    ax.set_title(f"cluster {clu}: {len(up)} up / {len(down)} down in KO at FDR < {FDR}  (bold = marker genes)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(out / "volcano.png", dpi=110); plt.close(fig)

    top = list(up.sort_values("PValue").head(30).index) + list(down.sort_values("PValue").head(30).index)
    if top:
        z = lcpm.loc[top]
        z = z.sub(z.mean(1), axis=0).div(z.std(1).replace(0, 1), axis=0)
        fig, ax = plt.subplots(figsize=(4.8, 0.2 * len(top) + 1.5))
        im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-1.6, vmax=1.6)
        ax.set_xticks(range(4)); ax.set_xticklabels(SAMPLE_ORDER, rotation=45, ha="right")
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top, fontsize=7)
        ax.axhline(len(up.head(30)) - 0.5, c="black", lw=0.8)
        fig.colorbar(im, ax=ax, fraction=0.08, label="logCPM z-score")
        ax.set_title("top up (above line) / down in KO", fontsize=10)
        fig.tight_layout(); fig.savefig(out / "heatmap_top.png", dpi=110); plt.close(fig)

    # ---- html ----
    show = ["logFC", "FDR", "logCPM_ctr_rep1", "logCPM_ctr_rep2", "logCPM_KO_rep1", "logCPM_KO_rep2",
            "pct_ctr", "pct_KO", "reps_separate"]

    def fmt(df):
        df = df[show].copy()
        df["logFC"] = df["logFC"].round(2); df["FDR"] = df["FDR"].map(lambda v: f"{v:.2e}")
        for c in show[2:6]:
            df[c] = df[c].round(2)
        return df.to_html(border=0, classes="t")

    mk = tab.loc[[g for g in MARKERS if g in tab.index]]
    mk_missing = [g for g in MARKERS if g not in tab.index]
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{day} Cluster {clu} DE</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1200px;margin:0 auto;padding:24px 16px 60px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:30px 0 8px}}
.sub{{color:var(--mut);margin:0 0 14px}} .tw{{overflow-x:auto}} pre{{font-size:13px;color:var(--mut)}}
table.t{{border-collapse:collapse;font-size:13px}}
table.t th,table.t td{{padding:3px 9px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px}}
.row{{display:flex;flex-wrap:wrap;gap:16px;align-items:flex-start}}
</style></head><body><main>
<h1>{day}, {NAMES[meth]} r={res}, cluster {clu}: KO vs ctr</h1>
<p class=sub>{len(cells):,} cells · out_all/day_analysis/{day}/{out.name}</p>
<div class=tw>{meta.to_html(border=0, classes="t")}</div>
<p class=sub>Pseudobulk: raw Cell Ranger counts summed per sample over this cluster's cells, edgeR quasi-likelihood,
~ genotype, 2 ctr vs 2 KO. With 2 replicates per group the test is conservative, so a short list is expected.
logCPM_* = edgeR normalised log2 CPM per sample. pct_ctr / pct_KO = % of cells with any count (cell level).
reps_separate = both KO samples lie on the same side of both ctr samples.</p>
<pre>{rlog}</pre>
<img src="data:image/png;base64,{b64(out / 'umap_cluster.png')}" alt="cluster on UMAP">
<h2>Volcano</h2>
<div class=row><img style="max-width:62%" src="data:image/png;base64,{b64(out / 'volcano.png')}" alt="volcano">
{f'<img style="max-width:34%" src="data:image/png;base64,{b64(out / "heatmap_top.png")}" alt="heatmap">' if top else ''}</div>
<h2>Up in KO (FDR &lt; {FDR}, top 60 by p)</h2><div class=tw>{fmt(up.head(60))}</div>
<h2>Down in KO (FDR &lt; {FDR}, top 60 by p)</h2><div class=tw>{fmt(down.head(60))}</div>
<h2>Marker genes</h2>
{"<p class=sub>Filtered out as too lowly expressed in this cluster: " + ", ".join(mk_missing) + "</p>" if mk_missing else ""}
<div class=tw>{fmt(mk)}</div>
<p class=sub>Full table (all tested genes, plus Wilcoxon cell-level columns): de_table.csv</p>
</main></body></html>"""
    (out / "de_report.html").write_text(html)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
