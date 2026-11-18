"""Step 7.1: one time point on its own -- merge its 4 samples, process to UMAP. No batch correction.

Decisions (user, 2026-11-17):
  - each day on its own, the 4 samples merged; no batch correction (one differentiation batch and
    OCM channel per day);
  - damaged cells removed here, after doublet removal: nuclear fraction > MAX_NF (0.6). Steps 1.0-1.3
    (NF < 0.10, 4 MAD + 2,000 genes + 10% MT, SoupX, scDblFinder) are upstream and unchanged;
  - log1p normalization + scaling on the merged object (not per sample: per-sample scaling would act
    as a batch correction and remove real ctr vs KO mean differences);
  - 3,000 HVGs, ranked within each sample and combined (batch_key="sample");
  - cell-cycle: regress out the S - G2M difference (keeps cycling vs non-cycling apart, removes the
    S vs G2M split), on the merged object, before scaling. Regression and scaling feed PCA only;
    X stays unregressed log1p for genes, markers and DE;
  - stop at UMAP; clustering is done day by day afterwards.

Inputs: ${OUT_DIR}/detect_doublets/h5ad/{sample}.doublets_removed.h5ad (X = SoupX counts; layers counts,
soupx_counts; obs carries nuclear_fraction from step 1.0).
Processing: concat -> drop NF > MAX_NF -> normalize_total(1e4) -> log1p -> cell-cycle scores (Seurat
cc.genes, Tirosh 2016) -> 3,000 HVGs (flavor seurat, batch_key sample) -> on HVGs: regress_out(S - G2M)
-> scale(max_value 10) -> PCA (50) -> neighbours (30 PCs, 15 neighbours) -> UMAP. Seed 1234.
Writes ${OUT_DIR}/day_analysis/{day}/ (reruns overwrite):
  {day}_merged.h5ad        X = log1p (unregressed), layers counts / soupx_counts, var highly_variable,
                           obsm X_pca + X_umap, obs QC + nuclear_fraction + S_score/G2M_score/phase
  per_sample.csv           cells in, removed (NF > MAX_NF), kept, QC medians, % per phase
  umap_samples.png, umap_genotype_split.png, umap_qc.png, umap_cellcycle.png, umap_genes.png,
  {day}_summary.html
Run: sbatch analysis/scripts/7.1_day_merged.sb D15
"""
import base64
import os
import sys
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

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = Path(os.environ["OUT_DIR"])
DAY = sys.argv[1]
RES_DIR = OUT_DIR / "day_analysis" / DAY
SEED = 1234
MAX_NF = 0.6
N_HVG, N_PCS_STORE, N_PCS, N_NEIGHBORS = 3000, 50, 30, 15
GENES = ["PDX1", "PTF1A", "NKX6-1", "SOX9", "ARX", "GCG", "PAX4", "PAX6", "INS", "SLC18A1", "FEV", "SST",
         "NKX2-2", "NEUROG3", "CHGA", "RFX6", "ISL1", "MKI67"]
SAMPLE_ORDER = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"]
SAMPLE_PAL = {"ctr_rep1": "#2f5d9e", "ctr_rep2": "#7fa7dd", "KO_rep1": "#b23b2a", "KO_rep2": "#eb9a84"}
GENO_PAL = {"ctr": "#3b6fb6", "KO": "#d0573f"}
PHASE_PAL = {"G1": "#c9c9c9", "S": "#3a8fc1", "G2M": "#d1495b"}
GREY = "#e3e3e3"


def cc_genes():
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    importr("Seurat")
    return list(ro.r("Seurat::cc.genes$s.genes")), list(ro.r("Seurat::cc.genes$g2m.genes"))


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def load_day():
    rename = pd.read_csv(REPO / "study_config" / "rename_manifest.csv")
    rows = rename[rename["new_name"].str.contains(f"_{DAY}-")]
    assert len(rows) == 4, rows
    parts = []
    for old, new in zip(rows["old_name"], rows["new_name"]):
        a = ad.read_h5ad(OUT_DIR / "detect_doublets" / "h5ad" / f"{old}.doublets_removed.h5ad")
        tag = new.split("-", 1)[1]                      # ctr_rep1 ...
        a.obs["sample"] = tag
        a.obs["genotype"] = tag.split("_")[0]
        a.obs["replicate"] = tag.split("_")[1]
        a.obs_names = [f"{new}#{b}" for b in a.obs_names]
        a.obs = a.obs[["sample", "genotype", "replicate", "total_counts", "n_genes_by_counts", "pct_counts_mt",
                       "pct_counts_ribo", "scDblFinder_doublet_score", "nuclear_fraction"]]
        a.obsm, a.uns, a.obsp = {}, {}, {}
        parts.append(a)
    var = parts[0].var[["gene_ids", "feature_types"]].copy()
    for p in parts:
        assert (p.var_names == var.index).all()
    m = ad.concat(parts, join="inner", merge=None)
    m.var = var
    m.obs["sample"] = pd.Categorical(m.obs["sample"], SAMPLE_ORDER)
    m.obs["genotype"] = pd.Categorical(m.obs["genotype"], ["ctr", "KO"])
    return m


def clean(ax):
    ax.set_xticks([]); ax.set_yticks([])


def cat_scatter(ax, u, values, title, pal, mask=None):
    """Categorical colours; cells outside `mask` drawn in light grey underneath."""
    vals = pd.Series(np.asarray(values).astype(str))
    if mask is None:
        mask = np.ones(len(vals), bool)
    ax.scatter(u[~mask, 0], u[~mask, 1], c=GREY, s=1, linewidths=0, rasterized=True)
    idx = np.random.default_rng(SEED).permutation(np.where(mask)[0])
    ax.scatter(u[idx, 0], u[idx, 1], c=vals.map(pal).to_numpy()[idx], s=1.5, linewidths=0, rasterized=True)
    for k, c in pal.items():
        if (vals[mask] == k).any():
            ax.scatter([], [], c=c, s=25, label=k)
    if (~mask).any():
        ax.scatter([], [], c=GREY, s=25, label="other")
    ax.legend(fontsize=8, frameon=False, loc="best")
    ax.set_title(title, fontsize=11); clean(ax)


def num_scatter(ax, u, v, title, cmap="viridis", vmin=None, vmax=None):
    v = np.asarray(v, float); o = np.argsort(v)
    s_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=1.5, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True)
    plt.colorbar(s_, ax=ax, fraction=0.04)
    ax.set_title(title, fontsize=11); clean(ax)


def main():
    RES_DIR.mkdir(parents=True, exist_ok=True)
    for p in RES_DIR.rglob("*"):                      # rerun replaces the day's outputs
        if p.is_file():
            p.unlink()
    a = load_day()
    n_in = a.obs["sample"].value_counts()
    damaged = a.obs["nuclear_fraction"].to_numpy() > MAX_NF
    n_rm = a.obs.loc[damaged, "sample"].value_counts()
    a = a[~damaged].copy()
    print(f"{DAY}: {int(damaged.sum())} cells with nuclear fraction > {MAX_NF} removed; {a.n_obs} cells kept")

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    s_genes, g2m_genes = cc_genes()
    s_use, g2m_use = [g for g in s_genes if g in a.var_names], [g for g in g2m_genes if g in a.var_names]
    print(f"cell cycle genes found: S {len(s_use)}/{len(s_genes)}, G2M {len(g2m_use)}/{len(g2m_genes)}")
    sc.tl.score_genes_cell_cycle(a, s_genes=s_use, g2m_genes=g2m_use, random_state=SEED)
    a.obs["phase"] = pd.Categorical(a.obs["phase"], ["G1", "S", "G2M"])
    a.obs["S_minus_G2M"] = a.obs["S_score"] - a.obs["G2M_score"]

    sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat", batch_key="sample")
    hv = a[:, a.var["highly_variable"]].copy()
    print(f"HVGs: {hv.n_vars}")
    sc.pp.regress_out(hv, ["S_minus_G2M"], n_jobs=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=N_PCS_STORE, random_state=SEED)
    a.obsm["X_pca"] = hv.obsm["X_pca"]
    a.uns["pca_variance_ratio"] = hv.uns["pca"]["variance_ratio"]
    sc.pp.neighbors(a, use_rep="X_pca", n_pcs=N_PCS, n_neighbors=N_NEIGHBORS, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    u = a.obsm["X_umap"]

    # ---- per-sample table
    ps = pd.DataFrame({"cells_in": n_in, f"removed_NF>{MAX_NF}": n_rm}).reindex(SAMPLE_ORDER).fillna(0).astype(int)
    ps[f"pct_removed"] = (100 * ps[f"removed_NF>{MAX_NF}"] / ps["cells_in"]).round(1)
    g = a.obs.groupby("sample", observed=True)
    ps["cells_kept"] = g.size()
    ps["median_counts"] = g["total_counts"].median().round(0)
    ps["median_genes"] = g["n_genes_by_counts"].median().round(0)
    ps["median_pct_mt"] = g["pct_counts_mt"].median().round(2)
    ps["median_pct_ribo"] = g["pct_counts_ribo"].median().round(2)
    ps["median_nuclear_fraction"] = g["nuclear_fraction"].median().round(3)
    ph = pd.crosstab(a.obs["sample"], a.obs["phase"], normalize="index").mul(100).round(1)
    for c in ["G1", "S", "G2M"]:
        ps[f"pct_{c}"] = ph[c]
    tot = ps[["cells_in", f"removed_NF>{MAX_NF}", "cells_kept"]].groupby(ps.index.str.split("_").str[0]).sum()
    tot.index = [f"{i} (both reps)" for i in tot.index]
    ps = pd.concat([ps, tot])
    ps.loc[tot.index, "pct_removed"] = (100 * ps.loc[tot.index, f"removed_NF>{MAX_NF}"] / ps.loc[tot.index, "cells_in"]).round(1)
    ps.to_csv(RES_DIR / "per_sample.csv")
    print(ps.to_string())

    title = f"{DAY}: 4 samples merged, no batch correction; S-G2M regressed for PCA ({a.n_obs:,} cells)"
    sample = a.obs["sample"].astype(str).to_numpy(); geno = a.obs["genotype"].astype(str).to_numpy()

    # ---- samples / genotype
    fig, axs = plt.subplots(1, 3, figsize=(21, 6.8))
    cat_scatter(axs[0], u, sample, "sample", SAMPLE_PAL)
    cat_scatter(axs[1], u, geno, "genotype", GENO_PAL)
    cat_scatter(axs[2], u, a.obs["replicate"], "replicate", {"rep1": "#555555", "rep2": "#e0a030"})
    fig.suptitle(title); fig.tight_layout(); fig.savefig(RES_DIR / "umap_samples.png", dpi=110); plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(16, 7.5))
    for ax, gt, lab in zip(axs, ("ctr", "KO"), ("Control", "ISL1 KO")):
        m = geno == gt
        cat_scatter(ax, u, sample, f"{lab} ({int(m.sum()):,} cells); other genotype in grey",
                    {k: v for k, v in SAMPLE_PAL.items() if k.startswith(gt)}, mask=m)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(RES_DIR / "umap_genotype_split.png", dpi=120); plt.close(fig)

    # ---- QC
    fig, axs = plt.subplots(2, 3, figsize=(18, 10.5))
    for ax, (col, lab, vmin, vmax) in zip(axs.flat, [
            ("total_counts", "log10 counts", None, None), ("n_genes_by_counts", "log10 genes", None, None),
            ("pct_counts_mt", "% MT", 0, 10), ("pct_counts_ribo", "% ribo", None, None),
            ("scDblFinder_doublet_score", "scDblFinder score", 0, None), ("nuclear_fraction", "nuclear fraction", 0.1, MAX_NF)]):
        v = a.obs[col].to_numpy(float)
        num_scatter(ax, u, np.log10(v) if lab.startswith("log10") else v, lab, vmin=vmin, vmax=vmax)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(RES_DIR / "umap_qc.png", dpi=110); plt.close(fig)

    # ---- cell cycle (after S-G2M regression: S and G2M cells should no longer separate)
    def gx(gname):
        x = a[:, gname].X
        return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.ravel(x)
    fig, axs = plt.subplots(2, 3, figsize=(18, 10.5))
    cat_scatter(axs[0, 0], u, a.obs["phase"], "phase", PHASE_PAL)
    num_scatter(axs[0, 1], u, a.obs["S_score"], "S score", cmap="magma")
    num_scatter(axs[0, 2], u, a.obs["G2M_score"], "G2M score", cmap="magma")
    for ax, gname in zip(axs[1], ["MKI67", "TOP2A", "PCNA"]):
        num_scatter(ax, u, gx(gname), f"{gname} (log1p)", cmap="magma")
    fig.suptitle(title); fig.tight_layout(); fig.savefig(RES_DIR / "umap_cellcycle.png", dpi=110); plt.close(fig)

    # ---- genes
    genes = [gn for gn in GENES if gn in a.var_names]
    ncol = 6; nrow = int(np.ceil(len(genes) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.8 * nrow))
    for ax, gn in zip(axs.flat, genes):
        num_scatter(ax, u, gx(gn), gn, cmap="magma")
    for ax in list(axs.flat)[len(genes):]:
        ax.axis("off")
    fig.suptitle(f"{DAY}: log1p expression"); fig.tight_layout(); fig.savefig(RES_DIR / "umap_genes.png", dpi=100); plt.close(fig)

    a.obs.drop(columns=["S_minus_G2M"]).to_csv(RES_DIR / "cells_obs.tsv.gz", sep="\t")
    a.write_h5ad(RES_DIR / f"{DAY}_merged.h5ad", compression="gzip")

    # ---- HTML
    def table(df):
        return df.to_html(border=0, classes="t")
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{DAY} Merged UMAP</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1180px;margin:0 auto;padding:24px 16px 60px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:30px 0 8px}}
.sub{{color:var(--mut)}} .tw{{overflow-x:auto}}
table.t{{border-collapse:collapse;font-size:13px}} .t th,.t td{{padding:3px 8px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>{DAY}: four samples merged, no batch correction</h1>
<p class=sub>QC upstream: nuclear fraction &ge; 0.10, 4 MAD + at least 2,000 genes + %MT &le; 10, SoupX, scDblFinder.
Here: damaged cells (nuclear fraction &gt; {MAX_NF}) removed; normalize_total (10,000) &rarr; log1p &rarr;
{N_HVG:,} HVGs (ranked within each sample) &rarr; S&minus;G2M regressed &rarr; scaled (merged object) &rarr; PCA &rarr;
{N_PCS} PCs, {N_NEIGHBORS} neighbours &rarr; UMAP. Seed {SEED}. No clustering yet. Gene plots use unregressed log1p values.</p>
<h2>Cells per sample</h2><div class=tw>{table(ps)}</div>
<h2>Samples and genotype</h2><img src="data:image/png;base64,{b64(RES_DIR / 'umap_samples.png')}" alt="UMAP by sample, genotype and replicate">
<img src="data:image/png;base64,{b64(RES_DIR / 'umap_genotype_split.png')}" alt="UMAP split by genotype, other genotype in grey">
<h2>QC metrics</h2><img src="data:image/png;base64,{b64(RES_DIR / 'umap_qc.png')}" alt="UMAP coloured by QC metrics">
<h2>Cell cycle</h2><p class=sub>Phase from scores on log1p data. After regressing S&minus;G2M, S and G2M cells
should no longer separate; cycling and non-cycling cells still can.</p>
<img src="data:image/png;base64,{b64(RES_DIR / 'umap_cellcycle.png')}" alt="UMAP coloured by cell-cycle phase and scores">
<h2>Genes</h2><img src="data:image/png;base64,{b64(RES_DIR / 'umap_genes.png')}" alt="UMAP coloured by marker genes">
</main></body></html>"""
    (RES_DIR / f"{DAY}_summary.html").write_text(html)
    print("wrote", RES_DIR)


if __name__ == "__main__":
    main()
