"""Step 7.2d: compare the D15 + D22 endocrine integrations (after 7.2b and 7.2c). No clustering (user, 2026-11-19).

Methods (30 dimensions each): uncorrected (SCT B PCA), Harmony, RPCA, scVI. For each: 15 neighbours, UMAP (seed
1234), exactly as the per-day steps. Checks on each embedding's 30-nearest-neighbour sets (exact, Euclidean):
  stage mixing   % of a cell's neighbours from the other stage, by label (all D15 cells would score
                 ~71% and all D22 cells ~29% if stages were fully mixed: the other stage's share of all cells)
  label kept     D22 cells: % of their D22 neighbours with the same D22 label; D15 cells the same with the D15
                 draft labels (high = cell types still separate after correction)
  genotype kept  within each stage: % of a cell's same-stage neighbours with the same genotype (genotype is not
                 corrected, so this should stay near the uncorrected value)
  D15 -> D22     for each D15 cell, the most common D22 label among its D22 neighbours: a readout of what each D15
                 draft state is placed next to (cells with no D22 neighbour: "none")
Gene plots: log1p of SoupX counts normalised to 10,000 per cell (unregressed), as in the per-day reports.
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/ (reruns overwrite these files):
  umap_{method}.csv, metrics_*.csv, d15_to_d22_{method}.csv, umap_overview.png, genes_{method}.png,
  D15_D22_endocrine_integration.html
With INT_SUBSET=<name>: the subset folder integration_D15_D22_endocrine_<name> (see 7.2a).
Run: sbatch analysis/scripts/7.2d_integrate_report.sb [uncorrected,harmony,rpca,scvi]   (default: all four)
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
from sklearn.neighbors import NearestNeighbors

OUT_DIR, REPO_DIR = Path(os.environ["OUT_DIR"]), Path(os.environ["REPO_DIR"])
SUBSET = os.environ.get("INT_SUBSET", "")
OUT = OUT_DIR / "day_analysis" / ("integration_D15_D22_endocrine" + (f"_{SUBSET}" if SUBSET else ""))
METHODS = [("uncorrected", "Uncorrected (SCT B PCA)"), ("harmony", "Harmony"), ("rpca", "RPCA"), ("scvi", "scVI")]
SEED, K = 1234, 30
GREY = "#e6e6e6"
STAGE_COL = {"D15": "#e08214", "D22": "#5e3c99"}
GENO_COL = {"ctr": "#4393c3", "KO": "#d6604d"}
GENES = ["CHGA", "NEUROG3", "NEUROD1", "SOX4", "INS", "NKX6-1", "GCG", "ARX", "IRX2", "SST", "HHEX", "GHRL",
         "FEV", "MNX1", "TPH1", "LMX1A", "PLCL1", "RUNX1T1", "PHOX2A", "ISL1", "PDX1", "MKI67"]
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", ["#d4d4d4", "#fdbb84", "#e34a33", "#7f0000"])


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def style(ax, title):
    ax.set_title(title, fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def cat_panel(ax, u, values, pal, title, legend=True, order=None):
    v = values.astype(str).to_numpy()
    rng = np.random.default_rng(SEED)
    bg = ~np.isin(v, list(pal))
    ax.scatter(u[bg, 0], u[bg, 1], c=GREY, s=1.5, linewidths=0, rasterized=True)
    idx = rng.permutation(np.where(~bg)[0])
    ax.scatter(u[idx, 0], u[idx, 1], c=[pal[x] for x in v[idx]], s=2, linewidths=0, rasterized=True)
    if legend:
        for k in (order or pal):
            if k in pal and (v == k).any():
                ax.scatter([], [], color=pal[k], s=30, label=f"{k} ({(v == k).sum():,})")
        ax.legend(fontsize=7, frameon=False, markerscale=1.2, loc="center left", bbox_to_anchor=(1.0, 0.5))
    style(ax, title)


def main():
    global METHODS
    if len(sys.argv) > 1:   # optional subset, e.g. before scVI has finished; a full run overwrites the report
        keep = sys.argv[1].split(",")
        METHODS = [x for x in METHODS if x[0] in keep]
    a = ad.read_h5ad(OUT / "merged_counts.h5ad")
    obs = a.obs.copy()
    n = len(obs)
    stage = obs.stage.astype(str).to_numpy()
    geno = obs.genotype.astype(str).to_numpy()
    lab = obs.label.astype(str).to_numpy()
    is15 = stage == "D15"
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    genes = [g for g in GENES if g in a.var_names]
    X = a[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)

    cols = pd.read_csv(REPO_DIR / "study_config" / "annotation" / "label_colours.tsv", sep="\t", dtype=str)
    order22 = [l for l in cols.label if l in set(lab[~is15])]
    pal22 = {l: c for l, c in zip(cols.label, cols.colour) if l in order22}
    lab15 = pd.read_csv(REPO_DIR / "study_config" / "annotation" / "D15_endocrine_r0.6.tsv", sep="\t", dtype=str)
    order15 = ["D15 " + l.split(" (")[0] for l in lab15.label]
    pal15 = {l: matplotlib.colors.to_hex(plt.get_cmap("tab20")(i * 2 % 20 + (i * 2 // 20))) for i, l in enumerate(order15)}
    share = {"D15": (~is15).mean(), "D22": is15.mean()}   # other stage's share of all cells

    emb, um, mix_rows, keep_rows, geno_rows, d15 = {}, {}, [], [], [], {}
    for m, name in METHODS:
        e = pd.read_csv(OUT / f"emb_{m}.csv.gz", index_col=0).reindex(obs.index).to_numpy()
        assert not np.isnan(e).any() and e.shape[1] == 30, m
        b = ad.AnnData(np.zeros((n, 1), np.float32), obs=obs[["stage"]].copy())
        b.obsm["E"] = e
        sc.pp.neighbors(b, n_neighbors=15, use_rep="E", random_state=SEED)
        sc.tl.umap(b, random_state=SEED)
        um[m] = b.obsm["X_umap"]
        pd.DataFrame(um[m], index=obs.index, columns=["UMAP1", "UMAP2"]).to_csv(OUT / f"umap_{m}.csv")

        nn = NearestNeighbors(n_neighbors=K + 1).fit(e).kneighbors(e, return_distance=False)[:, 1:]
        other = (stage[nn] != stage[:, None]).mean(1) * 100
        for l in list(dict.fromkeys(order15 + order22)):
            s = lab == l
            if s.any():
                mix_rows.append({"method": name, "label": l, "cells": int(s.sum()), "other_stage_pct": other[s].mean()})
        for st, s in [("D15 (all)", is15), ("D22 (all)", ~is15)]:
            mix_rows.append({"method": name, "label": st, "cells": int(s.sum()), "other_stage_pct": other[s].mean()})
        same_stage = stage[nn] == stage[:, None]
        same_lab = (lab[nn] == lab[:, None]) & same_stage
        kept = np.where(same_stage.sum(1) > 0, same_lab.sum(1) / np.maximum(same_stage.sum(1), 1) * 100, np.nan)
        same_geno = (geno[nn] == geno[:, None]) & same_stage
        gk = np.where(same_stage.sum(1) > 0, same_geno.sum(1) / np.maximum(same_stage.sum(1), 1) * 100, np.nan)
        for st, s in [("D15", is15), ("D22", ~is15)]:
            keep_rows.append({"method": name, "stage": st, "same_label_pct": np.nanmean(kept[s])})
            for g in ["ctr", "KO"]:
                geno_rows.append({"method": name, "stage": st, "genotype": g,
                                  "same_genotype_pct": np.nanmean(gk[s & (geno == g)])})
        # D15 -> D22: most common D22 label among each D15 cell's D22 neighbours
        best = []
        for i in np.where(is15)[0]:
            nb = nn[i][~is15[nn[i]]]
            best.append(pd.Series(lab[nb]).value_counts().index[0] if len(nb) else "none")
        ct = pd.crosstab(pd.Series(lab[is15], name="D15 draft label"), pd.Series(best, name="nearest D22 label"),
                         normalize="index") * 100
        ct = ct.reindex(index=[l for l in order15 if l in ct.index],
                        columns=[c for c in order22 + ["none"] if c in ct.columns]).round(0)
        n15 = pd.Series(lab[is15]).value_counts()
        ct.insert(0, "cells", n15.reindex(ct.index).astype(int))
        ct.to_csv(OUT / f"d15_to_d22_{m}.csv"); d15[m] = ct

    mix = pd.DataFrame(mix_rows).pivot(index="label", columns="method", values="other_stage_pct")
    mix = mix.reindex(index=["D15 (all)", "D22 (all)"] + [l for l in order15 + order22 if l in mix.index],
                      columns=[x for _, x in METHODS]).round(1)
    cells = pd.Series({r["label"]: r["cells"] for r in mix_rows})
    mix.insert(0, "cells", cells.reindex(mix.index).astype(int))
    mix.to_csv(OUT / "metrics_stage_mixing.csv")
    keep = pd.DataFrame(keep_rows).pivot(index="stage", columns="method", values="same_label_pct")[[x for _, x in METHODS]].round(1)
    keep.to_csv(OUT / "metrics_label_kept.csv")
    gk = pd.DataFrame(geno_rows).pivot_table(index=["stage", "genotype"], columns="method",
                                             values="same_genotype_pct")[[x for _, x in METHODS]].round(1)
    gk.to_csv(OUT / "metrics_genotype_kept.csv")

    # ---- overview: rows = methods; stage | genotype | D22 labels | D15 labels
    fig, axs = plt.subplots(len(METHODS), 4, figsize=(27, 5.6 * len(METHODS)),
                            gridspec_kw={"wspace": 0.55, "hspace": 0.12})
    for r, (m, name) in enumerate(METHODS):
        u = um[m]
        cat_panel(axs[r, 0], u, obs.stage, STAGE_COL, f"{name}: stage")
        cat_panel(axs[r, 1], u, obs.genotype, GENO_COL, f"{name}: genotype")
        cat_panel(axs[r, 2], u, obs.label, pal22, f"{name}: D22 labels (D15 grey)", order=order22)
        cat_panel(axs[r, 3], u, obs.label, pal15, f"{name}: D15 draft labels (D22 grey)", order=order15)
    fig.savefig(OUT / "umap_overview.png", dpi=70, bbox_inches="tight"); plt.close(fig)

    # ---- marker genes on each method's UMAP
    ncol = 6
    nrow = int(np.ceil(len(genes) / ncol))
    for m, name in METHODS:
        u = um[m]
        fig, axs = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 4.1 * nrow), squeeze=False)
        for k, (ax, g) in enumerate(zip(axs.flat, genes)):
            v = X[:, k]
            vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1
            o = np.argsort(v, kind="stable")
            s_ = ax.scatter(u[o, 0], u[o, 1], c=v[o], s=1.5, cmap=FEAT_CMAP, vmin=0, vmax=vmax, linewidths=0,
                            rasterized=True)
            plt.colorbar(s_, ax=ax, fraction=0.04, pad=0.01)
            style(ax, f"{g}  (D15 {100 * (v[is15] > 0).mean():.0f}% / D22 {100 * (v[~is15] > 0).mean():.0f}% > 0)")
        for ax in list(axs.flat)[len(genes):]:
            ax.axis("off")
        fig.suptitle(f"{name}: marker genes (log1p, all cells)", fontsize=14)
        fig.tight_layout(); fig.savefig(OUT / f"genes_{m}.png", dpi=65); plt.close(fig)

    comp = pd.read_csv(OUT / "composition.csv", index_col=[0, 1])
    sections = "".join(
        f"<h3>{name}</h3><div class=tw>{d15[m].to_html(border=0)}</div>" for m, name in METHODS)
    genes_html = "".join(
        f"<h3>{name}</h3><img src='data:image/png;base64,{b64(OUT / f'genes_{m}.png')}' alt='{name} marker genes'>"
        for m, name in METHODS)
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>D15 D22 Integration</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1400px;margin:0 auto;padding:24px 16px 60px}} h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:34px 0 8px}}
.sub{{color:var(--mut)}} .tw{{overflow-x:auto;margin:8px 0}}
table{{border-collapse:collapse;font-size:12.5px}} th,td{{padding:3px 7px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
th:first-child,td:first-child{{text-align:left}} img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>D15 + D22 endocrine cells{f" (subset: {SUBSET})" if SUBSET else ""}: batch correction on stage</h1>
{f"<p class=sub>Subset {SUBSET}: only the labels in study_config/annotation/integration_subset_{SUBSET}.txt; SCT B, Harmony, RPCA and scVI all refitted on these cells.</p>" if SUBSET else ""}
<p class=sub>{n:,} cells: D15 endocrine subset (7.1d, {int(is15.sum()):,} cells, draft r0.6 labels) and D22 endocrine
cells (7.1h, every label except non-endocrine, {int((~is15).sum()):,} cells). SoupX counts; SCT B on the merged cells
(SCTransform v2, S &minus; G2M regressed, 3,000 variable genes without MT- genes). Batch = stage only; genotype not
corrected. Harmony on the merged model's PCs 1&ndash;30; RPCA with one SCT B model per stage and the merged model's
3,000 genes as anchor features; scVI on the counts of the same genes (batch = stage, S &minus; G2M covariate).
30 dimensions, 15 neighbours, UMAP seed 1234 for every method. <b>No clustering.</b></p>
<h2>Cells</h2><div class=tw>{comp.to_html(border=0)}</div>
<h2>UMAPs</h2><p class=sub>Columns: stage, genotype, D22 labels (shared colours), D15 draft labels.</p>
<img src="data:image/png;base64,{b64(OUT / 'umap_overview.png')}" alt="UMAPs of the four methods by stage, genotype and labels">
<h2>Stage mixing</h2><p class=sub>% of each cell's {K} nearest neighbours from the other stage, mean per label. Fully
mixed: about {100 * share['D15']:.0f}% for D15 cells and {100 * share['D22']:.0f}% for D22 cells. Real stage-specific
states should stay low.</p><div class=tw>{mix.to_html(border=0)}</div>
<h2>Cell types kept apart</h2><p class=sub>% of a cell's same-stage neighbours with the same label (D22 labels;
D15 draft labels).</p><div class=tw>{keep.to_html(border=0)}</div>
<h2>Genotype kept</h2><p class=sub>% of a cell's same-stage neighbours with the same genotype (not corrected, so it
should stay near the uncorrected value).</p><div class=tw>{gk.to_html(border=0)}</div>
<h2>Where each D15 state lands</h2><p class=sub>For each D15 cell, the most common D22 label among its D22 neighbours
(of {K}); % of the D15 label's cells. "none" = no D22 neighbour.</p>{sections}
<h2>Marker genes</h2><p class=sub>log1p (SoupX counts normalised to 10,000), unregressed; one colour scale per gene,
capped at the 99th percentile of expressing cells.</p>{genes_html}
</main></body></html>"""
    (OUT / f"D15_D22_endocrine_integration{'_' + SUBSET if SUBSET else ''}.html").write_text(html)
    print(mix.to_string()); print(keep.to_string()); print(gk.to_string())
    for m, _ in METHODS:
        print(m); print(d15[m].to_string())
    print("wrote", OUT)


if __name__ == "__main__":
    main()
