"""Step 7.1j: cluster overview of one day split by genotype, to tell which cluster is which (after 7.1b).

Request (user, 2026-11-19, D40): the r0.2 UMAP split by genotype; feature plots of the marker genes found at the
highest resolution, split by genotype in the 7.1i style (control | KO on one colour scale per gene; all cells drawn,
the other genotype in a very light grey, lighter than the shown genotype's non-expressing cells); then the UMAP at
every resolution, split by genotype, to pick a resolution.

Marker genes: the top N_TOP Wilcoxon markers of each cluster at the highest 7.1b resolution (markers_r{max}.csv),
duplicates dropped; each fine cluster is grouped under the r0.2 cluster that holds most of its cells.
Inputs: ${OUT_DIR}/day_analysis/{day}/clustering/ (clusters.csv, markers_r*.csv, composition_r*.csv; 7.1b),
        cells_obs.tsv.gz (genotype; 7.1), {day}_merged.h5ad (X = unregressed log1p), sct/sct_umaps.csv.gz (SCT B UMAP).
Writes ${OUT_DIR}/day_analysis/{day}/clustering/overview/ (reruns overwrite):
  umap_r0.2_split.png        r0.2 clusters: all cells | control | KO, numbers at centres
  feature_c{k}.png           marker genes of the fine clusters inside r0.2 cluster k, control | KO
  umap_r{res}_split.png      every resolution: all cells | control | KO
  marker_genes.csv           gene, fine cluster, r0.2 cluster, % control / KO cells > 0
  feature_requested.png      genes asked for by name (optional second argument), control | KO, same style
  {day}_cluster_overview.html
Run: sbatch analysis/scripts/7.1j_cluster_overview.sb D40 [GENE,GENE,...]
     (D40 2026-11-19: SST,GHRL, user request, to see the delta-like and GHRL+ cells)
"""
import base64
import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR = Path(os.environ["OUT_DIR"])
DAY = sys.argv[1]
EXTRA = sys.argv[2].split(",") if len(sys.argv) > 2 and sys.argv[2] else []
DD = OUT_DIR / "day_analysis" / DAY
CL_DIR = DD / "clustering"
OUT = CL_DIR / "overview"
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
BASE, FINE = RESOLUTIONS[0], RESOLUTIONS[-1]
N_TOP = 2
SEED = 1234
GREY = "#dedede"
ZERO_GREY = "#c4c4c4"     # shown genotype, gene not detected (FEAT_CMAP at 0)
OTHER_GREY = "#f0f0f0"    # other genotype: lighter than ZERO_GREY (7.1i style)
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", [ZERO_GREY, "#fdbb84", "#e34a33", "#7f0000"])


def key(r):
    return "leiden_r" + str(r).replace(".", "p")


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def palette(cats):
    a, b = plt.get_cmap("tab20"), plt.get_cmap("tab20b")
    return {c: matplotlib.colors.to_hex((a if i < 20 else b)(i % 20)) for i, c in enumerate(cats)}


def umap_split(u, cl, ko, path, title):
    cats = sorted(cl.unique(), key=int)
    pal = palette(cats)
    col = cl.map(pal).to_numpy()
    fig, axs = plt.subplots(1, 3, figsize=(24, 8))
    for ax, (name, m) in zip(axs, [("all cells", np.ones(len(cl), bool)), ("control", ~ko), ("ISL1 KO", ko)]):
        ax.scatter(u[~m, 0], u[~m, 1], c=GREY, s=3, linewidths=0, rasterized=True)
        idx = np.random.default_rng(SEED).permutation(np.where(m)[0])
        ax.scatter(u[idx, 0], u[idx, 1], c=col[idx], s=4, linewidths=0, rasterized=True)
        for c in cats:
            sel = m & (cl == c).to_numpy()
            if sel.sum() >= 15:
                xy = np.median(u[sel], axis=0)
                ax.text(xy[0], xy[1], c, fontsize=11, weight="bold", ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
        ax.set_title(f"{title}: {name} ({m.sum():,} cells)" + ("" if name == "all cells" else "; other genotype grey"),
                     fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    for c in cats:
        sel = (cl == c).to_numpy()
        axs[0].scatter([], [], color=pal[c], s=50, label=f"{c}  (ctr {int((sel & ~ko).sum()):,} / KO {int((sel & ko).sum()):,})")
    axs[0].legend(loc="lower left", fontsize=9 if len(cats) <= 12 else 7.5, frameon=False,
                  ncol=1 if len(cats) <= 12 else 2)
    fig.tight_layout(); fig.savefig(path, dpi=90); plt.close(fig)


def feature_split(u, X, rows, ko, path, title):
    """rows: (gene, panel label); two genes per figure row: [g1 ctr, g1 KO, g2 ctr, g2 KO]."""
    nrow = int(np.ceil(len(rows) / 2))
    H = 4.1 * nrow + 0.8
    fig, axs = plt.subplots(nrow, 4, figsize=(19, H), squeeze=False,
                            gridspec_kw={"wspace": 0.05, "hspace": 0.38, "top": 1 - 0.8 / H, "bottom": 0.01,
                                         "left": 0.01, "right": 0.97})
    for k, (g, lab) in enumerate(rows):
        v = X[g]
        vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1
        for j, (gname, m) in enumerate([("control", ~ko), ("KO", ko)]):
            ax = axs[k // 2, 2 * (k % 2) + j]
            ax.scatter(u[~m, 0], u[~m, 1], c=OTHER_GREY, s=2, linewidths=0, rasterized=True)
            idx = np.where(m)[0]
            idx = idx[np.argsort(v[idx], kind="stable")]
            s_ = ax.scatter(u[idx, 0], u[idx, 1], c=v[idx], s=3, cmap=FEAT_CMAP, vmin=0, vmax=vmax,
                            linewidths=0, rasterized=True)
            ax.set_title(f"{g} ({lab}) — {gname}\n{100 * (v[m] > 0).mean():.1f}% of {gname} cells > 0", fontsize=10.5)
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            if j == 1:
                plt.colorbar(s_, ax=ax, fraction=0.035, pad=0.01)
    for ax in list(axs.flat)[2 * len(rows):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=13, y=1 - 0.25 / H)
    fig.savefig(path, dpi=80); plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in OUT.glob("*"):
        p.unlink()
    cl = pd.read_csv(CL_DIR / "clusters.csv", index_col=0, dtype=str)
    geno = pd.read_csv(DD / "cells_obs.tsv.gz", sep="\t", index_col=0, usecols=[0, 2]).reindex(cl.index)
    assert geno.notna().all().all()
    ko = geno.iloc[:, 0].astype(str).str.upper().str.startswith("KO").to_numpy()
    um = pd.read_csv(DD / "sct" / "sct_umaps.csv.gz", index_col=0).reindex(cl.index)
    assert um.notna().all().all()
    u = um[["SCT_B_UMAP1", "SCT_B_UMAP2"]].to_numpy()

    # fine clusters -> their main r0.2 cluster; top markers
    base, fine = cl[key(BASE)], cl[key(FINE)]
    parent = pd.crosstab(fine, base).idxmax(axis=1)
    share = (pd.crosstab(fine, base).max(axis=1) / fine.value_counts()).round(2)
    mk = pd.read_csv(CL_DIR / f"markers_r{FINE}.csv", dtype={"group": str})
    rows, seen = [], set()
    for b in sorted(base.unique(), key=int):
        for f in sorted(parent.index[parent == b], key=int):
            for g in mk.loc[mk.group == f, "names"].head(N_TOP):
                if g not in seen:
                    seen.add(g); rows.append((g, f, b))

    a = ad.read_h5ad(DD / f"{DAY}_merged.h5ad", backed="r")
    pos = pd.Index(a.obs_names).get_indexer(cl.index)
    assert (pos >= 0).all()
    missing = [g for g in EXTRA if g not in a.var_names]
    if missing:
        sys.exit(f"genes not found: {missing}")
    genes = list(dict.fromkeys([g for g, _, _ in rows] + EXTRA))
    Xm = a[:, [a.var_names.get_loc(g) for g in genes]].X
    Xm = Xm.toarray() if hasattr(Xm, "toarray") else np.asarray(Xm)
    a.file.close()
    X = {g: Xm[pos, i] for i, g in enumerate(genes)}
    pd.DataFrame([{"gene": g, f"r{FINE}_cluster": f, f"r{BASE}_cluster": b,
                   "pct_ctr_gt0": round(100 * (X[g][~ko] > 0).mean(), 1),
                   "pct_KO_gt0": round(100 * (X[g][ko] > 0).mean(), 1)} for g, f, b in rows]
                 ).to_csv(OUT / "marker_genes.csv", index=False)

    umap_split(u, base, ko, OUT / f"umap_r{BASE}_split.png", f"{DAY} r{BASE}")
    comp = {r: pd.read_csv(CL_DIR / f"composition_r{r}.csv", index_col=0) for r in RESOLUTIONS}
    cols = ["ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2", "cells", "pct_of_ctr_cells", "pct_of_KO_cells", "KO_over_ctr"]

    sec = []
    for b in sorted(base.unique(), key=int):
        rb = [(g, f"r{FINE} c{f}") for g, f, bb in rows if bb == b]
        fines = sorted(parent.index[parent == b], key=int)
        nb = int((base == b).sum())
        head = (f"<h3>r{BASE} cluster {b} <span class=sub>({nb:,} cells: control {int(((base == b).to_numpy() & ~ko).sum()):,}"
                f" / KO {int(((base == b).to_numpy() & ko).sum()):,}); r{FINE} clusters inside: "
                + ", ".join(f"{f} ({share[f]:.0%} of its cells here)" for f in fines) + "</span></h3>")
        if rb:
            p = OUT / f"feature_c{b}.png"
            feature_split(u, X, rb, ko, p, f"{DAY}: markers of the r{FINE} clusters inside r{BASE} cluster {b} "
                                           "(log1p; one colour scale per gene; other genotype very light grey)")
            sec.append(head + f'<img src="data:image/png;base64,{b64(p)}" alt="r{BASE} cluster {b} marker genes by genotype">')
        else:
            sec.append(head + "<p class=sub>no fine cluster has most of its cells here</p>")

    req, req_nav = "", ""
    if EXTRA:
        p = OUT / "feature_requested.png"
        feature_split(u, X, [(g, "requested") for g in EXTRA], ko, p,
                      f"{DAY}: requested genes (log1p; one colour scale per gene; other genotype very light grey)")
        req = (f'<h2 id="requested">Requested genes: {", ".join(EXTRA)}</h2>'
               f'<img src="data:image/png;base64,{b64(p)}" alt="{DAY} requested genes by genotype">')
        req_nav = '<a href="#requested">requested genes</a>'
    res_sec = []
    for r in RESOLUTIONS:
        p = OUT / f"umap_r{r}_split.png"
        umap_split(u, cl[key(r)], ko, p, f"{DAY} r{r}")
        res_sec.append(f'<h3 id="r{r}">r{r} <span class=sub>({cl[key(r)].nunique()} clusters)</span></h3>'
                       f'<img src="data:image/png;base64,{b64(p)}" alt="{DAY} clusters at resolution {r} by genotype">'
                       f'<div class=tw>{comp[r][cols].to_html(border=0)}</div>')

    n_ctr, n_ko = int((~ko).sum()), int(ko.sum())
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{DAY} Cluster Overview</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1240px;margin:0 auto;padding:24px 16px 60px}} h1{{font-size:24px;margin:0 0 4px}}
h2{{font-size:20px;margin:40px 0 8px;border-top:1px solid var(--line);padding-top:16px}} h3{{font-size:16px;margin:26px 0 4px}}
.sub{{color:var(--mut);font-weight:normal}} .tw{{overflow-x:auto;margin:8px 0}} nav a{{margin-right:12px}}
table{{border-collapse:collapse;font-size:13px}} th,td{{padding:3px 8px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>{DAY}: which cluster is which</h1>
<p class=sub>{len(cl):,} cells ({n_ctr:,} control, {n_ko:,} ISL1 KO); per-day merge, no batch correction; SCT B graph and UMAP
(7.1c / 7.1b). Gene values are unregressed log1p; each gene has one colour scale for control and KO, capped at the
99th percentile of expressing cells. Marker genes: the top {N_TOP} Wilcoxon markers of each r{FINE} cluster (duplicates
dropped), grouped under the r{BASE} cluster that holds most of that cluster's cells.</p>
<nav><a href="#base">r{BASE} by genotype</a>{req_nav}<a href="#markers">marker genes</a><a href="#res">all resolutions</a></nav>
<h2 id="base">r{BASE} clusters by genotype</h2>
<img src="data:image/png;base64,{b64(OUT / f'umap_r{BASE}_split.png')}" alt="{DAY} r{BASE} clusters by genotype">
<div class=tw>{comp[BASE][cols].to_html(border=0)}</div>
{req}
<h2 id="markers">Marker genes (top {N_TOP} of each r{FINE} cluster), control | KO</h2>
<p class=sub>Other genotype in very light grey; grey = cells of the shown genotype without the gene.</p>
{''.join(sec)}
<h2 id="res">All resolutions, by genotype</h2>
<nav>{''.join(f'<a href="#r{r}">r{r}</a>' for r in RESOLUTIONS)}</nav>
{''.join(res_sec)}
</main></body></html>"""
    out = OUT / f"{DAY}_cluster_overview.html"
    out.write_text(html)
    print(f"{len(rows)} marker genes; wrote", out, f"{out.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
