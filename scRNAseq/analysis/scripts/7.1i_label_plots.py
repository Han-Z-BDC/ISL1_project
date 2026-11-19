"""Step 7.1i: UMAPs of one day's final labels and of key marker genes, split by genotype (after 7.1h).

Request (user, 2026-11-19): D22 UMAP with the new annotation split by genotype, and key marker genes split by
genotype.

Inputs: ${OUT_DIR}/day_analysis/{day}/labels/labels.csv (7.1h), {day}_merged.h5ad (X = unregressed log1p, 7.1),
        sct/sct_umaps.csv.gz (SCT B UMAP, 7.1c; the UMAP 7.1b and 7.1h use);
        label order and marker genes from study_config/annotation/{day}_annotation.tsv (+ MARKERS below).
Writes ${OUT_DIR}/day_analysis/{day}/labels/ (reruns overwrite these files only):
  umap_labels_split.png     labels on the UMAP: all cells | control | KO (other genotype grey), names at centres
  feature_split_markers_own.png  key marker genes, control | KO panels with one colour scale per gene (log1p,
                            capped at the 99th percentile of expressing cells of both genotypes together);
                            each panel shows only its genotype's cells
  feature_split_markers_bg.png   the same, with the other genotype's cells drawn in a very light grey (lighter
                            than the shown genotype's non-expressing cells)
  {day}_labels.html         both figures + label composition table
Run: sbatch analysis/scripts/7.1i_label_plots.sb D22 [GENE,GENE,...]   (genes default to MARKERS[day])
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

OUT_DIR, REPO_DIR = Path(os.environ["OUT_DIR"]), Path(os.environ["REPO_DIR"])
SEED = 1234
GREY = "#dedede"
ZERO_GREY = "#c4c4c4"     # shown genotype, gene not detected (FEAT_CMAP at 0)
OTHER_GREY = "#f0f0f0"    # other genotype in the "bg" version: lighter than ZERO_GREY
FEAT_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("feat", [ZERO_GREY, "#fdbb84", "#e34a33", "#7f0000"])
# Key markers per label (user decisions in DECISIONS.md / the annotation file), in label order.
MARKERS = {
    "D22": [("SC_alpha", ["ARX", "GCG", "IRX2"]), ("SC_beta", ["INS", "ISL1", "PCP4", "NKX6-1"]),
            ("SC_delta", ["SST", "HHEX"]), ("ENP", ["NEUROG3", "GHRL", "CDKN1C"]),
            ("SC_EC1", ["PLCL1", "RUNX1T1"]), ("SC_EC2", ["TPH1", "LMX1A", "FEV"]),
            ("non-endocrine", ["AFP", "SOX9", "PDX1"]), ("all endocrine", ["CHGA"])],
}


def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def main():
    day = sys.argv[1]
    dd = OUT_DIR / "day_analysis" / day
    out = dd / "labels"
    lab = pd.read_csv(out / "labels.csv", index_col=0)
    rules = pd.read_csv(REPO_DIR / "study_config" / "annotation" / f"{day}_annotation.tsv", sep="\t", dtype=str)
    order = list(dict.fromkeys(rules.label))
    um = pd.read_csv(dd / "sct" / "sct_umaps.csv.gz", index_col=0).reindex(lab.index)
    assert um.notna().all().all()
    u = um[["SCT_B_UMAP1", "SCT_B_UMAP2"]].to_numpy()
    ko = lab.genotype.astype(str).str.upper().str.startswith("KO").to_numpy()
    n_ctr, n_ko = int((~ko).sum()), int(ko.sum())
    labels = lab.label.astype(str)
    pal = {l: matplotlib.colors.to_hex(plt.get_cmap("tab10")(i % 10)) for i, l in enumerate(order)}
    col = labels.map(pal).to_numpy()

    # ---- labels: all | control | KO
    fig, axs = plt.subplots(1, 3, figsize=(24, 8))
    for ax, (name, m) in zip(axs, [("all cells", np.ones(len(lab), bool)), ("control", ~ko), ("ISL1 KO", ko)]):
        ax.scatter(u[~m, 0], u[~m, 1], c=GREY, s=3, linewidths=0, rasterized=True)
        idx = np.random.default_rng(SEED).permutation(np.where(m)[0])
        ax.scatter(u[idx, 0], u[idx, 1], c=col[idx], s=4, linewidths=0, rasterized=True)
        for l in order:
            sel = m & (labels == l).to_numpy()
            if sel.sum() >= 20:
                xy = np.median(u[sel], axis=0)
                ax.text(xy[0], xy[1], l, fontsize=11, weight="bold", ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75))
        ax.set_title(f"{day} {name} ({m.sum():,} cells)" + ("" if name == "all cells" else "; other genotype grey"),
                     fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    for l in order:
        n_c = int(((labels == l).to_numpy() & ~ko).sum()); n_k = int(((labels == l).to_numpy() & ko).sum())
        axs[0].scatter([], [], color=pal[l], s=50, label=f"{l}  (ctr {n_c:,} / KO {n_k:,})")
    axs[0].legend(loc="lower left", fontsize=10, frameon=False, markerscale=1.2)
    fig.tight_layout(); fig.savefig(out / "umap_labels_split.png", dpi=100); plt.close(fig)

    # ---- marker genes: one row per gene, control | KO, shared scale
    if len(sys.argv) > 2:
        groups = [("requested", sys.argv[2].split(","))]
    else:
        groups = MARKERS[day]
    a = ad.read_h5ad(dd / f"{day}_merged.h5ad", backed="r")
    pos = pd.Index(a.obs_names).get_indexer(lab.index)
    assert (pos >= 0).all()
    rows = [(g, grp) for grp, gs in groups for g in gs if g in a.var_names]
    missing = [g for _, gs in groups for g in gs if g not in a.var_names]
    if missing:
        print("not in the data:", missing)
    gidx = [a.var_names.get_loc(g) for g, _ in rows]
    X = a[:, gidx].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    X = X[pos]
    a.file.close()
    pct = {}
    for k, (g, grp) in enumerate(rows):
        v = X[:, k]
        pct[g] = (grp, round(100 * (v[~ko] > 0).mean(), 1), round(100 * (v[ko] > 0).mean(), 1))
    # same axis limits in every panel, so the two genotypes line up whether or not the other one is drawn
    pad = 0.03 * (u.max(axis=0) - u.min(axis=0))
    lim = (u.min(axis=0) - pad, u.max(axis=0) + pad)
    # Two versions (user, 2026-11-19):
    #   own:  only the shown genotype's cells
    #   bg:   all cells; the other genotype in a very light grey, lighter than the shown genotype's
    #         non-expressing cells (FEAT_CMAP's 0 = ZERO_GREY)
    for mode, fname, note in [("own", "feature_split_markers_own.png", "only the shown genotype's cells"),
                              ("bg", "feature_split_markers_bg.png",
                               "other genotype in very light grey; grey = shown genotype, not detected")]:
        # 2 genes per figure row: [gene1 ctr, gene1 KO, gene2 ctr, gene2 KO]
        nrow = int(np.ceil(len(rows) / 2))
        H = 4.1 * nrow + 0.8
        fig, axs = plt.subplots(nrow, 4, figsize=(19, H), squeeze=False,
                                gridspec_kw={"wspace": 0.05, "hspace": 0.35, "top": 1 - 0.8 / H, "bottom": 0.01,
                                             "left": 0.01, "right": 0.97})
        for k, (g, grp) in enumerate(rows):
            v = X[:, k]
            vmax = np.percentile(v[v > 0], 99) if (v > 0).any() else 1
            for j, (gname, m) in enumerate([("control", ~ko), ("KO", ko)]):
                ax = axs[k // 2, 2 * (k % 2) + j]
                if mode == "bg":
                    ax.scatter(u[~m, 0], u[~m, 1], c=OTHER_GREY, s=2, linewidths=0, rasterized=True)
                idx = np.where(m)[0]
                idx = idx[np.argsort(v[idx], kind="stable")]
                sc_ = ax.scatter(u[idx, 0], u[idx, 1], c=v[idx], s=3, cmap=FEAT_CMAP, vmin=0, vmax=vmax,
                                 linewidths=0, rasterized=True)
                ax.set_xlim(lim[0][0], lim[1][0]); ax.set_ylim(lim[0][1], lim[1][1])
                ax.set_title(f"{g} ({grp}) — {gname}\n{pct[g][1 + j]}% of {gname} cells > 0", fontsize=10.5)
                ax.set_xticks([]); ax.set_yticks([])
                for s in ax.spines.values():
                    s.set_visible(False)
                if j == 1:
                    plt.colorbar(sc_, ax=ax, fraction=0.035, pad=0.01)
        for ax in list(axs.flat)[2 * len(rows):]:
            ax.axis("off")
        fig.suptitle(f"{day}: key marker genes by genotype (log1p; same colour scale for control and KO of each gene; "
                     f"{note})", fontsize=13, y=1 - 0.25 / H)
        fig.savefig(out / fname, dpi=85); plt.close(fig)
    old = out / "feature_split_markers.png"          # single version from the first 7.1i run
    if old.exists():
        old.unlink()

    # ---- page
    comp = pd.read_csv(out / "composition.csv", index_col=0)
    pt = pd.DataFrame(pct, index=["label", "% control cells > 0", "% KO cells > 0"]).T
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{day} Labels</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1240px;margin:0 auto;padding:24px 16px 60px}} h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:34px 0 8px}}
.sub{{color:var(--mut)}} .tw{{overflow-x:auto;margin:8px 0}}
table{{border-collapse:collapse;font-size:13px}} th,td{{padding:4px 8px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
th:first-child,td:first-child{{text-align:left}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px;margin:6px 0}}
</style></head><body><main>
<h1>{day}: final labels, split by genotype</h1>
<p class=sub>{len(lab):,} cells ({n_ctr:,} control, {n_ko:,} ISL1 KO). Labels from <code>study_config/annotation/{day}_annotation.tsv</code>
(step 7.1h). SCT B UMAP (7.1c/7.1b); no batch correction. Gene values are unregressed log1p; each gene uses one colour
scale for both genotypes, capped at the 99th percentile of expressing cells.</p>
<h2>Labels</h2>
<img src="data:image/png;base64,{b64(out / 'umap_labels_split.png')}" alt="{day} UMAP coloured by label: all cells, control, KO">
<div class=tw>{comp.to_html(border=0)}</div>
<h2>Key marker genes: control vs KO</h2>
<h3>Version 1: each panel shows only its genotype's cells</h3>
<img src="data:image/png;base64,{b64(out / 'feature_split_markers_own.png')}" alt="{day} marker genes, control and KO panels, own genotype only">
<h3>Version 2: other genotype in very light grey</h3>
<img src="data:image/png;base64,{b64(out / 'feature_split_markers_bg.png')}" alt="{day} marker genes, control and KO panels, other genotype light grey">
<div class=tw>{pt.to_html(border=0)}</div>
</main></body></html>"""
    (out / f"{day}_labels.html").write_text(html)
    print(comp.to_string()); print(pt.to_string()); print("wrote", out)


if __name__ == "__main__":
    main()
