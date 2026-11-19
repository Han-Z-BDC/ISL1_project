"""Step 7.1k: do the labels of two days match? (after 7.1h for both days)

Request (user, 2026-11-19): D22 and D40 are fully labelled; check that the labels given to the same cell type match
between the two days. The days were clustered separately (no shared embedding, no integration), so two checks that
need none:
  1. Marker overlap: top 50 Wilcoxon markers of each label vs the rest of its day (log2 fold change > 0.5,
     adjusted p < 0.05, genes detected in >= 0.5% of the day's cells), counted shared between each D22 x D40 pair.
  2. Profile correlation: mean unregressed log1p per label, on genes highly variable in both days (scanpy "seurat"
     flavour, 3,000 per day, MT- genes removed); each gene centred on the mean of that day's label means, then
     Pearson r between every D22 label and every D40 label. Run on all cells, control cells only and KO cells
     only (a label enters a genotype's run only with >= 20 cells of that genotype).
Same cell type on both days -> high overlap and highest r with each other.

Writes ${OUT_DIR}/day_analysis/compare_{dayA}_{dayB}/ (reruns overwrite):
  markers_{day}.csv              the top-50 marker lists used
  marker_overlap.csv             shared top-50 markers, dayA labels (rows) x dayB labels (columns)
  corr_{all,ctr,KO}.csv          Pearson r, dayA labels x dayB labels
  best_match.csv                 for each label of either day: best partner by overlap and by r (all cells)
  heatmaps.png                   overlap + the three correlation matrices
  {dayA}_vs_{dayB}_labels.html   report
Run: sbatch analysis/scripts/7.1k_label_match.sb D22 D40
"""
import base64
import os
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

OUT_DIR, REPO_DIR = Path(os.environ["OUT_DIR"]), Path(os.environ["REPO_DIR"])
N_TOP, MIN_LFC, MAX_PADJ, MIN_PCT = 50, 0.5, 0.05, 0.005
N_HVG, MIN_CELLS_GT = 3000, 20


def load(day):
    dd = OUT_DIR / "day_analysis" / day
    lab = pd.read_csv(dd / "labels" / "labels.csv", index_col=0)
    a = sc.read_h5ad(dd / f"{day}_merged.h5ad")
    a = a[lab.index].copy()
    a.obs["label"] = pd.Categorical(lab.label.astype(str))
    a.obs["geno"] = np.where(lab.genotype.astype(str).str.upper().str.startswith("KO"), "KO", "ctr")
    rules = pd.read_csv(REPO_DIR / "study_config" / "annotation" / f"{day}_annotation.tsv", sep="\t", dtype=str)
    order = [l for l in dict.fromkeys(rules.label) if l in set(a.obs.label)]
    return a, order


def markers(a, day):
    keep = np.asarray((a.X > 0).sum(axis=0)).ravel() >= MIN_PCT * a.n_obs
    b = a[:, keep].copy()
    sc.tl.rank_genes_groups(b, "label", method="wilcoxon", use_raw=False)
    df = sc.get.rank_genes_groups_df(b, None)
    df = df[(df.logfoldchanges > MIN_LFC) & (df.pvals_adj < MAX_PADJ)]
    df = df.sort_values(["group", "scores"], ascending=[True, False], kind="stable")
    df = df.groupby("group", observed=True, group_keys=False).head(N_TOP)
    df.to_csv(OUT / f"markers_{day}.csv", index=False)
    return {g: set(d.names) for g, d in df.groupby("group", observed=True)}


def hvg(a):
    b = a.copy()
    sc.pp.highly_variable_genes(b, flavor="seurat", n_top_genes=N_HVG)
    return set(b.var_names[b.var.highly_variable & ~b.var_names.str.startswith("MT-")])


def centred_means(a, genes, order, geno=None):
    m = np.ones(a.n_obs, bool) if geno is None else (a.obs.geno == geno).to_numpy()
    out = {}
    for l in order:
        s = m & (a.obs.label == l).to_numpy()
        if s.sum() >= (1 if geno is None else MIN_CELLS_GT):
            x = a[s, genes].X
            out[l] = np.asarray(x.mean(axis=0)).ravel()
    df = pd.DataFrame(out, index=genes)
    return df.sub(df.mean(axis=1), axis=0)


def corr(ma, mb):
    za = (ma - ma.mean()) / ma.std()
    zb = (mb - mb.mean()) / mb.std()
    return (za.T @ zb) / (len(ma) - 1)


def heat(ax, df, title, fmt, cmap, vmin, vmax):
    im = ax.imshow(df.to_numpy(float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(df.shape[1])); ax.set_xticklabels(df.columns, rotation=50, ha="right", fontsize=9)
    ax.set_yticks(range(df.shape[0])); ax.set_yticklabels(df.index, fontsize=9)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            v = df.iat[i, j]
            if pd.notna(v):
                best = v == df.iloc[i].max()
                ax.text(j, i, format(v, fmt), ha="center", va="center", fontsize=7.5,
                        weight="bold" if best else "normal",
                        color="white" if abs(v) > 0.6 * max(abs(vmin), abs(vmax)) else "black")
    ax.set_title(title, fontsize=11)
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)


def main():
    global OUT
    A, B = sys.argv[1], sys.argv[2]
    OUT = OUT_DIR / "day_analysis" / f"compare_{A}_{B}"
    OUT.mkdir(parents=True, exist_ok=True)
    for p in OUT.glob("*"):
        p.unlink()
    a, oa = load(A)
    b, ob = load(B)
    print(A, a.n_obs, oa); print(B, b.n_obs, ob)

    mk_a, mk_b = markers(a, A), markers(b, B)
    ov = pd.DataFrame([[len(mk_a.get(x, set()) & mk_b.get(y, set())) for y in ob] for x in oa], index=oa, columns=ob)
    ov.to_csv(OUT / "marker_overlap.csv")

    genes = sorted((hvg(a) & hvg(b)) & set(a.var_names) & set(b.var_names))
    print("shared HVGs:", len(genes))
    cors = {}
    for name, g in [("all", None), ("ctr", "ctr"), ("KO", "KO")]:
        ma, mb = centred_means(a, genes, oa, g), centred_means(b, genes, ob, g)
        c = corr(ma, mb).reindex(index=oa, columns=ob).round(2)
        c.to_csv(OUT / f"corr_{name}.csv"); cors[name] = c

    rows = []
    for day, other, labels, o_ov, o_c in [(A, B, oa, ov, cors["all"]), (B, A, ob, ov.T, cors["all"].T)]:
        for l in labels:
            r_ov, r_c = o_ov.loc[l], o_c.loc[l]
            rows.append({"day": day, "label": l, "markers": len((mk_a if day == A else mk_b).get(l, set())),
                         "same label in other day": l in (ob if day == A else oa),
                         "best by overlap": r_ov.idxmax() if r_ov.max() > 0 else "-", "shared": int(r_ov.max()),
                         "overlap with same label": int(r_ov[l]) if l in r_ov.index else None,
                         "best by r": r_c.idxmax(), "r": r_c.max(),
                         "r with same label": r_c[l] if l in r_c.index else None})
    bm = pd.DataFrame(rows)
    bm.to_csv(OUT / "best_match.csv", index=False)

    h = max(5, 0.42 * len(oa) + 2)
    fig, axs = plt.subplots(2, 2, figsize=(2 * (0.55 * len(ob) + 3.5), 2 * h))
    heat(axs[0, 0], ov, f"shared top-{N_TOP} markers ({A} rows x {B} columns)", "d", "Greens", 0, N_TOP)
    for ax, (name, t) in zip([axs[0, 1], axs[1, 0], axs[1, 1]],
                             [("all", "all cells"), ("ctr", "control cells"), ("KO", "KO cells")]):
        heat(ax, cors[name], f"Pearson r of centred label means, {t}", ".2f", "RdBu_r", -1, 1)
    fig.tight_layout(); fig.savefig(OUT / "heatmaps.png", dpi=100); plt.close(fig)

    b64 = base64.b64encode((OUT / "heatmaps.png").read_bytes()).decode()
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{A} vs {B} Labels</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1300px;margin:0 auto;padding:24px 16px 60px}} h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:30px 0 8px}}
.sub{{color:var(--mut)}} .tw{{overflow-x:auto;margin:8px 0}}
table{{border-collapse:collapse;font-size:13px}} th,td{{padding:3px 8px;border-bottom:1px solid var(--line);text-align:right;font-variant-numeric:tabular-nums}}
th:first-child,td:first-child{{text-align:left}} img{{max-width:100%;height:auto;background:#fff;border-radius:4px}}
</style></head><body><main>
<h1>{A} vs {B}: do the labels match?</h1>
<p class=sub>Days clustered separately, no integration. (1) Shared top-{N_TOP} markers: Wilcoxon, each label vs the rest of its
day, log2FC &gt; {MIN_LFC}, adjusted p &lt; {MAX_PADJ}. (2) Pearson r between label mean profiles (unregressed log1p) on
{len(genes):,} genes highly variable in both days (MT- removed), each gene centred within its day; control-only and
KO-only runs use labels with &ge; {MIN_CELLS_GT} cells of that genotype. Bold = best match in the row.</p>
<h2>Best match for each label</h2><div class=tw>{bm.to_html(index=False, border=0)}</div>
<h2>Matrices</h2><img src="data:image/png;base64,{b64}" alt="marker overlap and correlation heatmaps">
<h2>Marker overlap table</h2><div class=tw>{ov.to_html(border=0)}</div>
<h2>Correlation, all cells</h2><div class=tw>{cors['all'].to_html(border=0)}</div>
</main></body></html>"""
    (OUT / f"{A}_vs_{B}_labels.html").write_text(html)
    print(bm.to_string()); print(ov.to_string()); print(cors["all"].to_string()); print("wrote", OUT)


if __name__ == "__main__":
    main()
