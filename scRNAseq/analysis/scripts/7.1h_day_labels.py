"""Step 7.1h: per-cell labels for one day from its annotation rules (after 7.1b and any 7.1g splits).

Rules: study_config/annotation/{day}_annotation.tsv (user decisions), columns rule, source, cluster, label.
  rule "base":     source "r<res>" = a 7.1b resolution; every cell gets the label of its cluster there.
  rule "override": source "r<res>_c<cluster>_s<sub_res>" = a 7.1g split; cells of the listed subclusters take this
                   label instead of their earlier label (overrides are applied in file order; an override may list
                   only some subclusters of a split, e.g. a finer split of a group labelled by an earlier override).
Every cell must end up with a label.

Writes ${OUT_DIR}/day_analysis/{day}/labels/ (reruns overwrite):
  labels.csv        barcode, sample, genotype, base cluster, label
  composition.csv   per label: cells per sample, % of all ctr / KO cells, KO_over_ctr
  umap_labels.png   SCT B UMAP coloured by label (all cells, control, KO)
Run: sbatch analysis/scripts/7.1h_day_labels.sb D22
"""
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT_DIR, REPO_DIR = Path(os.environ["OUT_DIR"]), Path(os.environ["REPO_DIR"])
DAY = sys.argv[1]
DAY_DIR = OUT_DIR / "day_analysis" / DAY
OUT = DAY_DIR / "labels"
SEED = 1234
GREY = "#d9d9d9"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in OUT.glob("*"):
        p.unlink()
    rules = pd.read_csv(REPO_DIR / "study_config" / "annotation" / f"{DAY}_annotation.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    cl = pd.read_csv(DAY_DIR / "clustering" / "clusters.csv", index_col=0, dtype=str)
    base = rules[rules.rule == "base"]
    assert base.source.nunique() == 1, "one base resolution"
    bkey = "leiden_r" + base.source.iloc[0][1:].replace(".", "p")
    lab = cl[bkey].map(dict(zip(base.cluster, base.label)))
    for src, g in rules[rules.rule == "override"].groupby("source", sort=False):
        res, c, sr = re.fullmatch(r"r([\d.]+)_c(\d+)_s([\d.]+)", src).groups()
        sub = pd.read_csv(DAY_DIR / "clustering" / f"split_r{res}_c{c}" / "subclusters.csv", index_col=0, dtype=str)
        assert set(g.cluster) <= set(sub[f"s{sr}"]), f"{src}: subcluster not in the split"
        s = sub[f"s{sr}"].map(dict(zip(g.cluster, g.label))).dropna()
        lab.loc[s.index] = s
    assert lab.notna().all(), "cells without a label"

    obs = pd.read_csv(DAY_DIR / "cells_obs.tsv.gz", sep="\t", index_col=0, usecols=[0, 1, 2]).reindex(cl.index)
    assert obs.notna().all().all()
    out = obs.assign(**{bkey: cl[bkey], "label": lab})
    out.rename_axis("barcode").to_csv(OUT / "labels.csv")

    order = list(dict.fromkeys(rules.label))
    comp = pd.crosstab(out.label, out["sample"]).reindex(order)
    comp["cells"] = comp.sum(axis=1)
    ko = out.genotype.astype(str).str.upper().str.startswith("KO")
    n_ctr, n_ko = (~ko).sum(), ko.sum()
    comp["pct_of_ctr_cells"] = (out[~ko].label.value_counts().reindex(order).fillna(0) / n_ctr * 100).round(2)
    comp["pct_of_KO_cells"] = (out[ko].label.value_counts().reindex(order).fillna(0) / n_ko * 100).round(2)
    comp["KO_over_ctr"] = (comp.pct_of_KO_cells / comp.pct_of_ctr_cells).round(2)
    comp.to_csv(OUT / "composition.csv")
    print(comp.to_string())

    um = pd.read_csv(DAY_DIR / "sct" / "sct_umaps.csv.gz", index_col=0).reindex(cl.index)
    u = um[["SCT_B_UMAP1", "SCT_B_UMAP2"]].to_numpy()
    pal = {l: matplotlib.colors.to_hex(plt.get_cmap("tab10")(i % 10)) for i, l in enumerate(order)}
    col = out.label.map(pal).to_numpy()
    fig, axs = plt.subplots(1, 3, figsize=(22, 7.5))
    for ax, (name, m) in zip(axs, [("all cells", np.ones(len(out), bool)), ("control", ~ko.to_numpy()),
                                    ("KO", ko.to_numpy())]):
        ax.scatter(u[~m, 0], u[~m, 1], c=GREY, s=2.5, linewidths=0, rasterized=True)
        idx = np.random.default_rng(SEED).permutation(np.where(m)[0])
        ax.scatter(u[idx, 0], u[idx, 1], c=col[idx], s=3.5, linewidths=0, rasterized=True)
        ax.set_title(f"{DAY} labels: {name} ({m.sum():,} cells)"); ax.set_xticks([]); ax.set_yticks([])
    for l in order:
        axs[0].scatter([], [], color=pal[l], s=40, label=f"{l} ({(out.label == l).sum():,})")
    axs[0].legend(loc="lower left", fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(OUT / "umap_labels.png", dpi=100); plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
