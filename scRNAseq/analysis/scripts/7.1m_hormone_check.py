"""Step 7.1m (hormone check, one-off; user, 2026-11-22: "is INS higher in KO EC at D22?", then "how about GCG?").

Per-cell hormone reads in SC_EC1 / SC_EC2 by genotype, D22 and D40, with the hormone's own cell type as the reference
(INS: SC_beta, GCG: SC_alpha). Per cell: raw counts (layer counts) and SoupX counts (layer soupx_counts) of INS and
GCG, total SoupX counts, genes detected, scDblFinder score, and log1p (X) of program genes used to tell a real
EC cell from a doublet with the hormone cell type: beta (IAPP, G6PC2, MAFA, ERO1B, ADCYAP1, HADH), alpha (ARX, IRX2,
DPP4, ALDH1A1, LOXL4, GC), EC (TPH1, FEV), SST.
Needs 7.1h. Plot: fig_hormone_check.R (local).
Writes ${OUT_DIR}/day_analysis/ec_similarity/hormone_check/cells.csv.gz
Run: sbatch --wrap with ${TOOLS}/python (acpu, 6 CPUs, 22G; OUT_DIR exported)
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

OUT_DIR = Path(os.environ["OUT_DIR"])
OUT = OUT_DIR / "day_analysis" / "ec_similarity" / "hormone_check"
OUT.mkdir(parents=True, exist_ok=True)
HORMONES = ["INS", "GCG"]
GENES = ["INS", "GCG", "IAPP", "G6PC2", "MAFA", "ERO1B", "ADCYAP1", "HADH", "ARX", "IRX2", "DPP4", "ALDH1A1", "LOXL4",
         "GC", "SST", "TPH1", "FEV"]
LABELS = ["SC_EC1", "SC_EC2", "SC_beta", "SC_alpha"]


def col(m, j):
    v = m[:, j]
    return np.asarray(v.toarray()).ravel() if sp.issparse(v) else np.asarray(v).ravel()


rows = []
for day in ["D22", "D40"]:
    dd = OUT_DIR / "day_analysis" / day
    lab = pd.read_csv(dd / "labels" / "labels.csv", index_col=0)
    lab = lab[lab.label.isin(LABELS)]
    a = ad.read_h5ad(dd / f"{day}_merged.h5ad", backed="r")
    idx = a.obs_names.get_indexer(lab.index)
    assert (idx >= 0).all()
    order = np.argsort(idx)
    idx, lab = idx[order], lab.iloc[order]
    sub = a[idx].to_memory()
    gi = {g: sub.var_names.get_loc(g) for g in GENES if g in sub.var_names}
    raw, sx, X = sub.layers["counts"], sub.layers["soupx_counts"], sub.X
    d = pd.DataFrame({"barcode": lab.index, "day": day, "sample": lab["sample"].values,
                      "genotype": np.where(lab.genotype.astype(str).str.upper().str.startswith("KO"), "KO", "ctr"),
                      "label": lab.label.values,
                      "total_soupx": np.asarray(sx.sum(axis=1)).ravel(),
                      "n_genes": np.asarray((sx > 0).sum(axis=1)).ravel(),
                      "doublet_score": sub.obs["scDblFinder_doublet_score"].values})
    for h in HORMONES:
        d[f"{h}_raw"], d[f"{h}_soupx"] = col(raw, gi[h]), col(sx, gi[h])
    for g, j in gi.items():
        d[f"log1p_{g}"] = col(X, j)
    rows.append(d)
    print(day, d.groupby(["label", "genotype"]).size().to_dict())
out = pd.concat(rows)
out.to_csv(OUT / "cells.csv.gz", index=False)
print("wrote", OUT / "cells.csv.gz", len(out))
