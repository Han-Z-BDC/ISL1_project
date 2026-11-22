"""Figure data export (user, 2026-11-22): per-cell log1p expression of chosen genes for one day, for the local figure
scripts (fig_feature_dot.R). Values = X of {day}_merged.h5ad (unregressed log1p of the SoupX counts, as 7.1i), cells =
the day's 7.1h labels. Genes missing from the data are reported and skipped.
Writes ${OUT_DIR}/day_analysis/figure_data/{day}_genes.csv.gz (barcode + one column per gene)
Run on Alpine: sbatch --wrap with ${TOOLS}/python fig_export_genes.py D22 GCG ARX INS ...
"""
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

OUT_DIR = Path(os.environ["OUT_DIR"])
day, genes = sys.argv[1], sys.argv[2:]
dd = OUT_DIR / "day_analysis" / day
out = OUT_DIR / "day_analysis" / "figure_data"
out.mkdir(parents=True, exist_ok=True)
lab = pd.read_csv(dd / "labels" / "labels.csv", index_col=0)
a = ad.read_h5ad(dd / f"{day}_merged.h5ad", backed="r")
miss = [g for g in genes if g not in a.var_names]
if miss:
    print("not in the data, skipped:", miss)
genes = [g for g in genes if g in a.var_names]
cols = [a.var_names.get_loc(g) for g in genes]
rows = a.obs_names.get_indexer(lab.index)
assert (rows >= 0).all()
X = a.X[np.sort(rows)][:, cols] if not hasattr(a.X, "to_memory") else a.X[np.sort(rows), :][:, cols]
X = X.toarray() if sp.issparse(X) else np.asarray(X)
d = pd.DataFrame(X, columns=genes, index=a.obs_names[np.sort(rows)]).reindex(lab.index)
d.rename_axis("barcode").round(4).to_csv(out / f"{day}_genes.csv.gz")
print(day, d.shape, "->", out / f"{day}_genes.csv.gz")
