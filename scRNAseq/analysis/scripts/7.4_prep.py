"""Step 7.4 (prep): inputs for the trajectory of the D15 + D22 endocrine integration (after 7.2).

Request (user, 2026-11-19): run the 7.3 trajectory set (Slingshot, velocity, CellRank, tradeSeq) on the D15 + D22
endocrine integration; scVI is the chosen integration and its UMAP is the primary UMAP.
  cells      the cells of 7.2 (full run: 21,063 (D15 endocrine 5,994 + D22 endocrine 15,069), labels as in 7.2)
             (D15 states as "D15 <name>", D22 labels as they are)
  embedding  scVI latent (30 dimensions, 7.2c)
  UMAP 2D    the scVI UMAP of the 7.2 report (umap_scvi.csv; 15 neighbours on the scVI latent, seed 1234)
  UMAP 3D    same 15-neighbour graph on the scVI latent, UMAP with 3 components (seed 1234)
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/trajectory/ (reruns overwrite):
  cells.csv          barcode, stage, sample, genotype, label, UMAP1, UMAP2, UMAP3D_1-3
  scvi_latent.csv.gz barcode x scVI_1-30 (same rows as cells.csv)
With INT_SUBSET=<name> (user, 2026-11-19: alpha_beta_EC), everything reads and writes the 7.2 subset folder
integration_D15_D22_endocrine_<name>/ instead: its cells, scVI latent and scVI UMAP.
Run: sbatch analysis/scripts/7.4_prep.sb
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

OUT_DIR = Path(os.environ["OUT_DIR"])
SUBSET = os.environ.get("INT_SUBSET", "")   # 7.2 subset run (see 7.2a); empty = the full integration
INT = OUT_DIR / "day_analysis" / ("integration_D15_D22_endocrine" + (f"_{SUBSET}" if SUBSET else ""))
TR = INT / "trajectory"
SEED = 1234


def main():
    TR.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(INT / "cells.tsv", sep="\t", index_col=0)
    emb = pd.read_csv(INT / "emb_scvi.csv.gz", index_col=0).reindex(obs.index)
    u2 = pd.read_csv(INT / "umap_scvi.csv", index_col=0).reindex(obs.index)
    assert emb.notna().all().all() and u2.notna().all().all() and emb.shape[1] == 30
    a = ad.AnnData(obs=obs[[]].copy())
    a.obsm["E"] = emb.to_numpy()
    sc.pp.neighbors(a, n_neighbors=15, use_rep="E", random_state=SEED)     # as 7.2d
    sc.tl.umap(a, random_state=SEED, n_components=3)
    out = obs[["stage", "sample", "genotype", "label"]].copy()
    out[["UMAP1", "UMAP2"]] = u2[["UMAP1", "UMAP2"]].to_numpy()
    out[["UMAP3D_1", "UMAP3D_2", "UMAP3D_3"]] = a.obsm["X_umap"]
    out.rename_axis("barcode").to_csv(TR / "cells.csv")
    emb.rename_axis("barcode").to_csv(TR / "scvi_latent.csv.gz")
    print(out.groupby(["stage", "label", "genotype"]).size().unstack(fill_value=0).to_string())
    print("wrote", TR)


if __name__ == "__main__":
    main()
