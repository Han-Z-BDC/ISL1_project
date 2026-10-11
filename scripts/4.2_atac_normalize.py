#!/usr/bin/env python
"""
ATAC normalize + UMAP: SnapATAC2 spectral embedding -- NOT cross-sample/
replicate integration, that's a future step's scope.

Design (compared the lab's own notebook -- prepare.peakcall.snapatac2.
ipynb, the actual normalize/dim-reduction notebook -- against the Galaxy/
scverse SnapATAC2 tutorial as an IGVF-UCSD stand-in; confirmed with the
user -- see configs/normalize_atac.yaml for the full reasoning):
  spectral -> umap -> knn -> leiden, on the tile matrix + selected
  features step 3.1's {sample}.atac_final.h5ad already carries forward
  from step 2.2 -- NO fragment reloading, NO add_tile_matrix/
  select_features here, matching the user's own point that re-parsing
  fragments would just redo work 2.2 already did. Overwrites the stale
  spectral/umap/leiden embedding that rode along from 2.2's diagnostic
  pass (computed there on the pre-doublet-removal population) with one
  fit on the actual final, clean cell set.
  - umap/knn exclude spectral dim 0 (matches 2.2's own convention).
  - leiden's objective_function is set explicitly to "modularity"
    (Galaxy/scverse's own docs: CPM gives poorly-separated UMAP clusters)
    -- already this installed snapatac2's default, so this documents the
    choice rather than changing behavior.

Input is step 3.1's {sample}.atac_final.h5ad (barcode-intersected with
RNA, doublet-removed). One sample per invocation -- called by
4.2_atac_normalize.sb as a Slurm array task, same pattern as every other step.

Usage:
  python 4.2_atac_normalize.py --sample-id D22_rep1 \
      --atac-h5ad /path/to/merged/h5ad/D22_rep1.atac_final.h5ad \
      --config configs/normalize_atac.yaml \
      --outdir /path/to/sc_pipeline_out/atac_normalize
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import snapatac2 as snap


def run_normalize(adata, cfg: dict):
    """spectral -> umap -> knn -> leiden, in place. Overwrites
    obsm["X_spectral"]/["X_umap"] and obs["leiden_{resolution}"] -- same
    key names 2.2's diagnostic pass already used, by design (see module
    docstring): this is meant to replace that stale, pre-doublet-removal
    embedding, not create a second one alongside it."""
    snap.tl.spectral(adata, n_comps=cfg["n_comps"], random_state=cfg["random_state"])

    n_dims = adata.obsm["X_spectral"].shape[1]
    use_dims = list(range(1, n_dims)) if cfg["exclude_first_spectral_dim"] else list(range(n_dims))

    snap.tl.umap(adata, use_dims=use_dims, random_state=cfg["random_state"])
    snap.pp.knn(adata, use_rep="X_spectral", use_dims=use_dims, n_neighbors=cfg["n_neighbors"], random_state=cfg["random_state"])

    leiden_key = f"leiden_{cfg['clustering_resolution']}"
    snap.tl.leiden(
        adata, resolution=cfg["clustering_resolution"], objective_function="modularity",
        key_added=leiden_key, random_state=cfg["random_state"],
    )
    return leiden_key


def plot_umap(adata, leiden_key: str, sample_id: str, out_png: Path):
    coords = adata.obsm["X_umap"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    cats = pd.Categorical(adata.obs[leiden_key].astype(str))
    for cat in cats.categories:
        m = cats == cat
        axes[0].scatter(coords[m, 0], coords[m, 1], s=3, alpha=0.7, label=str(cat))
    axes[0].set_title(f"{sample_id}: Leiden ({leiden_key})")
    axes[0].legend(markerscale=4, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")

    for ax, color_key, title in [(axes[1], "tsse", "TSS enrichment"), (axes[2], "n_fragment", "log10(n_fragment + 1)")]:
        vals = adata.obs[color_key].to_numpy()
        if color_key == "n_fragment":
            vals = np.log10(vals + 1)
        sca = ax.scatter(coords[:, 0], coords[:, 1], s=3, alpha=0.7, c=vals, cmap="viridis")
        fig.colorbar(sca, ax=ax, shrink=0.8)
        ax.set_title(f"{sample_id}: {title}")

    for ax in axes:
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--atac-h5ad", required=True, help="3.1_merge_modalities.py output: {sample}.atac_final.h5ad")
    ap.add_argument("--config", required=True, help="Path to configs/normalize_atac.yaml")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "summary").mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    # Plain sc.read_h5ad, not snap.read -- 3.1_merge_modalities.py wrote
    # this file with plain scanpy (atac.write_h5ad), not SnapATAC2's own
    # writer. Confirmed snap.tl.spectral/umap/knn/leiden all accept a
    # regular in-memory anndata.AnnData directly (verified against this
    # exact file before writing this script) -- no snap.read needed.
    adata = sc.read_h5ad(args.atac_h5ad)
    print(f"[{args.sample_id}] loaded: {adata.shape}, {int(adata.var['selected'].sum())} selected features", flush=True)

    leiden_key = run_normalize(adata, cfg)
    print(f"[{args.sample_id}] normalize done: {adata.obs[leiden_key].nunique()} Leiden clusters", flush=True)

    if cfg["save_diagnostic_plot"]:
        plot_umap(adata, leiden_key, args.sample_id, outdir / "plots" / f"{args.sample_id}.atac_normalize_umap.png")

    adata.write_h5ad(outdir / "h5ad" / f"{args.sample_id}.atac_normalized.h5ad")

    summary = {
        "sample_id": args.sample_id,
        "n_cells": adata.n_obs,
        "n_selected_features": int(adata.var["selected"].sum()),
        "n_leiden_clusters": int(adata.obs[leiden_key].nunique()),
    }
    with open(outdir / "summary" / f"{args.sample_id}.atac_normalize_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"[{args.sample_id}] ATAC normalize done: {adata.n_obs} cells, "
          f"{int(adata.obs[leiden_key].nunique())} Leiden clusters. "
          f"Output written to {outdir}")


if __name__ == "__main__":
    main()
