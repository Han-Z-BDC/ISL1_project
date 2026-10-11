#!/usr/bin/env python
"""
RNA normalize + UMAP: two parallel single-sample methods (log-normalize,
SCTransform) -- NOT cross-sample/replicate integration, that's a future
step's scope.

Design (compared the lab's own notebook -- 3.single.timepoint.joint.
analysis.ipynb -- against current scanpy and current Seurat docs; every
fork confirmed with the user -- see configs/normalize_rna.yaml for the
full reasoning):
  - "log" path: scanpy-native normalize_total -> log1p -> HVG -> PCA ->
    UMAP -> Leiden. No scaling (matches current scanpy's own tutorial,
    which dropped it -- see config comments), no regressing out
    mt%/ribo%/cell-cycle either (same reasoning, confirmed with the user).
  - "sct" path: Seurat::SCTransform via rpy2 -- the one piece of math with
    no scanpy equivalent to what the notebook actually validated -- then
    RunPCA/RunUMAP in R (matching the notebook's own pipe chain exactly),
    with the results transferred back into this same AnnData object.
    Clustering on this path is still done in scanpy (on the transferred
    sct.pca), not Seurat's FindClusters, so both paths use one consistent
    clustering method.
  - Cell-cycle scores are computed ONCE, in scanpy, using cc.genes'
    s.genes/g2m.genes pulled from the R Seurat package at run time (not
    scanpy's own tutorial's gene list) -- so both paths agree on the exact
    same gene set and scores, whether or not a given path regresses on them.

Input is step 3.1's {sample}.rna_final.h5ad (barcode-intersected with ATAC,
SoupX-corrected, doublet-removed). One sample per invocation -- called by
4.1_rna_normalize.sb as a Slurm array task, same pattern as every other step.

Usage:
  python 4.1_rna_normalize.py --sample-id D22_rep1 \
      --rna-h5ad /path/to/merged/h5ad/D22_rep1.rna_final.h5ad \
      --config configs/normalize_rna.yaml \
      --outdir /path/to/sc_pipeline_out/rna_normalize
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

import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

anndata2ri.activate()
anndata2ri.scipy2ri.activate()


def get_cc_genes():
    """Pulls Seurat's own cc.genes$s.genes / cc.genes$g2m.genes from R --
    the exact gene list the notebook's CellCycleScoring call uses -- so the
    scanpy-side score_genes_cell_cycle call agrees with it exactly, rather
    than using scanpy's own tutorial's (slightly different) list."""
    importr("Seurat")
    s_genes = list(ro.r("Seurat::cc.genes$s.genes"))
    g2m_genes = list(ro.r("Seurat::cc.genes$g2m.genes"))
    return s_genes, g2m_genes


def run_log_path(adata, cfg: dict) -> None:
    """scanpy-native: normalize_total -> log1p -> HVG -> PCA -> UMAP ->
    Leiden. No scaling, no covariate regression -- see module docstring /
    configs/normalize_rna.yaml for why. Writes into adata in place:
    layers["log1p_norm"], obsm["X_pca_log"]/["X_umap_log"], obs["leiden_log"]."""
    log_cfg = cfg["log"]
    ad_log = adata.copy()
    sc.pp.normalize_total(ad_log)
    sc.pp.log1p(ad_log)
    adata.layers["log1p_norm"] = ad_log.X.copy()

    sc.pp.highly_variable_genes(ad_log, n_top_genes=log_cfg["n_top_genes"], flavor=log_cfg["hvg_flavor"])
    adata.var["highly_variable_log"] = ad_log.var["highly_variable"].to_numpy()

    # Explicit HVG subset before PCA (not scanpy's use_highly_variable= flag,
    # which is on a deprecation path) -- same pattern as 1.2_remove_background.py's
    # throwaway clustering scaffold. ad_log_hvg is throwaway; only its PCA
    # embedding gets transferred back to adata.
    ad_log_hvg = ad_log[:, ad_log.var["highly_variable"]].copy()
    sc.tl.pca(ad_log_hvg, random_state=cfg["random_state"])
    adata.obsm["X_pca_log"] = ad_log_hvg.obsm["X_pca"]
    ad_log = ad_log_hvg

    sc.pp.neighbors(
        ad_log, use_rep="X_pca", n_pcs=log_cfg["n_pcs"], n_neighbors=cfg["n_neighbors"],
        random_state=cfg["random_state"], key_added="log",
    )
    sc.tl.umap(ad_log, neighbors_key="log", random_state=cfg["random_state"])
    adata.obsm["X_umap_log"] = ad_log.obsm["X_umap"]

    sc.tl.leiden(
        ad_log, neighbors_key="log", key_added="leiden_log", resolution=cfg["leiden_resolution"],
        random_state=cfg["random_state"], flavor="igraph", n_iterations=2, directed=False,
    )
    adata.obs["leiden_log"] = ad_log.obs["leiden_log"].to_numpy()


def run_sct_path(adata, cfg: dict):
    """Seurat::SCTransform via rpy2, then RunPCA/RunUMAP in R -- matching
    the notebook's own pipe chain exactly -- with results transferred back
    into adata: layers["SCT_data"], obsm["X_pca_sct"]/["X_umap_sct"].
    Leiden clustering is done afterward, in scanpy, on the transferred PCA
    (see run_clustering_on_sct)."""
    sct_cfg = cfg["sct"]
    counts = adata.layers.get("soupx_counts", adata.X)
    gene_names = pd.Index(adata.var_names).str.replace("_", "-")  # Seurat-compatible, local copy only -- adata.var_names itself is untouched

    metadata = adata.obs[["pct_counts_mt", "pct_counts_ribo", "S_score", "G2M_score"]].copy()
    metadata.columns = ["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"]

    ro.globalenv["data_mat"] = counts.T.astype("float32")
    ro.globalenv["gene_names"] = gene_names
    ro.globalenv["cell_names"] = adata.obs_names
    ro.globalenv["metadata"] = metadata
    ro.globalenv["vars_to_regress"] = ro.StrVector(["pct_counts_mt", "pct_counts_ribo", "S.Score", "G2M.Score"])
    ro.r(f"set.seed({cfg['random_state']})")

    importr("Seurat")
    print("  [run_sct_path] building Seurat object...", flush=True)
    ro.r("""
        library(Matrix)
        counts <- as(data_mat, "dgCMatrix")
        rownames(counts) <- gene_names
        colnames(counts) <- cell_names
        seu <- CreateSeuratObject(counts = counts, meta.data = metadata)
    """)
    print("  [run_sct_path] Seurat object built. running SCTransform...", flush=True)
    ro.r(f"""
        seu <- SCTransform(seu, vst.flavor = "{sct_cfg['vst_flavor']}", assay = "RNA",
                            vars.to.regress = vars_to_regress)
    """)
    print("  [run_sct_path] SCTransform done. RunPCA...", flush=True)
    ro.r(f"""
        DefaultAssay(seu) <- "SCT"
        seu <- RunPCA(seu, npcs = {sct_cfg['n_pcs']}, reduction.name = "sct.pca")
    """)
    print("  [run_sct_path] RunPCA done. RunUMAP...", flush=True)
    ro.r(f"""
        seu <- RunUMAP(seu, dims = 1:{sct_cfg['n_pcs']}, reduction = "sct.pca", reduction.name = "sct.umap")
    """)
    print("  [run_sct_path] RunUMAP done.", flush=True)

    # Defensive: confirm R preserved cell order/count before any positional
    # assignment back into adata -- SCTransform/RunPCA aren't expected to
    # drop or reorder cells at this point in the pipeline (QC already
    # happened), but a silent mismatch here would misalign every cell's
    # embedding rather than raise anything. Fail loudly instead, same
    # pattern as 3.1_merge_modalities.py's strip_atac_prefix.
    r_cell_order = list(ro.r("Cells(seu)"))
    if r_cell_order != list(adata.obs_names):
        raise ValueError(
            f"Cell order/count changed inside R during SCTransform/RunPCA "
            f"({len(r_cell_order)} R cells vs. {adata.n_obs} in adata) -- "
            f"positional transfer back to adata would misalign cells. "
            f"First few R cells: {r_cell_order[:5]}; "
            f"first few adata cells: {list(adata.obs_names[:5])}"
        )

    sct_data = anndata2ri.rpy2py(ro.r('as(GetAssayData(seu, assay = "SCT", layer = "data"), "dgCMatrix")'))
    pca_emb = np.asarray(ro.r('Embeddings(seu, reduction = "sct.pca")'))
    umap_emb = np.asarray(ro.r('Embeddings(seu, reduction = "sct.umap")'))

    adata.layers["SCT_data"] = sct_data.T
    adata.obsm["X_pca_sct"] = pca_emb
    adata.obsm["X_umap_sct"] = umap_emb


def run_clustering_on_sct(adata, cfg: dict):
    """Leiden on the transferred sct.pca -- in scanpy, not Seurat's
    FindClusters, so both paths share one clustering method. A throwaway
    AnnData wrapping just the transferred embedding, since sc.pp.neighbors
    needs an AnnData to attach its graph to."""
    ad_sct = adata.copy()
    ad_sct.obsm["X_pca"] = adata.obsm["X_pca_sct"]
    sc.pp.neighbors(
        ad_sct, use_rep="X_pca", n_pcs=cfg["sct"]["n_pcs"], n_neighbors=cfg["n_neighbors"],
        random_state=cfg["random_state"], key_added="sct",
    )
    sc.tl.leiden(
        ad_sct, neighbors_key="sct", key_added="leiden_sct", resolution=cfg["leiden_resolution"],
        random_state=cfg["random_state"], flavor="igraph", n_iterations=2, directed=False,
    )
    adata.obs["leiden_sct"] = ad_sct.obs["leiden_sct"].to_numpy()


def plot_umaps(adata, sample_id: str, out_png: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    panels = [
        (axes[0, 0], "X_umap_log", "leiden_log", "log-normalize: Leiden"),
        (axes[0, 1], "X_umap_log", "Phase", "log-normalize: cell-cycle phase"),
        (axes[1, 0], "X_umap_sct", "leiden_sct", "SCTransform: Leiden"),
        (axes[1, 1], "X_umap_sct", "Phase", "SCTransform: cell-cycle phase"),
    ]
    for ax, obsm_key, color_key, title in panels:
        coords = adata.obsm[obsm_key]
        cats = pd.Categorical(adata.obs[color_key].astype(str))
        for cat in cats.categories:
            m = cats == cat
            ax.scatter(coords[m, 0], coords[m, 1], s=3, alpha=0.7, label=str(cat))
        ax.set_title(f"{sample_id}: {title}")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.legend(markerscale=4, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--rna-h5ad", required=True, help="3.1_merge_modalities.py output: {sample}.rna_final.h5ad")
    ap.add_argument("--config", required=True, help="Path to configs/normalize_rna.yaml")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "summary").mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    adata = sc.read_h5ad(args.rna_h5ad)
    sc.pp.filter_genes(adata, min_cells=cfg["min_cells"])
    print(f"[{args.sample_id}] loaded + gene-filtered: {adata.shape}", flush=True)

    s_genes, g2m_genes = get_cc_genes()
    print(f"[{args.sample_id}] got cc.genes from R: {len(s_genes)} s.genes, {len(g2m_genes)} g2m.genes", flush=True)
    s_genes = [g for g in s_genes if g in adata.var_names]
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]
    ad_for_cc = adata.copy()
    sc.pp.normalize_total(ad_for_cc)
    sc.pp.log1p(ad_for_cc)
    sc.tl.score_genes_cell_cycle(ad_for_cc, s_genes=s_genes, g2m_genes=g2m_genes)
    adata.obs["S_score"] = ad_for_cc.obs["S_score"].to_numpy()
    adata.obs["G2M_score"] = ad_for_cc.obs["G2M_score"].to_numpy()
    adata.obs["Phase"] = ad_for_cc.obs["phase"].to_numpy()

    print(f"[{args.sample_id}] starting log path...", flush=True)
    run_log_path(adata, cfg)
    print(f"[{args.sample_id}] log path done. starting SCT path...", flush=True)
    run_sct_path(adata, cfg)
    print(f"[{args.sample_id}] SCT path done. clustering sct.pca...", flush=True)
    run_clustering_on_sct(adata, cfg)
    print(f"[{args.sample_id}] clustering done.", flush=True)

    if cfg["save_diagnostic_plot"]:
        plot_umaps(adata, args.sample_id, outdir / "plots" / f"{args.sample_id}.rna_normalize_umap.png")

    adata.write_h5ad(outdir / "h5ad" / f"{args.sample_id}.rna_normalized.h5ad")

    summary = {
        "sample_id": args.sample_id,
        "n_cells": adata.n_obs,
        "n_genes_after_filter": adata.n_vars,
        "n_hvgs_log": int(adata.var["highly_variable_log"].sum()),
        "n_leiden_clusters_log": adata.obs["leiden_log"].nunique(),
        "n_leiden_clusters_sct": adata.obs["leiden_sct"].nunique(),
        "phase_counts": adata.obs["Phase"].value_counts().to_dict(),
    }
    with open(outdir / "summary" / f"{args.sample_id}.rna_normalize_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"[{args.sample_id}] RNA normalize done: {adata.n_obs} cells, "
          f"{int(adata.var['highly_variable_log'].sum())} log-path HVGs, "
          f"{adata.obs['leiden_log'].nunique()} log-path clusters, "
          f"{adata.obs['leiden_sct'].nunique()} SCT-path clusters. "
          f"Output written to {outdir}")


if __name__ == "__main__":
    main()
