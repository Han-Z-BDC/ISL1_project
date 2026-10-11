#!/usr/bin/env python
"""
RNA QC step 1.2: ambient RNA removal with SoupX.

Design (agreed lab default -- see configs/remove_background.yaml for the
full three-way comparison against the lab's own notebook, CellCommander/
IGVF-UCSD, and sc-best-practices.org):
  - Contamination estimation is SoupX's fully automatic autoEstCont() --
    NOT CellCommander/IGVF-UCSD's marker-gene-based route. Matches the
    lab's own notebook and sc-best-practices; avoids needing a curated
    marker-gene file for every project.
  - A throwaway Leiden clustering (on a normalized/log1p/HVG-subset copy)
    is used only to group cells for SoupX -- it is never used downstream.
  - roundToInt=TRUE and the raw/unfiltered CellRanger matrix as the soup
    source (tod) are universal across all three references.
  - A mandatory pre/post diagnostic UMAP plot, matching the lab's
    before/after-visualization requirement.

Input is step 1.1's *_qc_filtered.h5ad (the cells) plus the same raw h5
used in step 1 (the full droplet population, i.e. SoupX's "tod"). One
sample per invocation -- called by 1.2_remove_background.sb as a Slurm
array task, same pattern as 1.1_rna_qc.sb.

Usage:
  python 1.2_remove_background.py --sample-id D22_rep1 \
      --filtered-h5ad /path/to/rna_qc/h5ad/D22_rep1.qc_filtered.h5ad \
      --raw-h5-path /path/to/D22_rep1/outs/raw_feature_bc_matrix.h5 \
      --config configs/remove_background.yaml \
      --outdir /path/to/sc_pipeline_out/remove_background
"""
import argparse
import json
import sys
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


def load_config(config_path: str, sample_id: str) -> dict:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    cfg["exclude_preclusters"] = cfg.get("exclude_preclusters_by_sample", {}).get(sample_id)
    return cfg


def build_preclusters(adata_filtered, cfg: dict):
    """Normalize/HVG/PCA/neighbors/leiden/UMAP on a throwaway copy, purely
    to group cells for SoupX. Returns the leiden labels (aligned to
    adata_filtered.obs_names) and the UMAP coords (n_obs x 2)."""
    adata_pp = adata_filtered.copy()
    sc.pp.filter_genes(adata_pp, min_cells=cfg["clustering_min_cells"])
    sc.pp.normalize_total(adata_pp, target_sum=cfg["normalize_target_sum"])
    sc.pp.log1p(adata_pp)
    sc.pp.highly_variable_genes(adata_pp, n_top_genes=cfg["n_hvgs"])
    adata_pp = adata_pp[:, adata_pp.var.highly_variable].copy()
    sc.pp.pca(adata_pp, n_comps=cfg["n_pcs"], random_state=cfg["random_state"])
    sc.pp.neighbors(
        adata_pp, n_neighbors=cfg["n_neighbors"], n_pcs=cfg["n_pcs"],
        random_state=cfg["random_state"],
    )
    # flavor="igraph" is the current non-deprecated leiden call in scanpy --
    # a modernization, not a methods disagreement with the references.
    sc.tl.leiden(
        adata_pp, key_added="pre_soupx_leiden", resolution=cfg["leiden_resolution"],
        random_state=cfg["random_state"], flavor="igraph", n_iterations=2, directed=False,
    )
    sc.tl.umap(adata_pp, min_dist=cfg["umap_min_dist"], random_state=cfg["random_state"])
    return adata_pp.obs["pre_soupx_leiden"].copy(), adata_pp.obsm["X_umap"], adata_pp


def run_soupx(adata_filtered, adata_raw, leiden, umap_coords, cfg: dict):
    """Runs SoupX in R via rpy2 and returns the corrected counts matrix,
    plus the tfidfMin value (or manual contamination fraction, see below)
    that actually worked (for logging).

    Second-tier fallback (cfg["manual_contamination_fallbacks"], default
    [0.8, 0.75]): only reached if EVERY tfidfMin attempt still fails.
    Bypasses autoEstCont's estimation entirely via SoupX's own
    setContaminationFraction(sc, contFrac, forceAccept=TRUE). Confirmed
    real on a handful of batch1 T11/H1 samples after 1.1_rna_qc.py's
    adaptive %MT cutoff started retaining many more high-MT cells than
    the old fixed 3.0% cutoff did -- autoEstCont's own rho estimate on
    those samples came back implausibly high (0.75-0.80) and refused to
    proceed without forceAccept, a different failure mode than the
    tfidfMin ladder above was built for (that one is about gene
    selection, not the accept/reject gate). User-confirmed values, not
    derived from the notebook (this specific failure hadn't come up
    there).

    autoEstCont() is tried with SoupX's OWN default first (tfidfMin=1,
    not passed explicitly) -- NOT overridden as a blanket default. Only on
    a real failure does this fall back, in order, through
    cfg["tfidf_min_fallbacks"] (default [0.8, 0.5, 0.3]), matching the
    lab's own batch-processing notebook's fix (diff.pre.proc.ipynb:
    `autoEstCont(sc, tfidfMin = 0.8)`) -- SoupX's default requires a gene
    to look highly soup-specific to count as a marker, and some samples
    (across a batch this varied) don't have enough such genes to satisfy
    it, so autoEstCont fails outright rather than degrading gracefully.
    Confirmed by the user: don't change the default for every sample,
    only reach for a lower tfidfMin on the specific samples that actually
    need it.

    Deliberately does NOT try to extract a contamination-fraction
    percentage from soupc$fit -- the lab's own notebook never computes
    one either (its soupx_out_dict only has soup_total_counts_removed and
    top_10_genes), and a first real run here crashed on
    exp(coef(soupc$fit)) with "non-numeric argument to mathematical
    function" despite autoEstCont() itself succeeding (it printed
    "Estimated global rho of 0.02" immediately before). Matching the
    notebook's own approach avoids a fragile R-side extraction for a
    metric that was never actually needed -- total_counts_removed
    (computed on the Python side from before/after column sums, no R
    involved) already captures the same signal robustly."""
    metadata = pd.DataFrame(umap_coords, columns=["RD1", "RD2"])
    metadata["Cluster"] = leiden.to_numpy()
    metadata.index = adata_filtered.obs_names

    cells = adata_filtered.obs_names
    genes = adata_filtered.var_names
    data = adata_filtered.X.T.astype("float32")
    data_tod = adata_raw.X.T.astype("float32")
    raw_cells = adata_raw.obs_names
    raw_genes = adata_raw.var_names

    ro.globalenv["metadata"] = metadata
    ro.globalenv["data"] = data
    ro.globalenv["data_tod"] = data_tod
    ro.globalenv["genes"] = genes
    ro.globalenv["cells"] = cells
    ro.globalenv["raw_genes"] = raw_genes
    ro.globalenv["raw_cells"] = raw_cells
    ro.r(f"set.seed({cfg['random_state']})")
    importr("SoupX")

    round_to_int = "TRUE" if cfg["round_to_int"] else "FALSE"
    ro.r(f"""
        library(Matrix)
        data <- as(data, "dgCMatrix")
        data_tod <- as(data_tod, "dgCMatrix")
        rownames(data) <- genes
        colnames(data) <- cells
        rownames(data_tod) <- raw_genes
        colnames(data_tod) <- raw_cells
        soupc = SoupChannel(tod=data_tod, toc=data, calcSoupProfile=TRUE)
        soupc = setDR(soupc, metadata[colnames(soupc$toc), c("RD1", "RD2")])
        soupc = setClusters(soupc, setNames(metadata$Cluster, rownames(metadata)))
    """)

    # SoupX's own default (tfidfMin=1) first, only falling back to
    # progressively lower values -- and only for the samples that actually
    # need it -- on a real autoEstCont failure. See docstring above.
    attempts = [None] + list(cfg.get("tfidf_min_fallbacks", [0.8, 0.5, 0.3]))
    last_err = None
    for tfidf_min in attempts:
        tfidf_arg = f", tfidfMin = {tfidf_min}" if tfidf_min is not None else ""
        try:
            ro.r(f"""
                soupc = autoEstCont(soupc, doPlot = FALSE{tfidf_arg})
                out = adjustCounts(soupc, roundToInt = {round_to_int})
            """)
            if tfidf_min is not None:
                print(f"  [run_soupx] autoEstCont's own default failed; "
                      f"succeeded with tfidfMin={tfidf_min}", flush=True)
            out_py = ro.globalenv["out"]
            return out_py.T, tfidf_min
        except Exception as e:
            last_err = e
            print(f"  [run_soupx] autoEstCont failed with tfidfMin="
                  f"{tfidf_min if tfidf_min is not None else 'default'}: {e}", flush=True)

    # autoEstCont failed at every tfidfMin -- fall back to setting the
    # contamination fraction manually (bypassing autoEstCont's own
    # estimation entirely), at the fractions given here, forceAccept=TRUE.
    # Confirmed by the user for the specific samples that hit this: after
    # the adaptive %MT cutoff (see 1.1_rna_qc.py) started retaining many
    # more previously-excluded high-MT cells, autoEstCont on a handful of
    # batch1 T11/H1 samples started estimating implausibly high
    # contamination (rho 0.75-0.80) that it refuses to accept without
    # forceAccept -- not something a lower tfidfMin fixes, since that
    # ladder only changes which genes feed the estimate, not the
    # accept/reject gate itself.
    manual_fracs = cfg.get("manual_contamination_fallbacks", [0.8, 0.75])
    for frac in manual_fracs:
        try:
            ro.r(f"""
                soupc = setContaminationFraction(soupc, {frac}, forceAccept = TRUE)
                out = adjustCounts(soupc, roundToInt = {round_to_int})
            """)
            print(f"  [run_soupx] autoEstCont failed at every tfidfMin; "
                  f"succeeded with manual contamination fraction {frac} "
                  f"(forceAccept=TRUE)", flush=True)
            out_py = ro.globalenv["out"]
            return out_py.T, f"manual_{frac}"
        except Exception as e:
            last_err = e
            print(f"  [run_soupx] manual contamination fraction {frac} "
                  f"also failed: {e}", flush=True)

    raise RuntimeError(
        f"SoupX's autoEstCont failed at every tfidfMin tried ({attempts}) "
        f"AND every manual contamination fraction tried ({manual_fracs}) -- "
        f"last error: {last_err}"
    )


def plot_pre_post(adata_pp_pre, umap_post, leiden_post, sample_id, out_png: Path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, coords, labels, title in [
        (axes[0], adata_pp_pre.obsm["X_umap"], adata_pp_pre.obs["pre_soupx_leiden"], "Pre-SoupX"),
        (axes[1], umap_post, leiden_post, "Post-SoupX"),
    ]:
        cats = pd.Categorical(labels)
        for cat in cats.categories:
            m = cats == cat
            ax.scatter(coords[m, 0], coords[m, 1], s=3, label=str(cat))
        ax.set_title(f"{sample_id}: {title}")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--filtered-h5ad", required=True, help="Output of 1.1_rna_qc.py: {sample}.qc_filtered.h5ad")
    ap.add_argument("--raw-h5-path", required=True, help="The same raw feature_bc_matrix.h5 used in step 1")
    ap.add_argument("--config", required=True, help="Path to configs/remove_background.yaml")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "summary").mkdir(parents=True, exist_ok=True)

    cfg = load_config(args.config, args.sample_id)

    adata_filtered = sc.read_h5ad(args.filtered_h5ad)
    adata_raw = sc.read_10x_h5(args.raw_h5_path)
    adata_raw.var_names_make_unique()
    if "feature_types" in adata_raw.var.columns:
        adata_raw = adata_raw[:, adata_raw.var["feature_types"] == "Gene Expression"].copy()

    leiden, umap_coords, adata_pp_pre = build_preclusters(adata_filtered, cfg)

    excluded_n = 0
    if cfg["exclude_preclusters"]:
        keep = ~leiden.isin([str(c) for c in cfg["exclude_preclusters"]])
        excluded_n = int((~keep).sum())
        adata_filtered = adata_filtered[keep].copy()
        leiden = leiden[keep]
        umap_coords = umap_coords[keep.to_numpy()]

    corrected, tfidf_min_used = run_soupx(adata_filtered, adata_raw, leiden, umap_coords, cfg)

    adata_filtered.layers["counts"] = adata_filtered.X.copy()
    adata_filtered.layers["soupx_counts"] = corrected
    adata_filtered.X = adata_filtered.layers["soupx_counts"].copy()
    adata_filtered.obs["pre_soupx_leiden"] = leiden.to_numpy()
    adata_filtered.obsm["pre_soupx_umap"] = umap_coords

    cnt_before = adata_filtered.layers["counts"].sum(axis=0).A.squeeze()
    cnt_after = corrected.sum(axis=0).A.squeeze()
    zeroed = cnt_before - cnt_after
    adata_filtered.var["soupx_counts_removed"] = zeroed

    if cfg["save_diagnostic_plot"]:
        # Post-SoupX clustering, for the diagnostic plot only (not used
        # anywhere else -- this is a lightweight re-cluster on the
        # corrected counts, mirroring the lab's own notebook).
        adata_post = adata_filtered.copy()
        adata_post.X = adata_filtered.layers["soupx_counts"]
        sc.pp.filter_genes(adata_post, min_cells=cfg["clustering_min_cells"])
        sc.pp.normalize_total(adata_post, target_sum=cfg["normalize_target_sum"])
        sc.pp.log1p(adata_post)
        sc.pp.highly_variable_genes(adata_post, n_top_genes=cfg["n_hvgs"])
        adata_post = adata_post[:, adata_post.var.highly_variable].copy()
        sc.pp.pca(adata_post, n_comps=cfg["n_pcs"], random_state=cfg["random_state"])
        sc.pp.neighbors(adata_post, n_neighbors=cfg["n_neighbors"], n_pcs=cfg["n_pcs"], random_state=cfg["random_state"])
        sc.tl.leiden(
            adata_post, key_added="post_soupx_leiden", resolution=cfg["leiden_resolution"],
            random_state=cfg["random_state"], flavor="igraph", n_iterations=2, directed=False,
        )
        sc.tl.umap(adata_post, min_dist=cfg["umap_min_dist"], random_state=cfg["random_state"])
        plot_pre_post(
            adata_pp_pre, adata_post.obsm["X_umap"], adata_post.obs["post_soupx_leiden"],
            args.sample_id, outdir / "plots" / f"{args.sample_id}.soupx_pre_post_umap.png",
        )

    adata_filtered.write_h5ad(outdir / "h5ad" / f"{args.sample_id}.soupx_corrected.h5ad")

    top10 = (
        adata_filtered.var.sort_values("soupx_counts_removed", ascending=False)
        .head(10).index.tolist()
    )
    summary = {
        "sample_id": args.sample_id,
        "total_counts_removed": round(float(zeroed.sum()), 2),
        "n_cells_excluded_by_precluster_override": excluded_n,
        "top_10_genes_by_counts_removed": top10,
        "tfidf_min_used": tfidf_min_used if tfidf_min_used is not None else "default(1)",
    }
    with open(outdir / "summary" / f"{args.sample_id}.soupx_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[{args.sample_id}] SoupX: {zeroed.sum():.0f} counts removed. "
          f"Output written to {outdir}")


if __name__ == "__main__":
    main()
