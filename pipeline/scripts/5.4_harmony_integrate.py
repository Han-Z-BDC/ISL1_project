#!/usr/bin/env python
"""
Harmony integration, two strategies x two normalization paths = 4 runs,
all via harmonypy (Python-native port of the same algorithm as R's
harmony package -- NOT Seurat's HarmonyIntegration/IntegrateLayers
wrapper, which is hard-coded to a single per-sample "group" variable and
cannot accept explicit multi-covariate vars_use -- confirmed by reading
Seurat 5.3.0's actual HarmonyIntegration() source).

Strategy A ("explicit"): vars_use=["diff_batch","tech_batch"], a separate
theta per variable -- diff_batch corrected harder (it's the dominant real
batch effect here), tech_batch corrected lightly (confirmed near-negligible).
Strategy B ("per_sample"): vars_use=["sample_id"], harmonypy's own default
theta -- matches the lab's own notebook's actual behavior (Seurat's
HarmonyIntegration always groups by sample/layer identity; this
reproduces that with an equivalent, non-Seurat implementation).

Neither strategy ever includes "stage" as a correction variable -- that's
the biological signal to preserve, not batch noise to remove.

IMPORTANT (from 5.1_coverage_table.py): most stages in this project are
covered by only ONE diff_batch. Neither strategy can manufacture
cross-diff_batch correction evidence where none exists in the data --
look at the per-stage/per-diff_batch coverage before trusting correction
on a sparse stage, regardless of which strategy produced it.

Usage:
  python 5.4_harmony_integrate.py --normalized-rna-h5ad merged_rna_normalized.h5ad \
      --outdir out/ [--theta-diff-batch 3] [--theta-tech-batch 0.5]
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from harmonypy import run_harmony


def harmonize(pca_emb: np.ndarray, meta: pd.DataFrame, vars_use: list, theta, random_state: int):
    ho = run_harmony(pca_emb, meta, vars_use, theta=theta, random_state=random_state)
    # NOTE: this installed harmonypy (2.0.0, C++/pybind11-backed) already
    # returns Z_corr as (N, d) -- its property does the d x N -> N x d
    # transpose internally. This differs from the older, more commonly
    # documented pure-Python harmonypy (0.0.x) API where Z_corr is the raw
    # (d, N) array and callers transpose it themselves -- confirmed by
    # reading this env's actual harmony.py source (Z_corr property returns
    # `self._cpp.Z_corr.T`) and a toy run (`ho.Z_corr.shape == (200, 10)`
    # for N=200, d=10, no extra .T needed). Do NOT add `.T` here.
    corrected = np.asarray(ho.Z_corr)
    if corrected.shape != pca_emb.shape:
        raise ValueError(
            f"harmonypy Z_corr shape {corrected.shape} != input pca_emb shape {pca_emb.shape} "
            f"-- check the installed harmonypy version's Z_corr orientation convention"
        )
    return corrected


def umap_from_pca(obs_names, pca_emb, key, random_state):
    """neighbors/umap only ever touch obsm/obsp -- never .X -- so build a
    minimal AnnData holding just the PCA embedding rather than copying the
    full 68718 x 28843 expression matrix once per run (4 runs total). The
    original adata_template.copy() approach OOM'd at cpus=4/mem=14G on the
    very first run's neighbors/umap call; a full-object copy was the
    actual cause, not neighbors/umap's own cost."""
    a = ad.AnnData(X=np.empty((pca_emb.shape[0], 0), dtype=np.float32), obs=pd.DataFrame(index=obs_names))
    a.obsm["X_pca"] = pca_emb
    sc.pp.neighbors(a, use_rep="X_pca", random_state=random_state, key_added=key)
    sc.tl.umap(a, neighbors_key=key, random_state=random_state)
    return a.obsm["X_umap"]


def plot_comparison(adata, runs: dict, outpath: Path):
    # Every panel gets its own legend (each color_by has a different
    # category set -- diff_batch=4, stage=8, tech_batch=2 -- so a shared
    # legend on just the last column was silently unreadable for the
    # other two columns, not just inconvenient). Placed outside the axes
    # to the right so it never occludes points; figure width padded per
    # column to leave room for it.
    color_bys = ["diff_batch", "stage", "tech_batch"]
    fig, axes = plt.subplots(len(runs), len(color_bys), figsize=(6.5 * len(color_bys), 4.5 * len(runs)))
    if len(runs) == 1:
        axes = axes[None, :]
    for i, (run_name, umap_key) in enumerate(runs.items()):
        coords = adata.obsm[umap_key]
        for j, color_by in enumerate(color_bys):
            ax = axes[i, j]
            cats = pd.Categorical(adata.obs[color_by].astype(str))
            for cat in cats.categories:
                m = cats == cat
                ax.scatter(coords[m, 0], coords[m, 1], s=2, alpha=0.5, label=str(cat))
            ax.set_title(f"{run_name}: {color_by}")
            ax.legend(markerscale=6, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left", framealpha=0.9)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--normalized-rna-h5ad", required=True, help="5.3_normalize_merged.py output")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--theta-diff-batch", type=float, default=3.0, help="Strategy A: diff_batch theta (higher = harder correction)")
    ap.add_argument("--theta-tech-batch", type=float, default=0.5, help="Strategy A: tech_batch theta (lower = light correction)")
    ap.add_argument("--random-state", type=int, default=1234)
    ap.add_argument("--replot-only", action="store_true",
                     help="Skip Harmony/UMAP entirely -- reload this outdir's existing "
                          "h5ad/merged_rna_harmony_integrated.h5ad (from a prior run) and "
                          "just regenerate the comparison plot, e.g. after a plotting fix.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)

    if args.replot_only:
        out_h5ad = outdir / "h5ad" / "merged_rna_harmony_integrated.h5ad"
        adata = sc.read_h5ad(out_h5ad)
        print(f"[replot-only] loaded: {adata.shape}", flush=True)
        umap_keys = {
            k.removeprefix("X_umap_"): k for k in adata.obsm.keys() if k.startswith("X_umap_")
        }
        if not umap_keys:
            raise ValueError(f"No X_umap_* obsm keys found in {out_h5ad} -- was this written by a completed run?")
        print(f"[replot-only] regenerating plot for runs: {list(umap_keys.keys())}", flush=True)
        plot_comparison(adata, umap_keys, outdir / "plots" / "harmony_strategy_comparison.png")
        print(f"Written: {outdir / 'plots' / 'harmony_strategy_comparison.png'}")
        return

    adata = sc.read_h5ad(args.normalized_rna_h5ad)
    print(f"loaded: {adata.shape}", flush=True)
    meta = adata.obs[["diff_batch", "stage", "tech_batch", "sample_id"]].copy()

    strategies = {
        "harmonyA_explicit": dict(vars_use=["diff_batch", "tech_batch"], theta=[args.theta_diff_batch, args.theta_tech_batch]),
        "harmonyB_per_sample": dict(vars_use=["sample_id"], theta=None),
    }
    paths = {"log": "X_pca_log", "sct": "X_pca_sct"}
    if "X_pca_sct_persample" in adata.obsm:
        paths["sct_persample"] = "X_pca_sct_persample"
    print(f"paths found in this h5ad: {list(paths.keys())}", flush=True)

    umap_keys = {}
    for path_name, pca_key in paths.items():
        pca_emb = adata.obsm[pca_key]
        for strat_name, strat in strategies.items():
            run_name = f"{path_name}_{strat_name}"
            print(f"running {run_name}: vars_use={strat['vars_use']} theta={strat['theta']}", flush=True)
            corrected = harmonize(pca_emb, meta, strat["vars_use"], strat["theta"], args.random_state)
            adata.obsm[f"X_harmony_{run_name}"] = corrected
            umap = umap_from_pca(adata.obs_names, corrected, key=run_name, random_state=args.random_state)
            adata.obsm[f"X_umap_{run_name}"] = umap
            umap_keys[run_name] = f"X_umap_{run_name}"
            print(f"  {run_name} done", flush=True)

    out_h5ad = outdir / "h5ad" / "merged_rna_harmony_integrated.h5ad"
    adata.write_h5ad(out_h5ad)

    plot_comparison(adata, umap_keys, outdir / "plots" / "harmony_strategy_comparison.png")

    print(f"\nWritten: {out_h5ad}")
    print(f"Written: {outdir / 'plots' / 'harmony_strategy_comparison.png'}")
    print(f"4 runs complete: {list(umap_keys.keys())}")


if __name__ == "__main__":
    main()
