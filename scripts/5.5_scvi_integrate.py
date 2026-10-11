#!/usr/bin/env python
"""
scVI integration, two strategies, same comparison framing as
5.4_harmony_integrate.py -- but a genuinely different mechanism, not a
drop-in replacement:

  - Harmony corrects an EXISTING PCA embedding (5.3's log/SCT paths).
  - scVI models RAW COUNTS end-to-end (its own NB/ZINB likelihood learns
    its own latent space) -- so this script runs on 5.2's merged raw
    output directly, bypassing 5.3 entirely. There is no "log vs SCT"
    axis here, just the 2 strategies -> 2 runs, not 4.

Strategy A ("explicit"): batch_key="diff_batch" (the dominant real batch
effect, baked into scVI's encoder/decoder/library-size handling),
categorical_covariate_keys=["tech_batch"] (a lighter conditioning
covariate, confirmed near-negligible noise). IMPORTANT: this is NOT the
same correction mechanism as Harmony's theta -- there is no scalar
"strength" knob for categorical_covariate_keys; it's architecturally a
conditioning input, not a tunable-strength penalty. Don't read scVI's
strategy A as "the scVI equivalent of theta=3/0.5" -- it isn't one.
Strategy B ("per_sample"): batch_key="sample_id", no covariates -- matches
Harmony strategy B's "each sample is its own batch" framing.

Neither strategy ever uses "stage" as batch_key or covariate -- that's
the biological signal to preserve, not batch noise to remove. Same
IMPORTANT caveat as 5.4: most stages here are covered by only ONE
diff_batch (see 5.1_coverage_table.py) -- no strategy can manufacture
cross-diff_batch evidence where none exists.

HVG selection uses flavor="seurat_v3" (needs raw counts, not log-data)
with batch_key="sample_id" -- selects per-sample then combines, avoiding
the same batch-driven-gene-selection pitfall flagged for 5.3's log path.

GPU required (aa100 partition, gpu-normal QOS) -- confirmed working
CUDA-compiled torch + scvi-tools 1.2.2.post1 in the single.cell env
earlier this session.

Usage:
  python 5.5_scvi_integrate.py --merged-rna-h5ad merged_rna_raw.h5ad \
      --outdir out/ [--n-top-genes 2000] [--n-latent 30] [--max-epochs N]
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


def select_hvgs(adata, n_top_genes: int) -> ad.AnnData:
    """seurat_v3 flavor needs raw counts (not log-transformed) -- matches
    what merged_rna_raw.h5ad actually holds (5.2's docstring: raw
    SoupX-corrected counts). batch_key="sample_id" selects per-sample
    then combines, same reasoning as flagged for 5.3's log path (avoids
    picking sample-dominant genes from a pooled, compositionally-skewed
    variance calculation)."""
    counts_source = "soupx_counts layer" if "soupx_counts" in adata.layers else "X"
    print(f"  [select_hvgs] using {counts_source} as raw counts input", flush=True)
    a = adata.copy()
    if "soupx_counts" in a.layers:
        a.X = a.layers["soupx_counts"]
    sc.pp.highly_variable_genes(
        a, flavor="seurat_v3", n_top_genes=n_top_genes, batch_key="sample_id", subset=True,
    )
    print(f"  [select_hvgs] {a.n_vars} HVGs selected (batch_key=sample_id)", flush=True)
    return a


def umap_from_latent(obs_names, latent: np.ndarray, key, random_state):
    """Same lightweight pattern as 5.4_harmony_integrate.py's
    umap_from_pca() -- neighbors/umap only touch obsm/obsp, never .X, so
    build a minimal 0-gene AnnData rather than copying anything large."""
    a = ad.AnnData(X=np.empty((latent.shape[0], 0), dtype=np.float32), obs=pd.DataFrame(index=obs_names))
    a.obsm["X_scvi"] = latent
    sc.pp.neighbors(a, use_rep="X_scvi", random_state=random_state, key_added=key)
    sc.tl.umap(a, neighbors_key=key, random_state=random_state)
    return a.obsm["X_umap"]


def plot_comparison(adata, runs: dict, outpath: Path):
    """Legend on every panel from the start -- 5.4_harmony_integrate.py's
    first version only legended the last column, silently unreadable for
    the other two. Not repeating that here."""
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


def run_scvi(adata_hvg, batch_key, categorical_covariate_keys, n_latent, max_epochs, random_state):
    import scvi

    scvi.settings.seed = random_state
    a = adata_hvg.copy()
    scvi.model.SCVI.setup_anndata(
        a,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
        layer=None,  # already raw counts (see select_hvgs)
    )
    model = scvi.model.SCVI(a, n_layers=2, n_latent=n_latent, gene_likelihood="nb")
    model.train(max_epochs=max_epochs, accelerator="gpu", devices=1)
    return model.get_latent_representation()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-rna-h5ad", required=True, help="5.2_merge_rna.py output: merged_rna_raw.h5ad")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-top-genes", type=int, default=2000, help="matches 4.1/5.3's log-path HVG count")
    ap.add_argument("--n-latent", type=int, default=30, help="matches the pipeline's n_pcs convention")
    ap.add_argument("--max-epochs", type=int, default=None, help="default: scvi-tools' own heuristic")
    ap.add_argument("--random-state", type=int, default=1234)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.merged_rna_h5ad)
    print(f"loaded: {adata.shape}", flush=True)
    print(f"layers: {list(adata.layers.keys())}", flush=True)

    adata_hvg = select_hvgs(adata, args.n_top_genes)

    strategies = {
        "scviA_explicit": dict(batch_key="diff_batch", categorical_covariate_keys=["tech_batch"]),
        "scviB_per_sample": dict(batch_key="sample_id", categorical_covariate_keys=None),
    }

    umap_keys = {}
    for strat_name, strat in strategies.items():
        print(f"running {strat_name}: batch_key={strat['batch_key']} categorical_covariate_keys={strat['categorical_covariate_keys']}", flush=True)
        latent = run_scvi(
            adata_hvg, strat["batch_key"], strat["categorical_covariate_keys"],
            args.n_latent, args.max_epochs, args.random_state,
        )
        adata.obsm[f"X_scvi_{strat_name}"] = latent
        umap = umap_from_latent(adata.obs_names, latent, key=strat_name, random_state=args.random_state)
        adata.obsm[f"X_umap_{strat_name}"] = umap
        umap_keys[strat_name] = f"X_umap_{strat_name}"
        print(f"  {strat_name} done, latent shape {latent.shape}", flush=True)

    out_h5ad = outdir / "h5ad" / "merged_rna_scvi_integrated.h5ad"
    adata.write_h5ad(out_h5ad)

    plot_comparison(adata, umap_keys, outdir / "plots" / "scvi_strategy_comparison.png")

    print(f"\nWritten: {out_h5ad}")
    print(f"Written: {outdir / 'plots' / 'scvi_strategy_comparison.png'}")
    print(f"2 runs complete: {list(umap_keys.keys())}")


if __name__ == "__main__":
    main()
