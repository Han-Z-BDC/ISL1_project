"""Step 5.5c: scVI with batch handling chosen for this study's design (replaces 5.5's A/B).

Same design reasoning as 5.4c_harmony_by_design.py: tech_batch (rep1/rep2) is not paired across
samples, so it is never a batch variable or covariate. Three scVI models on the same HVGs:
  - scvi_nobatch:   batch_key=None (no correction; scVI's own latent space)
  - scvi_diffbatch: batch_key="diff_batch" (confounded with day: D40 vs D15/D22)
  - scvi_persample: batch_key="sample_id" (each of the 12 samples is its own batch)
"stage" is never a batch key. Reuses 5.5_scvi_integrate.py's select_hvgs() (2,000 seurat_v3
HVGs, batch_key=sample_id), run_scvi() (2 layers, 30 latent, NB) and plotting.

Writes ${OUT_DIR}/integration/h5ad/merged_rna_scvi_design.h5ad (obs + obsm only:
X_scvi_{strategy}, X_umap_{strategy}) and plots/scvi_design_comparison.png.
Run: sbatch analysis/scripts/5.5c_scvi_by_design.sb (GPU)
"""
import argparse
import importlib.util
from pathlib import Path

import anndata as ad
import scanpy as sc

REPO = Path(__file__).resolve().parents[2]


def load_pipeline_module(fname):
    spec = importlib.util.spec_from_file_location(fname.replace(".", "_"), REPO / "pipeline" / "scripts" / fname)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged-rna-h5ad", required=True, help="5.2 output: merged_rna_raw.h5ad")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n-top-genes", type=int, default=2000)
    ap.add_argument("--n-latent", type=int, default=30)
    ap.add_argument("--max-epochs", type=int, default=None)
    ap.add_argument("--random-state", type=int, default=1234)
    args = ap.parse_args()
    s = load_pipeline_module("5.5_scvi_integrate.py")
    out = Path(args.outdir)
    (out / "h5ad").mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.merged_rna_h5ad)
    print(f"loaded: {adata.shape}", flush=True)
    hvg = s.select_hvgs(adata, args.n_top_genes)
    strategies = {"scvi_nobatch": None, "scvi_diffbatch": "diff_batch", "scvi_persample": "sample_id"}

    res = ad.AnnData(obs=adata.obs.copy())
    umap_keys = {}
    for name, batch_key in strategies.items():
        print(f"running {name}: batch_key={batch_key}", flush=True)
        z = s.run_scvi(hvg, batch_key, None, args.n_latent, args.max_epochs, args.random_state)
        res.obsm[f"X_scvi_{name}"] = z
        res.obsm[f"X_umap_{name}"] = s.umap_from_latent(adata.obs_names, z, key=name, random_state=args.random_state)
        umap_keys[name] = f"X_umap_{name}"
    res.write_h5ad(out / "h5ad" / "merged_rna_scvi_design.h5ad")
    s.plot_comparison(res, umap_keys, out / "plots" / "scvi_design_comparison.png")
    print("written:", out / "h5ad" / "merged_rna_scvi_design.h5ad", flush=True)


if __name__ == "__main__":
    main()
