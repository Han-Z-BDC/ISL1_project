"""Step 5.4c: Harmony with batch variables chosen for this study's design (replaces 5.4's A/B).

Design (user, 2026-11-08): two differentiation batches -- Diff34 = D40 only, Diff37 = D15 + D22.
rep1/rep2 are biological replicates within a differentiation, NOT paired across days or
genotypes, so tech_batch (rep1/rep2) means nothing across samples and is never a batch variable.
Three ways to handle batch, run on every normalization path in the 5.3 output:
  - none:        no correction -- 5.3's PCA as is ({path}_unintegrated in 5.4b)
  - diffbatch:   vars_use=["diff_batch"], theta=HARMONY_THETA_DIFF_BATCH (config.sh). Confounded
                 with day (D40 vs D15/D22), so it can also remove real D40-vs-early differences.
  - persample:   vars_use=["sample_id"], harmonypy's default theta -- each of the 12 samples is
                 its own batch. Also pulls ctr/KO and days together.
"stage" (day-genotype) is never a correction variable. Reuses 5.4_harmony_integrate.py's
harmonize() / umap_from_pca() / plot_comparison() so the method matches the pipeline.

Writes ${OUT_DIR}/integration/h5ad/merged_rna_harmony_design.h5ad (obs + obsm only:
X_harmony_{path}_{strategy}, X_umap_{path}_{strategy}) and plots/harmony_design_comparison.png.
Run: sbatch analysis/scripts/5.4c_harmony_by_design.sb
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
    ap.add_argument("--normalized-rna-h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--theta-diff-batch", type=float, default=3.0)
    ap.add_argument("--random-state", type=int, default=1234)
    args = ap.parse_args()
    h = load_pipeline_module("5.4_harmony_integrate.py")
    out = Path(args.outdir)
    (out / "h5ad").mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.normalized_rna_h5ad)
    print(f"loaded: {adata.shape}", flush=True)
    meta = adata.obs[["diff_batch", "stage", "tech_batch", "sample_id"]].copy()
    strategies = {
        "harmony_diffbatch": dict(vars_use=["diff_batch"], theta=[args.theta_diff_batch]),
        "harmony_persample": dict(vars_use=["sample_id"], theta=None),
    }
    paths = {p: f"X_pca_{p}" for p in ("log", "sct", "sct_persample") if f"X_pca_{p}" in adata.obsm}

    res = ad.AnnData(obs=adata.obs.copy())
    umap_keys = {}
    for path, pca_key in paths.items():
        for name, st in strategies.items():
            run = f"{path}_{name}"
            print(f"running {run}: vars_use={st['vars_use']} theta={st['theta']}", flush=True)
            z = h.harmonize(adata.obsm[pca_key], meta, st["vars_use"], st["theta"], args.random_state)
            res.obsm[f"X_harmony_{run}"] = z
            res.obsm[f"X_umap_{run}"] = h.umap_from_pca(adata.obs_names, z, key=run, random_state=args.random_state)
            umap_keys[run] = f"X_umap_{run}"
    res.write_h5ad(out / "h5ad" / "merged_rna_harmony_design.h5ad")
    h.plot_comparison(res, umap_keys, out / "plots" / "harmony_design_comparison.png")
    print("written:", out / "h5ad" / "merged_rna_harmony_design.h5ad", flush=True)


if __name__ == "__main__":
    main()
