"""Step 7.2c: scVI on the merged D15 + D22 endocrine cells, batch = stage (after 7.2a and 7.2b).

Decisions (user, 2026-11-19): batch correction on stage (D15 vs D22) only; genotype not corrected. scVI models raw
counts, so it does not use SCT values: it gets the SoupX counts of the same 3,000 genes as SCT B (7.2b, MT- genes
removed), batch_key = stage, S - G2M as a continuous covariate (the cell-cycle regression of SCT B). Model as in the
pipeline's step 5.5: 2 layers, negative binomial; 30 latent dimensions (the 30 dimensions used for every method).
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/emb_scvi.csv.gz (cells x 30; reruns overwrite)
       and scvi_history.csv (training / validation ELBO per epoch).
Run on a GPU (7.2c_integrate_scvi.sb).
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch

OUT = Path(os.environ["OUT_DIR"]) / "day_analysis" / "integration_D15_D22_endocrine"
SEED, N_LATENT = 1234, 30


def main():
    scvi.settings.seed = SEED
    a = ad.read_h5ad(OUT / "merged_counts.h5ad")
    vf = pd.read_csv(OUT / "SCT_B_variable_features.csv").gene.tolist()
    assert len(vf) == 3000 and set(vf) <= set(a.var_names)
    a = a[:, vf].copy()
    a.layers["counts"] = a.X.copy()
    scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key="stage", continuous_covariate_keys=["S_minus_G2M"])
    model = scvi.model.SCVI(a, n_layers=2, n_latent=N_LATENT, gene_likelihood="nb")
    print(model, "cuda:", torch.cuda.is_available())
    model.train(accelerator="gpu", devices=1)
    z = pd.DataFrame(model.get_latent_representation(), index=a.obs_names,
                     columns=[f"scVI_{i + 1}" for i in range(N_LATENT)])
    z.to_csv(OUT / "emb_scvi.csv.gz")
    h = model.history
    pd.concat([h[k] for k in ["elbo_train", "elbo_validation"] if k in h], axis=1).to_csv(OUT / "scvi_history.csv")
    print("epochs:", len(h["elbo_train"]), "; wrote", OUT / "emb_scvi.csv.gz")


if __name__ == "__main__":
    main()
