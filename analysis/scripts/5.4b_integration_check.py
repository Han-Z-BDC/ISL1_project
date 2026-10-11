"""Step 5.4b: is integration needed at all? Unintegrated vs Harmony (and scVI, if present).

All 4 samples share one OCM chip channel and library, so there is no technical batch;
replicates (rep1/rep2) are separate differentiation wells. This compares, per
normalization path, the UNINTEGRATED embedding (pipeline 5.3's PCA, no correction) with
Harmony strategy A (light replicate correction; genotype kept in the never-corrected
"stage" field) and, when 5.5 has run, scVI strategy A. Harmony/scVI strategy B (per
sample) also pull ctr and KO together, so they are shown for reference only.

Metrics per embedding (k=30 nearest neighbours in the embedding itself, not the UMAP):
  - rep_mixing: within each genotype, fraction of a cell's same-genotype neighbours that
    come from the OTHER replicate, divided by that replicate's share of the genotype
    (1.0 = replicates fully mixed; ~0 = replicates separate).
  - genotype_purity: fraction of neighbours with the same genotype (ctr vs KO).
    Integration should not lower this much -- that would mean real biology being removed.

Writes to ${OUT_DIR}/integration_check/: integration_check.csv, embeddings_by_group.png,
markers_unintegrated_log.png, and h5ad/integration_check.h5ad (UMAPs added).

Run through Slurm (UMAP/kNN on 12.7k cells): sbatch analysis/scripts/5.4b_integration_check.sb
"""
import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors

MARKERS = ["INS", "GCG", "SST", "TPH1", "FEV", "CHGA", "NKX6-1", "ISL1", "ARX", "MAFB",
           "SOX9", "KRT19", "VIM", "COL1A1", "MKI67", "GHRL"]
K = 30


def knn_metrics(emb, genotype, rep):
    nn = NearestNeighbors(n_neighbors=K + 1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)[:, 1:]
    g_nb, r_nb = genotype[idx], rep[idx]
    purity = (g_nb == genotype[:, None]).mean()
    mix = []
    for g in np.unique(genotype):
        m = genotype == g
        share_other = {r: (rep[m] != r).mean() for r in np.unique(rep[m])}
        same_g = g_nb[m] == g
        other_rep = (r_nb[m] != rep[m][:, None]) & same_g
        frac = other_rep.sum(1) / np.maximum(same_g.sum(1), 1)
        expected = np.array([share_other[r] for r in rep[m]])
        mix.append((frac / expected).mean())
    return float(np.mean(mix)), float(purity)


def umap_of(emb, obs_names, seed=0):
    a = ad.AnnData(obs=pd.DataFrame(index=obs_names))
    a.obsm["E"] = emb
    sc.pp.neighbors(a, use_rep="E", n_neighbors=K, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    return a.obsm["X_umap"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--integration-dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    idir, out = Path(args.integration_dir), Path(args.outdir)
    (out / "h5ad").mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(idir / "merged_rna_normalized.h5ad")
    h = ad.read_h5ad(idir / "h5ad" / "merged_rna_harmony_integrated.h5ad", backed="r")
    assert (h.obs_names == a.obs_names).all()
    genotype = a.obs["stage"].astype(str).to_numpy()      # D40-ctr / D40-KO
    rep = a.obs["tech_batch"].astype(str).to_numpy()      # rep1 / rep2

    # embedding name -> (latent for metrics, existing UMAP or None)
    embs = {}
    for path, pca in [("log", "X_pca_log"), ("sct", "X_pca_sct"), ("sct_persample", "X_pca_sct_persample")]:
        embs[f"{path}_unintegrated"] = (a.obsm[pca], None)
        for strat in ("harmonyA_explicit", "harmonyB_per_sample"):
            run = f"{path}_{strat}"
            embs[run] = (np.asarray(h.obsm[f"X_harmony_{run}"]), np.asarray(h.obsm[f"X_umap_{run}"]))
    scvi_h5 = idir / "h5ad" / "merged_rna_scvi_integrated.h5ad"
    if scvi_h5.exists():
        s = ad.read_h5ad(scvi_h5, backed="r")
        s_idx = pd.Index(s.obs_names).get_indexer(a.obs_names)
        assert (s_idx >= 0).all()
        for strat in ("scviA_explicit", "scviB_per_sample"):
            if f"X_scvi_{strat}" in s.obsm:
                embs[strat] = (np.asarray(s.obsm[f"X_scvi_{strat}"])[s_idx], np.asarray(s.obsm[f"X_umap_{strat}"])[s_idx])

    rows, umaps = [], {}
    for name, (lat, um) in embs.items():
        mix, pur = knn_metrics(lat, genotype, rep)
        rows.append({"embedding": name, "rep_mixing": round(mix, 3), "genotype_purity": round(pur, 3)})
        umaps[name] = um if um is not None else umap_of(lat, a.obs_names)
        a.obsm[f"X_umap_{name}"] = umaps[name]
        print(f"{name}: rep_mixing={mix:.3f} genotype_purity={pur:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(out / "integration_check.csv", index=False)

    # embeddings coloured by genotype and replicate
    names = list(umaps)
    fig, axes = plt.subplots(len(names), 2, figsize=(11, 4.6 * len(names)))
    pal = {"D40-ctr": "#3b6fb6", "D40-KO": "#d0573f", "rep1": "#4c9a6a", "rep2": "#9b6bc0"}
    rng = np.random.default_rng(0)
    order = rng.permutation(a.n_obs)  # random draw order so no group hides another
    for i, name in enumerate(names):
        u = umaps[name][order]
        for j, (lab, vals) in enumerate((("genotype", genotype[order]), ("replicate", rep[order]))):
            ax = axes[i, j]
            ax.scatter(u[:, 0], u[:, 1], s=1.5, c=[pal[v] for v in vals], linewidths=0)
            r = next(x for x in rows if x["embedding"] == name)
            ax.set_title(f"{name}: {lab}\nrep_mixing={r['rep_mixing']}  genotype_purity={r['genotype_purity']}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
            for k in np.unique(vals):
                ax.scatter([], [], c=pal[k], s=20, label=k)
            ax.legend(fontsize=8, loc="best", markerscale=1.5)
    fig.tight_layout()
    fig.savefig(out / "embeddings_by_group.png", dpi=90)
    plt.close(fig)

    # markers on the unintegrated log UMAP -- high-expressing cells drawn on top
    lay = "soupx_counts" if "soupx_counts" in a.layers else "counts"
    b = ad.AnnData(X=a.layers[lay].copy(), obs=a.obs[[]], var=a.var[[]])
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    u = umaps["log_unintegrated"]
    genes = [g for g in MARKERS if g in b.var_names]
    ncol = 4
    fig, axes = plt.subplots(int(np.ceil(len(genes) / ncol)), ncol, figsize=(4.2 * ncol, 3.9 * np.ceil(len(genes) / ncol)))
    for ax, g in zip(axes.ravel(), genes):
        x = b[:, g].X
        x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
        o = np.argsort(x)
        im = ax.scatter(u[o, 0], u[o, 1], c=x[o], s=1.5, cmap="viridis", linewidths=0)
        ax.set_title(g); ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.04)
    for ax in axes.ravel()[len(genes):]:
        ax.axis("off")
    fig.suptitle("log path, unintegrated UMAP: log-normalized SoupX counts", y=1.0)
    fig.tight_layout()
    fig.savefig(out / "markers_unintegrated_log.png", dpi=90)
    plt.close(fig)

    keep = ad.AnnData(obs=a.obs, obsm={k: v for k, v in a.obsm.items() if k.startswith("X_umap_")})
    keep.write_h5ad(out / "h5ad" / "integration_check.h5ad")
    print("wrote", out)


if __name__ == "__main__":
    main()
