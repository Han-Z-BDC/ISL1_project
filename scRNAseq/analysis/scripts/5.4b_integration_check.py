"""Step 5.4b: is integration needed at all? Unintegrated vs Harmony (and scVI, if present).

Each day (D15, D22, D40) is one OCM chip channel of 4 samples (ctr/KO x rep1/rep2);
replicates are separate differentiation wells. Day is confounded with channel and with
differentiation (Diff37 = D15 + D22, Diff34 = D40), so integration on diff_batch can remove
real day differences -- day_purity tracks that. This compares, per
normalization path, the UNINTEGRATED embedding (pipeline 5.3's PCA, no correction) with
Harmony on diff_batch only and Harmony per sample (5.4c), and the three scVI models
(no batch / diff_batch / per sample, 5.5c). Replicates are not paired across samples, so
no embedding corrects on rep1/rep2 (2026-11-08 design). Per-sample correction also pulls
ctr/KO and days together.

Metrics per embedding (k=30 nearest neighbours in the embedding itself, not the UMAP):
  - rep_mixing: within each genotype, fraction of a cell's same-genotype neighbours that
    come from the OTHER replicate, divided by that replicate's share of the genotype
    (1.0 = replicates fully mixed; ~0 = replicates separate).
  - genotype_purity: among neighbours from the same day, fraction with the same genotype
    (ctr vs KO). Integration should not lower this much -- that would mean real biology
    being removed.
  - day_purity: fraction of neighbours from the same day. Should stay high for an
    embedding used to compare days.

Writes to ${OUT_DIR}/integration_check/: integration_check.csv, embeddings_by_group.png,
embeddings/{embedding}.png (merged / ctr only / KO only / replicate / day, shared axes),
markers_unintegrated_log.png, and h5ad/integration_check.h5ad (UMAPs added).

Run through Slurm (UMAP/kNN on ~60k cells): sbatch analysis/scripts/5.4b_integration_check.sb
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


def knn_metrics(emb, stage, rep):
    """stage = <day>-<genotype>; replicates are compared within each stage."""
    day = np.array([x.split("-")[0] for x in stage])
    geno = np.array([x.split("-")[1] for x in stage])
    nn = NearestNeighbors(n_neighbors=K + 1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)[:, 1:]
    g_nb, r_nb = stage[idx], rep[idx]
    same_day = day[idx] == day[:, None]
    day_purity = same_day.mean()
    geno_purity = ((geno[idx] == geno[:, None]) & same_day).sum() / max(same_day.sum(), 1)
    mix = []
    for g in np.unique(stage):
        m = stage == g
        share_other = {r: (rep[m] != r).mean() for r in np.unique(rep[m])}
        same_g = g_nb[m] == g
        other_rep = (r_nb[m] != rep[m][:, None]) & same_g
        frac = other_rep.sum(1) / np.maximum(same_g.sum(1), 1)
        expected = np.array([share_other[r] for r in rep[m]])
        mix.append((frac / expected).mean())
    return float(np.mean(mix)), float(geno_purity), float(day_purity)


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
    h = ad.read_h5ad(idir / "h5ad" / "merged_rna_harmony_design.h5ad")
    assert (h.obs_names == a.obs_names).all()
    stage = a.obs["stage"].astype(str).to_numpy()         # D15-ctr ... D40-KO
    genotype = np.array([x.split("-")[1] for x in stage])  # ctr / KO
    day = np.array([x.split("-")[0] for x in stage])       # D15 / D22 / D40
    rep = a.obs["tech_batch"].astype(str).to_numpy()      # rep1 / rep2

    # embedding name -> (latent for metrics, existing UMAP or None)
    embs = {}
    for path, pca in [("log", "X_pca_log"), ("sct", "X_pca_sct"), ("sct_persample", "X_pca_sct_persample")]:
        embs[f"{path}_unintegrated"] = (a.obsm[pca], None)
        for strat in ("harmony_diffbatch", "harmony_persample"):
            run = f"{path}_{strat}"
            embs[run] = (np.asarray(h.obsm[f"X_harmony_{run}"]), np.asarray(h.obsm[f"X_umap_{run}"]))
    scvi_h5 = idir / "h5ad" / "merged_rna_scvi_design.h5ad"
    if scvi_h5.exists():
        s = ad.read_h5ad(scvi_h5)
        s_idx = pd.Index(s.obs_names).get_indexer(a.obs_names)
        assert (s_idx >= 0).all()
        for strat in ("scvi_nobatch", "scvi_diffbatch", "scvi_persample"):
            if f"X_scvi_{strat}" in s.obsm:
                embs[strat] = (np.asarray(s.obsm[f"X_scvi_{strat}"])[s_idx], np.asarray(s.obsm[f"X_umap_{strat}"])[s_idx])

    rows, umaps = [], {}
    for name, (lat, um) in embs.items():
        mix, pur, dpur = knn_metrics(lat, stage, rep)
        rows.append({"embedding": name, "rep_mixing": round(mix, 3), "genotype_purity": round(pur, 3),
                     "day_purity": round(dpur, 3)})
        umaps[name] = um if um is not None else umap_of(lat, a.obs_names)
        a.obsm[f"X_umap_{name}"] = umaps[name]
        print(f"{name}: rep_mixing={mix:.3f} genotype_purity={pur:.3f} day_purity={dpur:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(out / "integration_check.csv", index=False)

    # embeddings coloured by genotype, replicate and day
    names = list(umaps)
    fig, axes = plt.subplots(len(names), 3, figsize=(16.5, 4.6 * len(names)))
    pal = {"ctr": "#3b6fb6", "KO": "#d0573f", "rep1": "#4c9a6a", "rep2": "#9b6bc0",
           "D15": "#e3a72f", "D22": "#2a9d8f", "D40": "#7a4fa3"}
    rng = np.random.default_rng(0)
    order = rng.permutation(a.n_obs)  # random draw order so no group hides another
    for i, name in enumerate(names):
        u = umaps[name][order]
        for j, (lab, vals) in enumerate((("genotype", genotype[order]), ("replicate", rep[order]), ("day", day[order]))):
            ax = axes[i, j]
            ax.scatter(u[:, 0], u[:, 1], s=1.5, c=[pal[v] for v in vals], linewidths=0)
            r = next(x for x in rows if x["embedding"] == name)
            ax.set_title(f"{name}: {lab}\nrep_mixing={r['rep_mixing']}  genotype_purity={r['genotype_purity']}  "
                         f"day_purity={r['day_purity']}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
            for k in np.unique(vals):
                ax.scatter([], [], c=pal[k], s=20, label=k)
            ax.legend(fontsize=8, loc="best", markerscale=1.5)
    fig.tight_layout()
    fig.savefig(out / "embeddings_by_group.png", dpi=90)
    plt.close(fig)

    # one row per embedding: merged / ctr only / KO only / replicate / day, same axes in every panel
    (out / "embeddings").mkdir(exist_ok=True)
    grey = "#d9d9d9"
    for name in names:
        u = umaps[name]
        r = next(x for x in rows if x["embedding"] == name)
        fig, axes = plt.subplots(1, 5, figsize=(23.5, 4.6), sharex=True, sharey=True)
        ax = axes[0]
        ax.scatter(u[order, 0], u[order, 1], s=1.5, c=[pal[v] for v in genotype[order]], linewidths=0)
        ax.set_title("merged (ctr + KO)")
        for k in ("ctr", "KO"):
            ax.scatter([], [], c=pal[k], s=20, label=k)
        ax.legend(fontsize=8, markerscale=1.5)
        for ax, g in zip(axes[1:3], ("ctr", "KO")):
            m = genotype == g
            ax.scatter(u[~m, 0], u[~m, 1], s=1.5, c=grey, linewidths=0)
            ax.scatter(u[m, 0], u[m, 1], s=1.5, c=pal[g], linewidths=0)
            ax.set_title(f"{g} only (n={m.sum():,}; other genotype in grey)")
        ax = axes[3]
        ax.scatter(u[order, 0], u[order, 1], s=1.5, c=[pal[v] for v in rep[order]], linewidths=0)
        for k in ("rep1", "rep2"):
            ax.scatter([], [], c=pal[k], s=20, label=k)
        ax.legend(fontsize=8, markerscale=1.5)
        ax.set_title("replicate")
        ax = axes[4]
        ax.scatter(u[order, 0], u[order, 1], s=1.5, c=[pal[v] for v in day[order]], linewidths=0)
        for k in sorted(np.unique(day)):
            ax.scatter([], [], c=pal[k], s=20, label=k)
        ax.legend(fontsize=8, markerscale=1.5)
        ax.set_title("day")
        for ax in axes:
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{name}   (rep_mixing={r['rep_mixing']}, genotype_purity={r['genotype_purity']}, "
                     f"day_purity={r['day_purity']})", fontsize=11)
        fig.tight_layout()
        fig.savefig(out / "embeddings" / f"{name}.png", dpi=80)
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
