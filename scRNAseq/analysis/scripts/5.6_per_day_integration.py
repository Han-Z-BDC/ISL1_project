"""Step 5.6: integration within each time point (D15, D22, D40 analysed separately).

Each day is one differentiation batch and one OCM channel, so within a day the only batch-like
factor is the sample (ctr/KO x rep1/rep2, replicates not paired). For each day, the 4 samples
are taken from the 12-sample merged object (5.2) and re-processed on their own:
  - none:      pipeline 5.3 log path on that day only (normalize_total -> log1p -> 2,000 HVGs
               -> 50-PC PCA), no correction
  - harmony:   Harmony on that PCA, vars_use=["sample_id"], harmonypy default theta
  - scvi:      scVI, batch_key="sample_id" (2,000 seurat_v3 HVGs per sample, 2 layers, 30 latent)
UMAPs: same neighbours/UMAP call as pipeline 5.4/5.5 (umap_from_pca / umap_from_latent).
Metrics (k=30 kNN in the embedding, from 5.4b): rep_mixing (replicates of the same genotype mix;
1 = fully) and genotype_purity (neighbours with the same genotype).

Writes ${OUT_DIR}/integration_per_day/: {day}.h5ad (obs + embeddings + UMAPs),
umap_{day}.png (rows = method; columns = sample, genotype, ctr only, KO only),
per_day_metrics.csv, per_day_integration.html.
Run: sbatch analysis/scripts/5.6_per_day_integration.sb (GPU, for scVI)
"""
import base64
import importlib.util
import os
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = Path(os.environ["OUT_DIR"])
RES = OUT_DIR / "integration_per_day"
SEED = 1234
METHODS = {"none": "no integration (log PCA)", "harmony": "Harmony, per sample", "scvi": "scVI, per sample"}
GENO_PAL = {"ctr": "#3b6fb6", "KO": "#d0573f"}
SAMPLE_PAL = {"ctr_rep1": "#2f5d9e", "ctr_rep2": "#8fb3e3", "KO_rep1": "#b23b2a", "KO_rep2": "#f0a38f"}


def load(rel, name):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def scatter_groups(ax, u, labels, pal, order, title):
    for k in pal:
        ax.scatter([], [], c=pal[k], s=25, label=k)
    ax.scatter(u[order, 0], u[order, 1], s=1.5, c=[pal[x] for x in labels[order]], linewidths=0)
    ax.legend(fontsize=8, markerscale=1.3, frameon=False, loc="best")
    ax.set_title(title, fontsize=10)


def main():
    RES.mkdir(parents=True, exist_ok=True)
    norm = load("pipeline/scripts/5.3_normalize_merged.py", "norm53")
    harm = load("pipeline/scripts/5.4_harmony_integrate.py", "harm54")
    scvi_m = load("pipeline/scripts/5.5_scvi_integrate.py", "scvi55")
    chk = load("analysis/scripts/5.4b_integration_check.py", "chk54b")
    cfg = yaml.safe_load((REPO / "pipeline" / "configs" / "normalize_rna.yaml").read_text())

    full = sc.read_h5ad(OUT_DIR / "integration" / "merged_rna_raw.h5ad")
    full.obs["day"] = full.obs["stage"].astype(str).str.split("-").str[0]
    full.obs["genotype"] = full.obs["stage"].astype(str).str.split("-").str[1]
    full.obs["sample_short"] = full.obs["genotype"] + "_" + full.obs["tech_batch"].astype(str)
    days = sorted(full.obs["day"].unique(), key=lambda d: int(d[1:]))

    rows, figs = [], {}
    for day in days:
        a = full[full.obs["day"] == day].copy()
        print(f"== {day}: {a.n_obs} cells, samples {sorted(a.obs['sample_id'].unique())}", flush=True)
        a.obs["sample_id"] = a.obs["sample_id"].astype(str).astype("category")
        meta = a.obs[["sample_id"]].copy()

        pca = norm.run_log_path(a, cfg)
        emb = {"none": pca,
               "harmony": harm.harmonize(pca, meta, ["sample_id"], None, SEED)}
        hvg = scvi_m.select_hvgs(a, 2000)
        emb["scvi"] = scvi_m.run_scvi(hvg, "sample_id", None, 30, None, SEED)

        res = ad.AnnData(obs=a.obs[["sample_id", "stage", "tech_batch", "genotype", "sample_short",
                                    "total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy())
        stage = a.obs["stage"].astype(str).to_numpy()
        rep = a.obs["tech_batch"].astype(str).to_numpy()
        for m, z in emb.items():
            res.obsm[f"X_{m}"] = np.asarray(z)
            res.obsm[f"X_umap_{m}"] = harm.umap_from_pca(a.obs_names, np.asarray(z), key=f"{day}_{m}", random_state=SEED)
            mix, gpur, _ = chk.knn_metrics(np.asarray(z), stage, rep)
            rows.append({"day": day, "method": m, "n_cells": a.n_obs, "rep_mixing": round(mix, 3),
                         "genotype_purity": round(gpur, 3)})
            print(f"  {m}: rep_mixing={mix:.3f} genotype_purity={gpur:.3f}", flush=True)
        res.write_h5ad(RES / f"{day}.h5ad")

        geno = res.obs["genotype"].astype(str).to_numpy()
        samp = res.obs["sample_short"].astype(str).to_numpy()
        n_by = res.obs["sample_short"].value_counts()
        order = np.random.default_rng(0).permutation(res.n_obs)
        fig, axes = plt.subplots(len(METHODS), 4, figsize=(21, 4.9 * len(METHODS)))
        for i, (m, mname) in enumerate(METHODS.items()):
            u = res.obsm[f"X_umap_{m}"]
            r = next(x for x in rows if x["day"] == day and x["method"] == m)
            spal = {f"{k} (n={n_by.get(k, 0):,})": v for k, v in SAMPLE_PAL.items()}
            slab = np.array([f"{x} (n={n_by.get(x, 0):,})" for x in samp])
            scatter_groups(axes[i, 0], u, slab, spal, order, f"{mname}\nsample")
            scatter_groups(axes[i, 1], u, geno, GENO_PAL, order,
                           f"genotype  (rep_mixing={r['rep_mixing']}, genotype_purity={r['genotype_purity']})")
            for ax, g in zip(axes[i, 2:], ("ctr", "KO")):
                msk = geno == g
                ax.scatter(u[~msk, 0], u[~msk, 1], s=1.5, c="#dcdcdc", linewidths=0)
                for rr, col in (("rep1", SAMPLE_PAL[f"{g}_rep1"]), ("rep2", SAMPLE_PAL[f"{g}_rep2"])):
                    mm = msk & (res.obs["tech_batch"].astype(str).to_numpy() == rr)
                    ax.scatter(u[mm, 0], u[mm, 1], s=1.5, c=col, linewidths=0, label=f"{g} {rr}")
                ax.legend(fontsize=8, markerscale=4, frameon=False, loc="best")
                ax.set_title(f"{g} only, by replicate (other genotype grey)", fontsize=10)
            for ax in axes[i]:
                ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{day}: {res.n_obs:,} cells, 4 samples", fontsize=14, y=1.0)
        fig.tight_layout()
        png = RES / f"umap_{day}.png"
        fig.savefig(png, dpi=80, bbox_inches="tight")
        plt.close(fig)
        figs[day] = base64.b64encode(png.read_bytes()).decode()

    met = pd.DataFrame(rows)
    met.to_csv(RES / "per_day_metrics.csv", index=False)
    print(met.to_string(index=False), flush=True)

    trows = "".join(f"<tr><td>{r.day}</td><td>{METHODS[r.method]}</td><td>{r.n_cells:,}</td>"
                    f"<td>{r.rep_mixing}</td><td>{r.genotype_purity}</td></tr>" for r in met.itertuples())
    sections = "".join(f'<h2>{d}</h2><img src="data:image/png;base64,{b}" alt="{d} UMAPs">' for d, b in figs.items())
    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Per-Day Integration</title>
<style>
:root{{--bg:#fbfbf9;--fg:#1f1f1c;--mut:#6b6b66;--line:#e2e2dc}}
@media (prefers-color-scheme:dark){{:root{{--bg:#1a1a18;--fg:#ececea;--mut:#a3a39e;--line:#3a3a36}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.55 system-ui,-apple-system,sans-serif;margin:0}}
main{{max-width:1240px;margin:0 auto;padding:24px 16px 60px}}
h1{{font-size:24px;margin:0 0 4px}} h2{{font-size:19px;margin:30px 0 8px}}
.sub{{color:var(--mut);margin:0 0 16px}} table{{border-collapse:collapse;font-size:14px}}
th,td{{padding:4px 10px;border-bottom:1px solid var(--line);text-align:left}}
td:nth-child(n+3){{text-align:right;font-variant-numeric:tabular-nums}}
img{{max-width:100%;height:auto;background:#fff;border-radius:4px}} ul{{padding-left:20px}}
</style></head><body><main>
<h1>Integration within each time point</h1>
<p class=sub>D15, D22, D40 processed separately · 4 samples each (ctr/KO x rep1/rep2) · out_all/integration_per_day</p>
<table><tr><th>Day</th><th>Method</th><th>Cells</th><th>rep_mixing</th><th>genotype_purity</th></tr>{trows}</table>
<ul>
<li><b>rep_mixing</b>: within each genotype, how much rep1 and rep2 cells neighbour each other (k=30),
relative to their share; 1 = fully mixed, near 0 = replicates separate.</li>
<li><b>genotype_purity</b>: fraction of a cell's neighbours with the same genotype. Correction should not
drop this much, since ctr vs KO is the biology.</li>
<li>Per-sample correction treats each sample as a batch. Because each sample is a single genotype, it can
also pull ctr and KO together; compare genotype_purity against "no integration".</li>
<li>Normalization is the pipeline's log path, recomputed on each day's cells alone (HVGs and PCA per day).</li>
</ul>
{sections}
</main></body></html>"""
    (RES / "per_day_integration.html").write_text(html)
    print("wrote", RES, flush=True)


if __name__ == "__main__":
    main()
