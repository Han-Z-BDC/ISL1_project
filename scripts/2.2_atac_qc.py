#!/usr/bin/env python
"""
ATAC QC: SnapATAC2 feature space/general QC + a 3-detector doublet union
(SnapATAC2's own scrublet + AMULET + scDblFinder).

Run 2.1_atac_amulet.py (via scripts/submit_2.1_atac_amulet.sh) for every sample
BEFORE this -- AMULET is slow and runs independently of everything here,
so this script reads its cached {sample_id}.amulet_result.csv rather than
calling AMULET itself. See 2.1_atac_amulet.py for why.

Design (compared the lab's own notebook, IGVF-UCSD/CellCommander, and
sc-best-practices.org, confirmed with the user at each fork -- see
configs/atac_qc.yaml for the full reasoning):
  - Feature space: SnapATAC2's genome-wide tile matrix built directly from
    fragments.tsv.gz -- NOT CellRanger's own called peaks. The one point
    here where the chosen option matches only 1 of 3 references
    (IGVF-UCSD/CellCommander), a deliberate choice.
  - General QC (min_tsse, fragment count bounds): SnapATAC2's own
    algorithm and thresholds, following the IGVF-UCSD/CellCommander recipe
    -- these don't transfer from the muon-based metrics the other two
    references use, which are a different computation on a different
    feature space.
  - Doublets: union of THREE detectors -- SnapATAC2's own built-in
    scrublet (ATAC-native), AMULET (fragment-coverage-based, matches all
    three references), and scDblFinder with aggregateFeatures=TRUE
    (sc-best-practices' explicit fix for applying an RNA-simulation-based
    method to ATAC's sparsity). scDblFinder runs with clusters=FALSE, NOT
    sc-best-practices' own tutorial default of TRUE -- their own text says
    cluster-based scoring isn't recommended for continuous-trajectory data
    like this lab's differentiation time-course.
  - Generic scanpy Scrublet (used in RNA QC and the lab's own ATAC
    notebook) is NOT used here -- SnapATAC2's built-in scrublet is a
    different, ATAC-adapted implementation.

One sample per invocation -- called by 2.2_atac_qc.sb as a Slurm array
task, same pattern as the RNA steps.

Usage:
  python 2.2_atac_qc.py --sample-id D22_rep1 \
      --frag-path /path/to/D22_rep1/outs/atac_fragments.tsv.gz \
      --blacklist-bed /path/to/blacklist.bed \
      --config configs/atac_qc.yaml \
      --outdir /path/to/sc_pipeline_out/atac_qc
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import snapatac2 as snap
import anndata2ri
from rpy2 import robjects as ro
from rpy2.robjects.packages import importr

anndata2ri.activate()
anndata2ri.scipy2ri.activate()


def load_amulet_result(amulet_csv: Path, sample_id: str) -> pd.DataFrame:
    """Loads 2.1_atac_amulet.py's cached output. Returns a DataFrame
    indexed by RAW barcode (no sample prefix)."""
    if not amulet_csv.exists():
        raise FileNotFoundError(
            f"AMULET result not found at {amulet_csv} -- run "
            f"scripts/submit_2.1_atac_amulet.sh for sample '{sample_id}' first "
            f"(2.1_atac_amulet.py). AMULET is run separately because it's "
            f"slow and independent of QC thresholds -- see that script's "
            f"docstring."
        )
    return pd.read_csv(amulet_csv, index_col="barcode")


def run_scdblfinder_atac(count_mat, barcodes, cfg: dict) -> pd.DataFrame:
    """Runs scDblFinder on the tile/selected-feature matrix with
    aggregateFeatures=TRUE (sc-best-practices' recommended ATAC
    adaptation). Returns a DataFrame indexed by barcode."""
    importr("scDblFinder")
    SingleCellExperiment = importr("SingleCellExperiment")
    ro.r(f"set.seed({cfg['random_state']})")

    data_mat = count_mat.T.astype("float32")
    ro.globalenv["data_mat"] = data_mat
    clusters_r = "TRUE" if cfg["scdblfinder"]["clusters"] else "FALSE"
    ro.r(f"""
        sce <- SingleCellExperiment(list(counts=data_mat))
        sce <- scDblFinder(
            sce,
            clusters={clusters_r},
            aggregateFeatures={"TRUE" if cfg["scdblfinder"]["aggregate_features"] else "FALSE"},
            nfeatures={cfg["scdblfinder"]["n_features"]},
            processing="{cfg["scdblfinder"]["processing"]}"
        )
    """)
    col_data = anndata2ri.rpy2py(ro.r("as.data.frame(colData(sce))"))
    col_data.index = barcodes
    return pd.DataFrame(
        {
            "scdblfinder_doublet_score": col_data["scDblFinder.score"].values,
            "scdblfinder_predicted_doublet": (col_data["scDblFinder.class"].values == "doublet"),
        },
        index=barcodes,
    )


def plot_threshold_diagnostics(adata, cfg: dict, sample_id: str, out_png: Path):
    """The notebook-style threshold review: histograms of fragment count
    and TSSe, plus their joint scatter, with the CURRENT config's
    threshold lines drawn on top -- computed on the full pre-filter
    population (every imported barcode), so it's possible to actually see
    whether these lines separate real cells from background on THIS
    sample, not just trust the imported defaults. Meant to be looked at
    and iterated on, same as the lab's own notebook's threshold plots."""
    n_frag = adata.obs["n_fragment"].to_numpy()
    tsse = adata.obs["tsse"].to_numpy()
    log_n_frag = np.log10(n_frag + 1)

    min_frag, max_frag = cfg["min_num_fragments"], cfg["max_num_fragments"]
    min_tsse, max_tsse = cfg["min_tsse"], cfg.get("max_tsse")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    axes[0].hist(log_n_frag, bins=100, color="#999999")
    axes[0].axvline(np.log10(min_frag + 1), color="red", linestyle="--")
    if max_frag is not None:
        axes[0].axvline(np.log10(max_frag + 1), color="red", linestyle="--")
    axes[0].set_xlabel("log10(fragments + 1)")
    axes[0].set_ylabel("n barcodes")
    axes[0].set_title("Fragment count per barcode")

    tsse_plot_max = np.quantile(tsse, 0.995) if len(tsse) else 20
    axes[1].hist(tsse[tsse < tsse_plot_max], bins=100, color="#999999")
    axes[1].axvline(min_tsse, color="red", linestyle="--")
    if max_tsse is not None:
        axes[1].axvline(max_tsse, color="red", linestyle="--")
    axes[1].set_xlabel("TSS enrichment (SnapATAC2)")
    axes[1].set_title(f"TSSe, up to 99.5th pct ({tsse_plot_max:.1f})")

    keep_for_scatter = tsse < tsse_plot_max
    axes[2].scatter(log_n_frag[keep_for_scatter], tsse[keep_for_scatter], s=2, alpha=0.3, color="black")
    axes[2].axvline(np.log10(min_frag + 1), color="red", linestyle="--")
    if max_frag is not None:
        axes[2].axvline(np.log10(max_frag + 1), color="red", linestyle="--")
    axes[2].axhline(min_tsse, color="red", linestyle="--")
    if max_tsse is not None:
        axes[2].axhline(max_tsse, color="red", linestyle="--")
    axes[2].set_xlabel("log10(fragments + 1)")
    axes[2].set_ylabel("TSS enrichment")
    axes[2].set_title("Joint distribution -- current cutoffs in red")

    fig.suptitle(
        f"{sample_id}: review these cutoffs before trusting them "
        f"(min_num_fragments={min_frag}, max_num_fragments={max_frag}, "
        f"min_tsse={min_tsse}, max_tsse={max_tsse})"
    )
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_doublet_breakdown(obs: pd.DataFrame, sample_id: str, out_png: Path):
    fig, ax = plt.subplots(figsize=(5, 4))
    cols = ["snapatac2_predicted_doublet", "amulet_predicted_doublet", "scdblfinder_predicted_doublet"]
    counts = obs[cols].sum()
    counts.index = ["SnapATAC2\nscrublet", "AMULET", "scDblFinder"]
    ax.bar(counts.index, counts.values, color=["#4299e1", "#ed8936", "#805ad5"])
    ax.set_ylabel("n cells flagged")
    ax.set_title(f"{sample_id}: doublet calls by tool\n(union = {int(obs['doublet_filter_union'].sum())} cells removed)")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--frag-path", required=True, help="Path to atac_fragments.tsv.gz")
    ap.add_argument("--blacklist-bed", required=True)
    ap.add_argument("--config", required=True, help="Path to configs/atac_qc.yaml")
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--amulet-result", default=None,
        help="Path to {sample_id}.amulet_result.csv from 2.1_atac_amulet.py. "
             "Not needed with --qc-only (stops before doublet detection); "
             "required otherwise -- defaults to <outdir>/amulet/{sample_id}.amulet_result.csv.",
    )
    ap.add_argument(
        "--qc-only", action="store_true",
        help="Stop after writing threshold_review.png -- skips filter_cells, the tile "
             "matrix, feature selection, and all three doublet detectors (the slow "
             "part, especially AMULET). For fast per-sample threshold iteration: run "
             "this, look at the plot, edit configs/atac_qc.yaml, rerun -- only drop "
             "--qc-only once the cutoffs look right. On this lab's real fragment "
             "files, import_fragments+tsse is still real compute (confirmed via "
             "2.1/2.2's own measured runtimes) -- run via sbatch/acompile, NOT the "
             "bare login node.",
    )
    ap.add_argument(
        "--sample-cutoffs-csv", default=None,
        help="Optional CSV (sample_id,min_tsse,min_num_fragments,max_num_fragments,max_tsse) "
             "with per-sample general-QC threshold decisions, made after reviewing each "
             "sample's threshold_review.png (from a --qc-only run -- see "
             "scripts/submit_2.2_atac_qc_only.sh). Any non-empty cell overrides "
             "configs/atac_qc.yaml's global default for that sample only; a sample not "
             "in the CSV (or the CSV itself omitted) falls back to the global config "
             "untouched.",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "summary").mkdir(parents=True, exist_ok=True)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    if args.sample_cutoffs_csv and Path(args.sample_cutoffs_csv).exists():
        cutoffs_df = pd.read_csv(args.sample_cutoffs_csv, dtype=str).set_index("sample_id")
        if args.sample_id in cutoffs_df.index:
            row = cutoffs_df.loc[args.sample_id]
            for key, caster in [
                ("min_tsse", float), ("min_num_fragments", int),
                ("max_num_fragments", int), ("max_tsse", float),
            ]:
                val = row.get(key)
                if pd.notna(val) and str(val).strip() != "":
                    cfg[key] = caster(val)
            print(f"[{args.sample_id}] applied per-sample cutoff overrides from {args.sample_cutoffs_csv}: "
                  f"min_tsse={cfg['min_tsse']}, min_num_fragments={cfg['min_num_fragments']}, "
                  f"max_num_fragments={cfg['max_num_fragments']}, max_tsse={cfg.get('max_tsse')}")

    genome = getattr(snap.genome, cfg["chrom_sizes"])

    # --- SnapATAC2 backbone: import, TSSe ---
    # This is the expensive part (parses the whole fragment file) -- cached
    # to a checkpoint h5ad so a rerun (--qc-only while tuning thresholds,
    # or the eventual full run once they're set) never repeats it. Delete
    # the checkpoint file to force a clean re-import (e.g. if the fragment
    # file itself changed).
    checkpoint_path = outdir / "h5ad" / f"{args.sample_id}.import_tsse_checkpoint.h5ad"
    if checkpoint_path.exists():
        print(f"[{args.sample_id}] Reusing existing checkpoint: {checkpoint_path} "
              f"(delete it to force a clean re-import from fragments)")
        adata = snap.read(str(checkpoint_path), backed=None)
    else:
        adata = snap.pp.import_fragments(
            args.frag_path,
            chrom_sizes=genome,
            min_num_fragments=cfg["min_load_num_fragments"],
            sorted_by_barcode=cfg["sorted_by_barcode"],
            chunk_size=cfg["chunk_size"],
        )
        adata.obs_names = [f"{args.sample_id}#{bc}" for bc in adata.obs_names]
        snap.metrics.tsse(adata, genome)
        adata.write(str(checkpoint_path))

    # raw (unprefixed) barcodes, for AMULET's barcode matching later --
    # re-derived rather than kept from import so it also works when adata
    # came from the checkpoint.
    prefix = f"{args.sample_id}#"
    raw_barcodes = [bc[len(prefix):] for bc in adata.obs_names]

    snap.pl.frag_size_distr(adata, interactive=False, out_file=str(outdir / "plots" / f"{args.sample_id}.frag_size_distr.png"))

    # Threshold review, notebook-style: the actual metric spread with the
    # CURRENT config's cutoff lines drawn on top. Look at this before
    # trusting the filter on a new dataset -- these thresholds are fixed
    # values (IGVF-UCSD's own, not yet validated on this lab's data), not
    # something derived per-sample, so there is no substitute for actually
    # looking at where the lines fall.
    plot_threshold_diagnostics(adata, cfg, args.sample_id, outdir / "plots" / f"{args.sample_id}.threshold_review.png")

    n_before_qc = adata.n_obs
    qc_mask = snap.pp.filter_cells(
        adata,
        min_counts=cfg["min_num_fragments"],
        min_tsse=cfg["min_tsse"],
        max_counts=cfg["max_num_fragments"],
        max_tsse=cfg.get("max_tsse"),
        inplace=not args.qc_only,
    )

    if args.qc_only:
        n_after_qc = int(qc_mask.sum())
        print(f"[{args.sample_id}] --qc-only: {n_before_qc} -> {n_after_qc} cells would "
              f"pass general QC at the current thresholds in {args.config}. "
              f"Review {args.sample_id}.threshold_review.png; if the lines look wrong, "
              f"edit the config and rerun with --qc-only again (reuses the checkpoint, "
              f"so this stays fast). Drop --qc-only once you're satisfied to run the "
              f"full pipeline (tile matrix, feature selection, doublet detection).")
        return

    n_after_qc = adata.n_obs

    snap.pp.add_tile_matrix(adata, bin_size=cfg["bin_size"])
    snap.pp.select_features(adata, n_features=cfg["num_features"], blacklist=args.blacklist_bed)

    # --- doublet detector 1: SnapATAC2's own built-in scrublet ---
    # features= passed explicitly as a numpy array, NOT left at the
    # default (the string "selected", which snap.pp.scrublet's own
    # internal code resolves via `adata.var["selected"]` -- a pandas
    # Series). Real crash on this lab's data: snapatac2 2.8.0's internal
    # `adata.X[:, features]` needs a numpy array for scipy 1.15's sparse
    # indexing; a pandas Series raises "'Series' object has no attribute
    # 'nonzero'" (scipy removed that path; snapatac2 hasn't caught up --
    # confirmed as an open, unresolved upstream bug, kaizhang/SnapATAC2#386,
    # workaround there is downgrading scipy, which isn't worth doing to the
    # shared env: every other completed step here already validated
    # against scipy 1.15). features= itself is typed 'str | np.ndarray |
    # None' in snap.pp.scrublet's own signature, so passing a numpy array
    # directly is a documented, intended path -- not a hack -- it just
    # skips the buggy string-to-Series lookup entirely. No notebook
    # precedent for this call at all (the lab's own notebook never uses
    # SnapATAC2's built-in scrublet, only generic scanpy's sce.pp.scrublet
    # -- see 1.3_detect_doublets.py / RNA QC), so this is a from-scratch
    # library-version fix, not a methods disagreement to reconcile.
    snap.pp.scrublet(adata, features=adata.var["selected"].to_numpy(), random_state=cfg["random_state"])
    adata.obs["snapatac2_predicted_doublet"] = adata.obs["doublet_probability"] >= 0.5

    # --- doublet detector 2: AMULET -- cached result from 2.1_atac_amulet.py, not run here ---
    amulet_csv = Path(args.amulet_result) if args.amulet_result else outdir / "amulet" / f"{args.sample_id}.amulet_result.csv"
    amulet_result = load_amulet_result(amulet_csv, args.sample_id)
    amulet_result.index = [f"{args.sample_id}#{bc}" for bc in amulet_result.index]
    amulet_result = amulet_result.reindex(adata.obs_names)
    adata.obs["amulet_qVal"] = amulet_result["q.value"].to_numpy()
    adata.obs["amulet_predicted_doublet"] = (
        adata.obs["amulet_qVal"] < cfg["amulet"]["q_value_threshold"]
    ).fillna(False)

    # --- doublet detector 3: scDblFinder, aggregateFeatures=TRUE, on the selected-feature tile matrix ---
    selected = adata.var["selected"].to_numpy() if "selected" in adata.var.columns else np.ones(adata.n_vars, dtype=bool)
    count_mat = adata.X[:, selected]
    scdbl = run_scdblfinder_atac(count_mat, adata.obs_names, cfg)
    adata.obs["scdblfinder_doublet_score"] = scdbl["scdblfinder_doublet_score"].to_numpy()
    adata.obs["scdblfinder_predicted_doublet"] = scdbl["scdblfinder_predicted_doublet"].to_numpy()

    # --- union ---
    doublet_cols = ["snapatac2_predicted_doublet", "amulet_predicted_doublet", "scdblfinder_predicted_doublet"]
    adata.obs["doublet_filter_union"] = adata.obs[doublet_cols].any(axis=1)
    adata.obs["doublet_filter_all_agree"] = adata.obs[doublet_cols].all(axis=1)

    if cfg["save_diagnostic_plot"]:
        plot_doublet_breakdown(adata.obs, args.sample_id, outdir / "plots" / f"{args.sample_id}.doublet_breakdown.png")

    # --- diagnostic clustering/embedding (not consumed by any later step) ---
    snap.tl.spectral(adata)
    n_dims = adata.obsm["X_spectral"].shape[1]
    snap.tl.umap(adata, use_dims=list(range(1, n_dims)))
    snap.pp.knn(adata, use_rep="X_spectral", use_dims=list(range(1, n_dims)))
    leiden_key = f"leiden_{cfg['clustering_resolution']}"
    snap.tl.leiden(adata, resolution=cfg["clustering_resolution"], key_added=leiden_key)

    with plt.rc_context({"figure.figsize": (5, 5)}):
        import scanpy as sc
        sc.pl.embedding(adata.to_memory() if hasattr(adata, "to_memory") else adata,
                         basis="X_umap", color=["doublet_filter_union", leiden_key], show=False)
        plt.savefig(outdir / "plots" / f"{args.sample_id}.umap_doublets_clusters.png", dpi=150)
        plt.close()

    if cfg["build_gene_activity_matrix"]:
        gene_matrix = snap.pp.make_gene_matrix(adata, gene_anno=genome)
        gene_matrix.write(str(outdir / "h5ad" / f"{args.sample_id}.gene_activity.h5ad"))

    adata.write(str(outdir / "h5ad" / f"{args.sample_id}.atac_qc_annotated.h5ad"))

    keep = ~adata.obs["doublet_filter_union"].to_numpy()
    adata_final = adata[keep, :]
    adata_final.write(str(outdir / "h5ad" / f"{args.sample_id}.atac_qc_final.h5ad"))
    n_final = int(keep.sum())

    summary = {
        "sample_id": args.sample_id,
        "n_cells_after_import": len(raw_barcodes),
        "n_cells_after_general_qc": n_after_qc,
        "n_cells_before_general_qc": n_before_qc,
        "n_snapatac2_scrublet_doublets": int(adata.obs["snapatac2_predicted_doublet"].sum()),
        "n_amulet_doublets": int(adata.obs["amulet_predicted_doublet"].sum()),
        "n_scdblfinder_doublets": int(adata.obs["scdblfinder_predicted_doublet"].sum()),
        "n_all_3_agree": int(adata.obs["doublet_filter_all_agree"].sum()),
        "n_removed_union": int((~keep).sum()),
        "n_cells_final": n_final,
    }
    pd.DataFrame([summary]).to_csv(outdir / "summary" / f"{args.sample_id}.atac_qc_summary.csv", index=False)

    print(f"[{args.sample_id}] {n_before_qc} -> {n_after_qc} cells after general QC -> "
          f"{n_final} cells after doublet union removal.")


if __name__ == "__main__":
    main()
