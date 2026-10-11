#!/usr/bin/env python
"""
Merge RNA objects across samples for cross-sample integration -- matches
the lab's own notebook's merge step (3.single.timepoint.joint.analysis.
ipynb), generalized from 2 replicates to N samples: concatenates each
sample's step-3.1 final (SoupX-corrected, doublet-removed, RNA/ATAC
barcode-intersected) raw counts as one object, attaching diff_batch/
stage/tech_batch parsed from each sample's (possibly renamed, see
--sample-map-csv) sample_id via split("_") -- same convention as
5.1_coverage_table.py.

Deliberately uses 3.1's rna_final.h5ad (not 1.3's doublets_removed.h5ad
directly) so any per-sample manual barcode exclusion already applied at
merge time (3.1's --exclude-barcodes-dir, e.g. this project's
HZ_03_HZ021_D22_rep1 GEM-barcoding-failure correction) carries through
here for free, rather than needing to be re-applied.

Input counts are RAW (SoupX-corrected, not normalized) -- integration
needs to normalize ONCE on the combined object (next step,
5.3_normalize_merged.py), not merge already-independently-normalized
per-sample data, since independently-fit per-sample HVGs/SCT models
don't share a common feature basis to integrate on.

--sample-map-csv is optional: a CSV with old_name,new_name columns (this
project's rename from a real subfolder name to a study-specific
{diff_batch}_{stage}_{tech_batch} alias) -- when given, old_name is used
to find each sample's h5ad on disk, new_name is what gets parsed into
diff_batch/stage/tech_batch and used as the sample_id in the merged
object's obs. Without it, sample_id is used directly for both.

Usage:
  python 5.2_merge_rna.py --sample-ids-file samples.txt \
      --merged-dir /path/to/sc_pipeline_out/merged/h5ad \
      --outdir /path/to/out \
      [--sample-map-csv rename_manifest.csv] [--min-cells 20]
"""
import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def parse_sample_id(sample_id: str):
    parts = sample_id.split("_")
    if len(parts) != 3:
        raise ValueError(
            f"Expected sample_id '{sample_id}' to split into exactly 3 "
            f"underscore-separated parts (diff_batch, stage, tech_batch) -- got {parts}"
        )
    return parts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-ids-file", required=True, help="One sample_id per line (the {diff_batch}_{stage}_{tech_batch} name)")
    ap.add_argument("--merged-dir", required=True, help="3.1_merge_modalities.py's output h5ad dir")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sample-map-csv", default=None, help="Optional old_name,new_name CSV -- see module docstring")
    ap.add_argument("--min-cells", type=int, default=20, help="Gene filter after merge, matches the notebook")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sample_ids = [l.strip() for l in Path(args.sample_ids_file).read_text().splitlines() if l.strip()]

    old_for = {s: s for s in sample_ids}
    if args.sample_map_csv:
        rename = pd.read_csv(args.sample_map_csv)
        by_new = dict(zip(rename["new_name"], rename["old_name"]))
        missing = [s for s in sample_ids if s not in by_new]
        if missing:
            raise ValueError(f"--sample-map-csv given but missing new_name entries for: {missing}")
        old_for = {s: by_new[s] for s in sample_ids}

    merged_dir = Path(args.merged_dir)
    adatas = []
    for sid in sample_ids:
        old_name = old_for[sid]
        p = merged_dir / f"{old_name}.rna_final.h5ad"
        if not p.exists():
            raise FileNotFoundError(f"{p} not found for sample_id '{sid}' (old_name '{old_name}')")
        a = sc.read_h5ad(p)
        diff_batch, stage, tech_batch = parse_sample_id(sid)
        a.obs["sample_id"] = sid
        a.obs["diff_batch"] = diff_batch
        a.obs["stage"] = stage
        a.obs["tech_batch"] = tech_batch
        a.obs_names = [f"{sid}#{bc}" for bc in a.obs_names]
        print(f"[{sid}] (from {old_name}): {a.n_obs} cells, {a.n_vars} genes")
        adatas.append(a)

    merged = ad.concat(adatas, join="outer", merge="unique")
    merged.obs_names_make_unique()
    print(f"\nMerged: {merged.n_obs} cells x {merged.n_vars} genes (union of gene sets across samples)")

    sc.pp.filter_genes(merged, min_cells=args.min_cells)
    print(f"After filter_genes(min_cells={args.min_cells}): {merged.n_vars} genes")

    out_path = outdir / "merged_rna_raw.h5ad"
    merged.write_h5ad(out_path)

    summary = merged.obs.groupby(["diff_batch", "stage", "tech_batch"]).size().reset_index(name="n_cells")
    summary.to_csv(outdir / "merge_cell_counts.csv", index=False)
    print(f"\nWritten: {out_path}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
