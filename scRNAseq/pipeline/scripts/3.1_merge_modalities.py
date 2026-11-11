#!/usr/bin/env python
"""
Merge RNA and ATAC QC: intersect each sample's QC-passing barcodes across
both modalities, and subset both objects down to the shared cells --
treated as just another filter step, same as every step before it: reads
two h5ad in, writes two h5ad out. No MuData bundling -- if a joint
container turns out to be useful later, build it then (one line:
mu.MuData({"rna": rna, "atac": atac})); not needed now.

Matches the lab's own notebook exactly
(3.single.timepoint.joint.analysis.ipynb, "Merging RNA and ATAC to keep
good quality cells from both modalities"): plain barcode-set intersection,
no translation table -- 10x Multiome shares one barcode whitelist between
the GEX and ATAC libraries of the same GEM well, so a real cell has the
identical barcode string in both modalities.

This pipeline's ATAC QC (2.2_atac_qc.py) prefixes barcodes with
"{sample_id}#" (a SnapATAC2/CellCommander convention, for safe
multi-sample concatenation later) -- that prefix is stripped before
intersecting, since the lab's own RNA QC output never uses one.

Run this AFTER both RNA QC (1.3_detect_doublets.py's
{sample}.doublets_removed.h5ad) and ATAC QC (2.2_atac_qc.py's
{sample}.atac_qc_final.h5ad) are done for a sample.

Optional --exclude-barcodes-dir: a directory of per-sample manual barcode
exclusion lists ({sample_id}.excluded_barcodes.txt, one barcode per line,
no header), applied AFTER the RNA/ATAC intersection so both modalities'
final output stay barcode-consistent. For a documented, one-off technical
problem identified downstream (e.g. an isolated, non-specifically-
co-expressing cluster on the normalize-step UMAP, traced to a known
wet-lab GEM barcoding failure) -- NOT a general QC step; that's what
1.1-2.2 already do. A sample with no file here is unaffected. Deliberately
applied HERE, not by hand-editing the h5ad in place, so this correction
survives a future rerun of this step instead of silently reverting.

Scope note: the lab's own notebook goes further from here -- merging
replicates, normalizing (log1p AND SCTransform), regressing cell-cycle +
MT%/Ribo%, and comparing four integration methods (CCA/RPCA/Harmony/
FastMNN). That's a future normalization + integration step's job, not
built yet -- this script only does the per-sample barcode reconciliation,
matching what was actually asked for.

Usage:
  python 3.1_merge_modalities.py --sample-id D22_rep1 \
      --rna-h5ad /path/to/detect_doublets/h5ad/D22_rep1.doublets_removed.h5ad \
      --atac-h5ad /path/to/atac_qc/h5ad/D22_rep1.atac_qc_final.h5ad \
      --outdir /path/to/sc_pipeline_out/merged
"""
import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def strip_atac_prefix(adata, sample_id: str):
    """This pipeline's own ATAC QC always prefixes barcodes with
    "{sample_id}#" -- strip it back off before comparing to RNA barcodes,
    which never carry one. Fails loudly if the assumption doesn't hold,
    rather than silently producing an empty intersection."""
    prefix = f"{sample_id}#"
    if not all(bc.startswith(prefix) for bc in adata.obs_names):
        raise ValueError(
            f"Expected every ATAC barcode to start with '{prefix}' (this "
            f"pipeline's own 2.2_atac_qc.py convention) but some don't -- "
            f"was --atac-h5ad produced by a different tool/version? "
            f"Example barcodes: {list(adata.obs_names[:5])}"
        )
    adata.obs_names = [bc[len(prefix):] for bc in adata.obs_names]
    return adata


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--rna-h5ad", required=True, help="1.3_detect_doublets.py output: {sample}.doublets_removed.h5ad")
    ap.add_argument("--atac-h5ad", required=True, help="2.2_atac_qc.py output: {sample}.atac_qc_final.h5ad")
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--exclude-barcodes-dir", default=None,
        help="Optional dir of {sample_id}.excluded_barcodes.txt manual exclusion lists -- "
             "see module docstring. A sample with no matching file is unaffected.",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "h5ad").mkdir(parents=True, exist_ok=True)
    (outdir / "summary").mkdir(parents=True, exist_ok=True)

    rna = sc.read_h5ad(args.rna_h5ad)
    atac = sc.read_h5ad(args.atac_h5ad)
    atac = strip_atac_prefix(atac, args.sample_id)

    n_rna_before = rna.n_obs
    n_atac_before = atac.n_obs
    print(f"[{args.sample_id}] RNA before merge: {n_rna_before}")
    print(f"[{args.sample_id}] ATAC before merge: {n_atac_before}")

    rna_barcodes = set(rna.obs_names)
    atac_barcodes = set(atac.obs_names)
    common = sorted(rna_barcodes & atac_barcodes)

    pct_of_rna = 100 * len(common) / n_rna_before if n_rna_before else 0.0
    pct_of_atac = 100 * len(common) / n_atac_before if n_atac_before else 0.0

    if pct_of_rna < 50 or pct_of_atac < 50:
        print(f"[{args.sample_id}] WARNING: low barcode overlap "
              f"({pct_of_rna:.1f}% of RNA, {pct_of_atac:.1f}% of ATAC retained). "
              f"This usually signals a barcode format mismatch, not genuinely "
              f"poor multiome capture -- compare examples before trusting this:")
        print(f"  RNA examples:  {sorted(rna_barcodes)[:5]}")
        print(f"  ATAC examples: {sorted(atac_barcodes)[:5]}")

    rna = rna[common].copy()
    atac = atac[common].copy()

    n_manually_excluded = 0
    manual_excl_file = None
    if args.exclude_barcodes_dir:
        manual_excl_file = Path(args.exclude_barcodes_dir) / f"{args.sample_id}.excluded_barcodes.txt"
        if manual_excl_file.exists():
            excluded_bcs = set(manual_excl_file.read_text().split())
            not_found = excluded_bcs - set(rna.obs_names)
            if not_found:
                print(f"[{args.sample_id}] WARNING: {len(not_found)} of {len(excluded_bcs)} "
                      f"manual-exclusion barcodes in {manual_excl_file} aren't in the merged "
                      f"set (already dropped by QC/intersection, or a stale list) -- ignoring those.")
            keep_mask = ~rna.obs_names.isin(excluded_bcs)
            n_manually_excluded = int((~keep_mask).sum())
            rna = rna[keep_mask].copy()
            atac = atac[keep_mask].copy()
            print(f"[{args.sample_id}] Manually excluded {n_manually_excluded} cells "
                  f"per {manual_excl_file} (documented technical problem, not general QC).")

    # Same convention as every step before this one -- plain h5ad, one per
    # modality. What a future step (normalization, per-modality analysis)
    # should read.
    out_rna = outdir / "h5ad" / f"{args.sample_id}.rna_final.h5ad"
    out_atac = outdir / "h5ad" / f"{args.sample_id}.atac_final.h5ad"
    rna.write_h5ad(out_rna)
    atac.write_h5ad(out_atac)

    summary = {
        "sample_id": args.sample_id,
        "n_rna_before_merge": n_rna_before,
        "n_atac_before_merge": n_atac_before,
        "n_shared_cells": len(common),
        "n_rna_only_dropped": n_rna_before - len(common),
        "n_atac_only_dropped": n_atac_before - len(common),
        "pct_of_rna_retained": round(pct_of_rna, 2),
        "pct_of_atac_retained": round(pct_of_atac, 2),
        "n_manually_excluded": n_manually_excluded,
        "manual_exclusion_file": str(manual_excl_file) if n_manually_excluded else "",
        "n_final_cells": rna.n_obs,
    }
    pd.DataFrame([summary]).to_csv(outdir / "summary" / f"{args.sample_id}.merge_summary.csv", index=False)

    print(f"[{args.sample_id}] {rna.n_obs} cells in final output "
          f"({len(common)} shared by RNA+ATAC, {n_manually_excluded} further manually excluded). "
          f"Written to {out_rna} and {out_atac}")


if __name__ == "__main__":
    main()
