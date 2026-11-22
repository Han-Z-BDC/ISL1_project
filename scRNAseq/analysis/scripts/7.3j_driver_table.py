"""Step 7.3j: one alpha / beta / EC driver table for D15, and the gene set for the subset tradeSeq (7.3k) (after 7.3d, 7.3h).

Decisions (user, 2026-11-21): alpha, beta and EC driver genes in one table, then which genes differ by genotype along
the alpha, beta and EC trajectories (7.3k). No CellRank rerun: the drivers come from the control-only runs of 7.3d,
which already computed every fixed terminal state (the all-cell runs only did alpha and beta).
  driver      r >= 0.2 and FDR < 0.05 toward the fate in the control early cells (ENP1, ENP2, ENP_FEV_high) of the
              control-only CellRank runs (7.3d cellrank_{kernel}_ctr_drivers.csv), as 7.3f
              beta: fate pre_beta2; alpha: fate pre_alpha; EC: fate pre_EC
              program = both (velocity and pseudotime kernels) / one (one kernel) / -
  gene set    any driver (both or one kernel, any of the three fates)
              + the 2,656 KO vs control trajectory genes classified in 7.3h (gene_patterns_k10.csv)
              + the endocrine markers of 7.3e
Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
  drivers_alpha_beta_EC.csv      gene x {alpha, beta, EC}: r and FDR per kernel, program; n_fates (fates with a
                                 driver call), fates (which)
  tradeseq_subset_genes.csv      gene, driver_alpha / driver_beta / driver_EC (program), ko_vs_ctr (7.3h cluster or
                                 -), marker (True / False)
Run: python 7.3j_driver_table.py   (TOOLS python; small tables, seconds)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

TR = Path(os.environ["OUT_DIR"]) / "day_analysis" / "D15" / "trajectory"
R_MIN, FDR_MAX = 0.2, 0.05
FATES = {"alpha": "pre_alpha", "beta": "pre_beta2", "EC": "pre_EC"}
KERNELS = ["velocity", "pseudotime"]
MARKERS = ["CHGA", "NEUROG3", "NKX2-2", "NEUROD1", "INSM1", "PAX4", "ARX", "IRX2", "PAX6", "ISL1", "MAFB", "MAFA",
           "INS", "GCG", "SST", "GHRL", "PDX1", "NKX6-1", "FEV", "LMX1A", "TPH1", "PHOX2A", "HHEX", "MNX1", "TOX"]

cols = {}
for k in KERNELS:
    d = pd.read_csv(TR / f"cellrank_{k}_ctr_drivers.csv")
    assert d.cells.nunique() == 1, d.cells.unique()
    for name, fate in FATES.items():
        s = d[d.lineage == fate].set_index("gene")
        assert len(s), (k, fate)
        cols[(name, f"r_{k}")] = s.r
        cols[(name, f"FDR_{k}")] = s.FDR
t = pd.DataFrame(cols)
out = pd.DataFrame(index=t.index)
for name in FATES:
    n = sum(((t[(name, f"r_{k}")] >= R_MIN) & (t[(name, f"FDR_{k}")] < FDR_MAX)).astype(int) for k in KERNELS)
    for k in KERNELS:
        out[f"{name}_r_{k}"] = t[(name, f"r_{k}")].round(4)
        out[f"{name}_FDR_{k}"] = t[(name, f"FDR_{k}")]
    out[f"{name}_program"] = np.select([n == 2, n == 1], ["both", "one"], "-")
prog = out[[f"{n}_program" for n in FATES]]
out["n_fates"] = (prog != "-").sum(1)
out["fates"] = prog.apply(lambda r: ",".join(n for n in FATES if r[f"{n}_program"] != "-"), axis=1)
out.index.name = "gene"
drv = out[out.n_fates > 0].copy()
drv["best_r"] = drv[[c for c in drv if "_r_" in c]].max(1)
drv = drv.sort_values(["n_fates", "best_r"], ascending=False).drop(columns="best_r")
drv.to_csv(TR / "drivers_alpha_beta_EC.csv")

print("drivers (r >= 0.2, FDR < 0.05, control-only runs):")
for n in FATES:
    print(f"  {n:5s} both kernels {(drv[f'{n}_program'] == 'both').sum():4d}   one kernel "
          f"{(drv[f'{n}_program'] == 'one').sum():4d}")
print("  genes driving 1 / 2 / 3 fates:", drv.n_fates.value_counts().sort_index().to_dict())
print("  fate combinations:", drv.fates.value_counts().to_dict())

pat = pd.read_csv(TR / "gene_patterns_k10.csv")
assert {"gene", "cluster"} <= set(pat.columns), pat.columns
pat = pat.drop_duplicates("gene").set_index("gene")
genes = sorted(set(drv.index) | set(pat.index) | (set(MARKERS) & set(out.index)))
gs = pd.DataFrame(index=pd.Index(genes, name="gene"))
for n in FATES:
    gs[f"driver_{n}"] = out[f"{n}_program"].reindex(gs.index).fillna("-")
gs["ko_vs_ctr"] = pat.cluster.reindex(gs.index).map(lambda c: "-" if pd.isna(c) else f"C{int(c)}")
gs["marker"] = gs.index.isin(MARKERS)
gs.to_csv(TR / "tradeseq_subset_genes.csv")
d_any = (gs[[f"driver_{n}" for n in FATES]] != "-").any(axis=1)
print(f"tradeSeq gene set: {len(gs)} genes = drivers {d_any.sum()} | KO vs control {(gs.ko_vs_ctr != '-').sum()} "
      f"| overlap {(d_any & (gs.ko_vs_ctr != '-')).sum()} | markers only {(gs.marker & ~d_any & (gs.ko_vs_ctr == '-')).sum()}")
print("wrote", TR / "drivers_alpha_beta_EC.csv", "and", TR / "tradeseq_subset_genes.csv")
