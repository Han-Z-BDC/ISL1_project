"""Step 7.2a: cells for the D15 + D22 endocrine integration (first of 7.2a-d).

Request (user, 2026-11-19): integrate the D15 endocrine cells (the 7.1d subset) with the D22 endocrine cells
(every 7.1h label except non-endocrine). Merge, SCT B with cell-cycle regression only, batch correction on stage
(D15 vs D22) with RPCA, Harmony and scVI; genotype is not corrected. No clustering; report with marker genes.

Cells:  D15  day_analysis/D15/subcluster_endocrine/cells.txt (7.1d; r0.2 clusters 0, 1, 5); draft labels from
             study_config/annotation/D15_endocrine_r0.6.tsv on subclusters.csv sub_r0p6
        D22  day_analysis/D22/labels/labels.csv (7.1h), every label except "non-endocrine"
Counts: layers/soupx_counts of each day's merged h5ad (7.1; QC, doublets and NF > 0.6 already removed).
Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/ (reruns overwrite):
  merged_counts.h5ad  X = SoupX counts (CSR), obs = stage, sample, genotype, label, S / G2M scores (from each
                      day's 7.1 object)
  cells.tsv           the same obs as a table (read by 7.2b)
  composition.csv     cells per stage x genotype x label
Subset (optional): with INT_SUBSET=<name> set, only the cells whose label is listed in
  study_config/annotation/integration_subset_<name>.txt are kept, and every 7.2a-d output goes to
  ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine_<name>/ instead (the full run is untouched).
Run: bash analysis/scripts/7.2_integrate_submit.sh [subset]   (submits 7.2a-d in order)
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

OUT_DIR, REPO_DIR = Path(os.environ["OUT_DIR"]), Path(os.environ["REPO_DIR"])
DA = OUT_DIR / "day_analysis"
SUBSET = os.environ.get("INT_SUBSET", "")
OUT = DA / ("integration_D15_D22_endocrine" + (f"_{SUBSET}" if SUBSET else ""))
OBS = ["sample", "genotype", "S_score", "G2M_score", "S_minus_G2M", "phase"]


def take(day, cells, labels):
    a = ad.read_h5ad(DA / day / f"{day}_merged.h5ad")
    pos = pd.Index(a.obs_names).get_indexer(cells)
    assert (pos >= 0).all(), f"{day}: cells missing from the merged object"
    x = a.layers["soupx_counts"][pos]
    x = sp.csr_matrix(x).astype(np.float32)
    obs = a.obs.iloc[pos][OBS].copy()
    obs["sample"] = day + "_" + obs["sample"].astype(str)
    obs["genotype"] = np.where(obs.genotype.astype(str).str.upper().str.startswith("KO"), "KO", "ctr")
    obs.insert(0, "stage", day)
    obs["label"] = labels
    return ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=a.var_names.copy()))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in OUT.glob("*"):
        if p.is_file():
            p.unlink()
    # D15: 7.1d subset with its draft r0.6 labels
    sub = DA / "D15" / "subcluster_endocrine"
    cells15 = [l.strip() for l in open(sub / "cells.txt") if l.strip()]
    sc15 = pd.read_csv(sub / "subclusters.csv", index_col=0, dtype=str).reindex(cells15)
    lab15 = pd.read_csv(REPO_DIR / "study_config" / "annotation" / "D15_endocrine_r0.6.tsv", sep="\t", dtype=str)
    m15 = dict(zip(lab15.subcluster, lab15.label.str.split(" \\(").str[0]))
    l15 = sc15.sub_r0p6.map(m15)
    assert l15.notna().all(), "D15 subcluster without a draft label"
    # D22: every label except non-endocrine
    lab22 = pd.read_csv(DA / "D22" / "labels" / "labels.csv", index_col=0)
    lab22 = lab22[lab22.label != "non-endocrine"]

    a15 = take("D15", cells15, ("D15 " + l15).to_numpy())
    a22 = take("D22", list(lab22.index), lab22.label.astype(str).to_numpy())
    assert list(a15.var_names) == list(a22.var_names), "the two days have different gene lists"
    a = ad.concat([a15, a22], join="inner", merge="same")
    if SUBSET:
        f = REPO_DIR / "study_config" / "annotation" / f"integration_subset_{SUBSET}.txt"
        keep = [l.strip() for l in open(f) if l.strip() and not l.startswith("#")]
        missing = set(keep) - set(a.obs.label)
        assert not missing, f"labels in {f.name} not found: {missing}"
        print(f"subset {SUBSET}: keeping {len(keep)} labels; dropped:",
              a.obs.label[~a.obs.label.isin(keep)].value_counts().to_dict())
        a = a[a.obs.label.isin(keep).to_numpy()].copy()
    for c in ["stage", "sample", "genotype", "label", "phase"]:
        a.obs[c] = a.obs[c].astype(str).astype("category")
    assert np.allclose(a.X.data, np.round(a.X.data))
    a.write_h5ad(OUT / "merged_counts.h5ad")
    a.obs.to_csv(OUT / "cells.tsv", sep="\t", index_label="barcode")
    comp = a.obs.groupby(["stage", "label", "genotype"], observed=True).size().unstack(fill_value=0)
    comp.to_csv(OUT / "composition.csv")
    print(a); print(comp.to_string()); print("wrote", OUT)


if __name__ == "__main__":
    main()
