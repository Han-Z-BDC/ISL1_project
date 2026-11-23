"""Step 7.1o (cells, layer 2): per-cell beta / EC identity scores of the D22 EC cells and control beta cells (after
7.1o_ec_beta_model.R model).

Question (user, 2026-11-23): do some KO EC keep beta features? A pseudobulk model only compares averages, so a minority
of beta-like KO EC would be diluted; this layer looks at single cells.
  cells     D22 SC_EC1 / SC_EC2 of both genotypes and control SC_beta (7.1h labels)
  scores    scanpy score_genes on the unregressed log1p X (seed 1234, ctrl_size 50): beta score from the 7.1o beta-high
            genes, EC score from the EC-high genes (score_genes.csv; INS left out); beta index = beta - EC score
  beta-leaning KO EC   KO cells of an EC type with a beta index above the 99th percentile (main) / 95th percentile
            (check) of the control cells of the same EC type, i.e. beyond the control EC range
  neighbours  share of each cell's neighbours (SCT B kNN graph of D22, 7.1b) that are control SC_beta
  INS       INS per 10,000 SoupX counts (as fig_ins_distribution.R), Spearman with the beta index within KO EC
  pseudobulk  SoupX counts summed per KO sample x EC type x (beta-leaning 99th / other) and per control sample for
            SC_beta, input of the R subset mode
  all cells  the same beta / EC scores computed on every D22 cell (user, 2026-11-23), for the whole-day UMAPs
Writes ${OUT_DIR}/day_analysis/ec_beta_D22/: cells.csv.gz, cells_all.csv.gz, leaning_counts.csv, pb_subset_counts.csv.gz,
  pb_subset_samples.csv, cells_summary.csv
Run: sbatch analysis/scripts/7.1o_ec_beta.sb
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr

OUT_DIR = Path(os.environ["OUT_DIR"])
DD = OUT_DIR / "day_analysis" / "D22"
OD = OUT_DIR / "day_analysis" / "ec_beta_D22"
SEED, Q_MAIN, Q_CHECK = 1234, 0.99, 0.95

lab = pd.read_csv(DD / "labels" / "labels.csv", index_col=0, keep_default_na=False)
a = ad.read_h5ad(DD / "D22_merged.h5ad")
a = a[lab.index].copy()
a.obs["label"] = lab.label.to_numpy()
keep = a.obs.label.isin(["SC_EC1", "SC_EC2"]) | ((a.obs.label == "SC_beta") & (a.obs.genotype == "ctr"))
# neighbours on the full graph before subsetting
G = a.obsp["SCT_B_connectivities"] if "SCT_B_connectivities" in a.obsp else a.obsp["connectivities"]
is_beta = ((a.obs.label == "SC_beta") & (a.obs.genotype == "ctr")).to_numpy().astype(float)
Gb = (G > 0).astype(float)
nb_beta = np.asarray(Gb @ is_beta).ravel() / np.maximum(np.asarray(Gb.sum(1)).ravel(), 1)
a.obs["beta_neighbours"] = nb_beta
sg = pd.read_csv(OD / "score_genes.csv")
bg, eg = [g for g in sg[sg.set == "beta"].gene if g in a.var_names], [g for g in sg[sg.set == "EC"].gene if g in a.var_names]
# all D22 cells scored too (user, 2026-11-23: scores on the whole-day UMAP); separate output, the subset scores below
# (and their thresholds) are computed on the subset as before
al = a.copy()
sc.tl.score_genes(al, bg, score_name="beta_score", ctrl_size=50, random_state=SEED)
sc.tl.score_genes(al, eg, score_name="EC_score", ctrl_size=50, random_state=SEED)
al.obs[["sample", "genotype", "label", "beta_score", "EC_score"]].to_csv(OD / "cells_all.csv.gz")
print("all D22 cells, median scores by label x genotype:\n",
      al.obs.groupby(["label", "genotype"], observed=True)[["beta_score", "EC_score"]].median().round(3).to_string())
del al
a = a[keep.to_numpy()].copy()
a.obs["group"] = a.obs.genotype.astype(str) + "_" + a.obs.label.str.replace("SC_", "", regex=False)
sc.tl.score_genes(a, bg, score_name="beta_score", ctrl_size=50, random_state=SEED)
sc.tl.score_genes(a, eg, score_name="EC_score", ctrl_size=50, random_state=SEED)
a.obs["beta_index"] = a.obs.beta_score - a.obs.EC_score
sx = a.layers["soupx_counts"]
tot = np.asarray(sx.sum(1)).ravel()
ins = np.asarray(sx[:, a.var_names.get_loc("INS")].todense()).ravel()
a.obs["INS_per10k"] = ins / np.maximum(tot, 1) * 1e4

o = a.obs[["sample", "genotype", "label", "group", "beta_score", "EC_score", "beta_index", "INS_per10k",
           "beta_neighbours", "scDblFinder_doublet_score", "n_genes_by_counts"]].copy()
counts, summ = [], []
for k in ["EC1", "EC2"]:
    c = o[o.group == f"ctr_{k}"].beta_index
    for q, tag in [(Q_MAIN, "lean99"), (Q_CHECK, "lean95")]:
        t = c.quantile(q)
        m = (o.group == f"KO_{k}") & (o.beta_index > t)
        o.loc[o.group == f"KO_{k}", tag] = m[o.group == f"KO_{k}"]
        for s_ in sorted(o[o.group == f"KO_{k}"]["sample"].unique()):
            n = ((o.group == f"KO_{k}") & (o["sample"] == s_)).sum()
            counts.append({"ec": k, "cut": tag, "threshold": round(t, 4), "sample": s_, "cells": int(n),
                           "beta_leaning": int((m & (o["sample"] == s_)).sum())})
    ko = o[o.group == f"KO_{k}"]
    for g_ in [f"ctr_{k}", f"KO_{k}", "ctr_beta"]:
        v = o[o.group == g_]
        summ.append({"group": g_, "cells": len(v), "beta_score_median": v.beta_score.median(),
                     "EC_score_median": v.EC_score.median(), "beta_index_median": v.beta_index.median(),
                     "beta_index_p95": v.beta_index.quantile(.95), "beta_neighbours_mean": v.beta_neighbours.mean(),
                     "INS_per10k_median": v.INS_per10k.median()})
    r = spearmanr(ko.beta_index, ko.INS_per10k).correlation
    summ.append({"group": f"KO_{k}: Spearman beta index vs INS", "cells": len(ko), "beta_index_median": r})
o.to_csv(OD / "cells.csv.gz")
lc = pd.DataFrame(counts); print(lc.to_string(index=False)); lc.to_csv(OD / "leaning_counts.csv", index=False)
sm = pd.DataFrame(summ).drop_duplicates("group"); print(sm.round(3).to_string(index=False))
sm.to_csv(OD / "cells_summary.csv", index=False)

# pseudobulks for the subset DE (main cut)
cols, meta = {}, []
for k in ["EC1", "EC2"]:
    for s_ in sorted(o[o.group == f"KO_{k}"]["sample"].unique()):
        for lean, nm in [(True, "beta_leaning"), (False, "other")]:
            m = ((o.group == f"KO_{k}") & (o["sample"] == s_) & (o.lean99.fillna(False) == lean)).to_numpy()
            pb = f"{s_}|KO_{k}_{nm}"
            cols[pb] = np.asarray(sx[m].sum(0)).ravel(); meta.append({"pb": pb, "sample": s_, "group": f"KO_{k}_{nm}", "cells": int(m.sum())})
for s_ in sorted(o[o.group == "ctr_beta"]["sample"].unique()):
    m = ((o.group == "ctr_beta") & (o["sample"] == s_)).to_numpy()
    pb = f"{s_}|ctr_beta"; cols[pb] = np.asarray(sx[m].sum(0)).ravel()
    meta.append({"pb": pb, "sample": s_, "group": "ctr_beta", "cells": int(m.sum())})
pd.DataFrame(cols, index=a.var_names).to_csv(OD / "pb_subset_counts.csv.gz")
pd.DataFrame(meta).to_csv(OD / "pb_subset_samples.csv", index=False)
print(pd.DataFrame(meta).to_string(index=False))
print("wrote", OD)
