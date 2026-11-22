"""Step 7.1m (prep): are KO SC_EC1 / SC_EC2 cells the same cell types as control SC_EC1 / SC_EC2? (after 7.1h)

Request (user, 2026-11-22): beyond the UMAP overlap, check per day (D22, D40) that the EC1 / EC2 cells of the two
genotypes are similar. This script writes the inputs of 7.1m_ec_similarity_edger.R and two checks:
  pseudobulk     SoupX counts summed per sample x label (all labels; the R step drops sample x label groups with
                 < MIN_CELLS cells)
  dotplot data   per label x genotype x gene: % of cells with log1p > 0 and mean log1p (unregressed X), for GENES
  kNN mixing     on the day's SCT B PCs (30, as 7.1b), within one label: genotype-balanced subsamples (the larger
                 genotype down-sampled to the smaller, N_REP seeds); each cell's K nearest neighbours among the
                 subsample; mixing = share of neighbours from the other group / share expected at random
                 (n_other / (n - 1)); 1 = fully mixed, 0 = separate. References with the same procedure:
                   replicates     rep1 vs rep2 within one genotype and label (what "no difference" looks like)
                   EC1 vs EC2     SC_EC1 vs SC_EC2 within one genotype (what "different cell types" looks like)
Writes ${OUT_DIR}/day_analysis/ec_similarity/ (reruns overwrite):
  pb_counts_{day}.csv.gz, pb_samples_{day}.csv, dotplot_{day}.csv, mixing.csv (day, comparison, label, groups,
  cells per group, mixing mean / sd over seeds)
Run: sbatch analysis/scripts/7.1m_ec_similarity.sb D22 D40
"""
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

OUT_DIR = Path(os.environ["OUT_DIR"])
DAYS = sys.argv[1:]
OUT = OUT_DIR / "day_analysis" / "ec_similarity"
N_PCS, K, N_REP, SEED = 30, 30, 20, 1234
EC = ["SC_EC1", "SC_EC2"]
GENES = ["TPH1", "FEV", "LMX1A", "SLC18A1", "DDC", "TAC1", "CHGA", "CHGB", "PLCL1", "RUNX1T1", "AFF3", "NKX6-1",
         "PRKCH", "ISL1", "NEUROD1", "PAX6", "INS", "GCG", "SST", "ARX", "PDX1", "MAFB"]


def mixing(pcs, group, rng):
    """share of the K neighbours from the other group / expected share, on a balanced subsample"""
    g = np.asarray(group)
    a, b = np.where(g == np.unique(g)[0])[0], np.where(g == np.unique(g)[1])[0]
    n = min(len(a), len(b))
    idx = np.concatenate([rng.choice(a, n, replace=False), rng.choice(b, n, replace=False)])
    x, gg = pcs[idx], g[idx]
    nn = NearestNeighbors(n_neighbors=K + 1).fit(x).kneighbors(x, return_distance=False)[:, 1:]
    other = (gg[nn] != gg[:, None]).mean()
    return other / (n / (2 * n - 1)), n


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for day in DAYS:
        dd = OUT_DIR / "day_analysis" / day
        lab = pd.read_csv(dd / "labels" / "labels.csv", index_col=0)
        a = ad.read_h5ad(dd / f"{day}_merged.h5ad")
        a = a[lab.index].copy()
        geno = np.where(lab.genotype.astype(str).str.upper().str.startswith("KO"), "KO", "ctr")
        a.obs["label"] = lab.label.astype(str).values
        a.obs["geno"] = geno
        a.obs["sample"] = lab["sample"].astype(str).values

        # pseudobulk: sample x label
        grp = (a.obs["sample"] + "|" + a.obs["label"]).values
        cats, inv = np.unique(grp, return_inverse=True)
        M = sp.csr_matrix((np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(cats), len(inv)))
        C = a.layers["soupx_counts"]
        pb = M @ C
        pb = pb.toarray() if sp.issparse(pb) else np.asarray(pb)
        pd.DataFrame(np.rint(pb).astype(np.int64).T, index=a.var_names, columns=cats).to_csv(
            OUT / f"pb_counts_{day}.csv.gz")
        smp = pd.DataFrame({"pb": cats})
        smp[["sample", "label"]] = smp.pb.str.split("|", n=1, expand=True)
        smp["genotype"] = np.where(smp["sample"].str.upper().str.startswith("KO"), "KO", "ctr")
        smp["cells"] = np.bincount(inv)
        smp.to_csv(OUT / f"pb_samples_{day}.csv", index=False)
        print(day, "pseudobulk:", len(cats), "sample x label groups")

        # dot plot data
        g = [x for x in GENES if x in a.var_names]
        X = a[:, g].X
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        d = pd.DataFrame(X, columns=g)
        d["label"], d["geno"] = a.obs["label"].values, geno
        pct = (d.groupby(["label", "geno"])[g].agg(lambda v: (v > 0).mean() * 100)).stack().rename("pct")
        mean = d.groupby(["label", "geno"])[g].mean().stack().rename("mean_log1p")
        n = d.groupby(["label", "geno"]).size().rename("cells")
        dp = pd.concat([pct, mean], axis=1).reset_index().rename(columns={"level_2": "gene"}).merge(
            n.reset_index(), on=["label", "geno"])
        dp.to_csv(OUT / f"dotplot_{day}.csv", index=False)

        # kNN mixing on the SCT B PCs
        pcs = pd.read_csv(dd / "sct" / "SCT_B_pca.csv.gz", index_col=0).reindex(a.obs_names).iloc[:, :N_PCS]
        assert pcs.notna().all().all()
        P = pcs.to_numpy()
        obs = a.obs.assign(rep=a.obs["sample"].str.extract(r"(rep\d)")[0].values)
        tests = []
        for l in EC:
            m = (obs.label == l).to_numpy()
            tests.append(("genotype (KO vs ctr)", l, m, obs.geno.to_numpy()))
            for gt in ["ctr", "KO"]:
                mm = m & (obs.geno == gt).to_numpy()
                tests.append((f"replicates ({gt} rep1 vs rep2)", l, mm, obs.rep.to_numpy()))
        for gt in ["ctr", "KO"]:
            mm = obs.label.isin(EC).to_numpy() & (obs.geno == gt).to_numpy()
            tests.append((f"EC1 vs EC2 ({gt})", "SC_EC1 + SC_EC2", mm, obs.label.to_numpy()))
        for comp, l, m, grp_ in tests:
            gsub = grp_[m]
            if len(np.unique(gsub)) < 2:
                continue
            vals = [mixing(P[m], gsub, np.random.default_rng(SEED + s)) for s in range(N_REP)]
            rows.append({"day": day, "comparison": comp, "label": l,
                         "groups": " / ".join(f"{k} {v}" for k, v in pd.Series(gsub).value_counts().items()),
                         "cells_per_group_used": vals[0][1], "mixing_mean": np.mean([v[0] for v in vals]),
                         "mixing_sd": np.std([v[0] for v in vals])})
            print(day, comp, l, round(rows[-1]["mixing_mean"], 3))
    pd.DataFrame(rows).to_csv(OUT / "mixing.csv", index=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
