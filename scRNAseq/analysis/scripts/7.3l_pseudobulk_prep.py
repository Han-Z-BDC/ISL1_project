"""Step 7.3l (prep): replicate pseudobulk counts for the D15 KO vs control tests (after 7.3a, 7.3d; tests in
7.3l_pseudobulk_edger.R).

Decisions (user, 2026-11-21): one simple, comprehensive KO vs control test on all genes, replicate-aware (pseudobulk,
edgeR), in three comparisons:
  1 curve       per lineage (beta = Slingshot lineage ending in pre_beta2, alpha = pre_alpha, EC = pre_EC; 7.3a PCs
                1-30; a cell is on a lineage at curve weight >= W_MIN): NWIN windows of equal cell numbers (pseudotime
                quantiles of the cells on it, as 7.3g); one sample = one replicate in one window
  2 past stop   beta and alpha: the 7.3f stop window (KO cells' 75th-95th percentile of pseudotime on the lineage) and
                the control path just past it (from the window's end to the median of the control cells beyond it);
                samples = KO replicates at the stop, control replicates at the stop, control replicates past the stop
  3 KO states   KO-diverted states vs the nearest control state, per replicate: KO KO_only vs control ENP_FEV_high,
                KO pre_alpha_KO vs control pre_alpha, KO ENP_PHOX2A+ vs control ENP_PHOX2A+ (control has its own
                ENP_PHOX2A+ cells, so that is the nearest state)
Counts = SoupX counts (summed per sample); genes detected (count > 0) in >= 10 of the cells (as 7.3d / 7.3f / 7.3g);
labels as 7.3f / 7.3g (cellrank_velocity_fates.csv). Samples with < MIN_REP_CELLS cells are written with their cell
number and dropped in the tests.
Writes ${OUT_DIR}/day_analysis/D15/trajectory/pseudobulk/ (reruns overwrite):
  pb_counts.csv.gz   gene x sample summed counts
  pb_samples.csv     sample, comparison, set (lineage or state pair), group, window, genotype, replicate, cells,
                     pt_from, pt_to
Run: sbatch analysis/scripts/7.3l_pseudobulk.sb   (CELLRANK_ENV python; the .sb runs prep, tests, re-clustering)
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

OUT_DIR = Path(os.environ["OUT_DIR"])
DD = OUT_DIR / "day_analysis" / "D15"
SUB, TR = DD / "subcluster_endocrine", DD / "trajectory"
PB = TR / "pseudobulk"
W_MIN, NWIN, Q_LO, Q_HI = 0.5, 8, 0.75, 0.95
LINEAGES = {"beta": "pre_beta2", "alpha": "pre_alpha", "EC": "pre_EC"}
STATES = {"KO_only": ("KO_only", "ENP_FEV_high"), "pre_alpha_KO": ("pre_alpha_KO", "pre_alpha"),
          "ENP_PHOX2A+": ("ENP_PHOX2A+", "ENP_PHOX2A+")}   # set: (KO label, control label)


def main():
    PB.mkdir(exist_ok=True)
    a = ad.read_h5ad(SUB / "D15_endocrine.h5ad")
    lab = pd.read_csv(TR / "cellrank_velocity_fates.csv", index_col=0, keep_default_na=False)
    a = a[lab.index].copy()
    C = a.layers["soupx_counts"]
    C = C.tocsr() if sparse.issparse(C) else sparse.csr_matrix(C)
    keep = np.asarray((C > 0).sum(0)).ravel() >= 10
    C, genes = C[:, keep], a.var_names[keep]
    label = lab.label.to_numpy()
    geno = a.obs.genotype.astype(str).to_numpy()
    rep = (a.obs.genotype.astype(str) + "_" + a.obs.replicate.astype(str)).to_numpy()
    reps = sorted(set(rep))
    print(a.n_obs, "cells,", len(genes), "genes; replicates:", pd.Series(rep).value_counts().to_dict())
    sl = pd.read_csv(TR / "slingshot_pc30_cells.csv", index_col=0).reindex(a.obs_names)
    ln = pd.read_csv(TR / "slingshot_pc30_lineages.csv")
    end = {r.lineage: r.clusters.split(" -> ")[-1] for r in ln.itertuples()}

    cols, rows = [], []

    def add(mask, **info):
        for r_ in reps:
            m = mask & (rep == r_)
            if not m.any() and info.get("comparison") != "curve":
                continue
            name = "|".join(str(info[k]) for k in ["comparison", "set", "group"]) + f"|{r_}"
            cols.append(np.asarray(C[m].sum(0)).ravel())
            rows.append(dict(sample=name, **info, genotype=r_.split("_")[0], replicate=r_, cells=int(m.sum())))

    for lin, fate in LINEAGES.items():
        L = [k for k, v in end.items() if v == fate]
        assert len(L) == 1, (fate, L)
        L = L[0]
        pt, w = sl[f"pt_{L}"].to_numpy(), sl[f"w_{L}"].fillna(0).to_numpy()
        on = (w >= W_MIN) & ~np.isnan(pt)
        # 1 curve: equal-count windows (as 7.3g)
        edges = np.quantile(pt[on], np.linspace(0, 1, NWIN + 1))
        for i in range(NWIN):
            m = on & (pt >= edges[i]) & ((pt <= edges[i + 1]) if i == NWIN - 1 else (pt < edges[i + 1]))
            for g in ["ctr", "KO"]:
                add(m & (geno == g), comparison="curve", set=lin, group=f"{g}_w{i + 1}", window=i + 1,
                    pt_from=round(edges[i], 2), pt_to=round(edges[i + 1], 2))
        # 2 past the stop (as 7.3f)
        if lin in ("beta", "alpha"):
            lo_, hi_ = np.quantile(pt[on & (geno == "KO")], [Q_LO, Q_HI])
            beyond = on & (geno == "ctr") & (pt > hi_)
            mid = np.median(pt[beyond])
            stop = on & (pt >= lo_) & (pt <= hi_)
            for grp, m, a_, b_ in [("KO_stop", stop & (geno == "KO"), lo_, hi_),
                                    ("ctr_stop", stop & (geno == "ctr"), lo_, hi_),
                                    ("ctr_after", beyond & (pt <= mid), hi_, mid)]:
                add(m, comparison="after_stop", set=lin, group=grp, window=np.nan, pt_from=round(a_, 2),
                    pt_to=round(b_, 2))
            print(f"{lin}: stop window {lo_:.1f}-{hi_:.1f}, control past the stop {hi_:.1f}-{mid:.1f}")
    # 3 KO-diverted states
    for s, (kl, cl) in STATES.items():
        add((geno == "KO") & (label == kl), comparison="ko_state", set=s, group=f"KO_{kl}", window=np.nan,
            pt_from=np.nan, pt_to=np.nan)
        add((geno == "ctr") & (label == cl), comparison="ko_state", set=s, group=f"ctr_{cl}", window=np.nan,
            pt_from=np.nan, pt_to=np.nan)

    st = pd.DataFrame(rows)
    st.to_csv(PB / "pb_samples.csv", index=False)
    cnt = pd.DataFrame(np.vstack(cols).T.astype(np.int64), index=genes, columns=st["sample"])
    cnt.index.name = "gene"
    cnt.to_csv(PB / "pb_counts.csv.gz")
    piv = st.pivot_table(index=["comparison", "set", "group"], columns="replicate", values="cells", aggfunc="sum")
    print(piv.fillna(0).astype(int).to_string())
    print("wrote", PB / "pb_counts.csv.gz", f"({cnt.shape[0]} genes x {cnt.shape[1]} samples)", "and pb_samples.csv")


if __name__ == "__main__":
    main()
