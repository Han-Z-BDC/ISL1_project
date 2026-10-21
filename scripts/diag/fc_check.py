import random, sys
import pyBigWig
G = "/scratch/alpine/hzhu10@xsede.org/ISL1_cutNrun/results_group"
random.seed(1)
sizes = [l.split() for l in open("/pl/active/ZhuLab_BDC/genome/hg38.canonical.chrom.sizes")]
rand = []
for _ in range(2000):
    c, n = random.choice(sizes[:22]); p = random.randrange(1000, int(n) - 1000); rand.append((c, p))
def mean_at(bw, pts):
    v = [bw.stats(c, p - 100, p + 100, type="mean")[0] or 0 for c, p in pts]
    return sum(v) / len(v)
ref = {}
for s in ["D12_Ctr", "D15_Ctr", "D22_Ctr"]:
    rows = [l.split("\t") for l in open(f"{G}/final_peaks/optimal/{s}.summits.bed")]
    rows.sort(key=lambda r: -float(r[4]))
    ref[s] = [(r[0], int(r[1])) for r in rows[:500]]
print("pooled FC track: mean fold change in 200 bp around the top 500 optimal summits of each Ctr group, and at 2000 random sites")
print("track\t" + "\t".join(f"top500 {k}" for k in ref) + "\trandom")
for g in ["D12_Ctr", "D15_Ctr", "D22_Ctr", "D12_KO", "D15_KO"]:
    bw = pyBigWig.open(f"{G}/tracks/pooled/{g}_pooled.fc.bw")
    print(g + "\t" + "\t".join(f"{mean_at(bw, ref[k]):.1f}" for k in ref) + f"\t{mean_at(bw, rand):.2f}")
