#!/usr/bin/env Rscript
# Step 05: ISL1 binding signal per site and day, for the sites linked (GREAT) to each KO-effect cluster's genes.
#   Signal     step 27 featureCounts reads per site divided by step 27's background size factors (DESeq2 on 15 kb
#              bins), averaged over the 3 libraries of each day (one batch 3 and two earlier libraries per day, so
#              batch is balanced across days).
#   Groups     clusters 1-6 and no KO effect (reference); sites from step 03 site_groups.tsv (a site linked to genes
#              of two groups counts in both).
#   Tests      paired Wilcoxon signed-rank test between days within each group (sites paired): D18 vs D15,
#              D22 vs D18, D22 vs D15; BH over all tests. Sites are not independent and binding rises at most sites
#              after D15, so these p values say that binding changes at a group's sites, not that it changes more
#              than elsewhere (that is step 03's camera test).
# usage: Rscript 05_binding_signal.R <step 27 isl1_diff dir> <integration results dir> <out_dir>
# Output in <out_dir>: site_signal.tsv (site, normalized signal per library and mean per day),
#   group_signal.tsv (group, site, signal per day), wilcoxon.tsv

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
cnr <- args[1]; res <- args[2]; od <- args[3]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
days <- c("D15", "D18", "D22")

cs <- fread(file.path(cnr, "counts_sites.tsv"))
m <- as.matrix(cs[, -1]); rownames(m) <- cs$name
smp <- fread(file.path(cnr, "samples.tsv"))
stopifnot(identical(colnames(m), smp$sample))
nf <- fread(file.path(cnr, "norm_factors.tsv"))[method == "deseq2" & normalization == "background"]
sf <- nf$norm_factor[match(colnames(m), nf$sample)]
stopifnot(!anyNA(sf))
nm <- sweep(m, 2, sf, "/")
sig <- data.table(site = rownames(nm), round(nm, 3))
for (d in days) sig[, (paste0("signal_", d)) := round(rowMeans(nm[, smp$day == d, drop = FALSE]), 3)]
fwrite(sig, file.path(od, "site_signal.tsv"), sep = "\t", quote = FALSE)

lv <- c("no_KO_effect", paste0("cluster", 1:6))
sg <- fread(file.path(res, "03_binding_change", "site_groups.tsv"))[group %in% lv]
gs <- merge(sg, sig[, c("site", paste0("signal_", days)), with = FALSE], by = "site")
gs[, group := factor(group, levels = lv)]
setorder(gs, group, site)
fwrite(gs, file.path(od, "group_signal.tsv"), sep = "\t", quote = FALSE)

pairs <- list(c("D18", "D15"), c("D22", "D18"), c("D22", "D15"))
wt <- rbindlist(lapply(lv, function(l) {
  x <- gs[group == l]
  rbindlist(lapply(pairs, function(p) {
    a <- x[[paste0("signal_", p[1])]]; b <- x[[paste0("signal_", p[2])]]
    data.table(group = l, comparison = paste(p[1], "vs", p[2]), n_sites = nrow(x),
               median_a = round(median(a), 2), median_b = round(median(b), 2),
               median_log2_ratio = round(median(log2((a + 1) / (b + 1))), 3),
               pct_sites_up = round(100 * mean(a > b), 1),
               p = signif(wilcox.test(a, b, paired = TRUE)$p.value, 3))
  }))
}))
wt[, fdr := signif(p.adjust(p, "BH"), 3)]
fwrite(wt, file.path(od, "wilcoxon.tsv"), sep = "\t", quote = FALSE)
print(wt)
