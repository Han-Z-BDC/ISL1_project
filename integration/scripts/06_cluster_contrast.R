#!/usr/bin/env Rscript
# Step 06: ISL1 binding at alpha-leaning against beta-leaning sites, per day and over time, libraries as replicates.
#   Site sets (sites linked by GREAT, step 01; a site in both sets of a comparison counts in both)
#     cluster4 / cluster2   sites of the KO-effect clusters (step 03 site_groups.tsv): cluster 4 leans alpha
#                           (GCG, ARX, IRX2, TTR, LOXL4, ...), cluster 2 leans beta (INS, MAFA, PCSK1, ABCC8, ...)
#     alpha / beta          sites of all alpha / beta marker genes (markers.tsv)
#     alpha_in_c4 / beta_in_c2  sites of the alpha markers in cluster 4 / the beta markers in cluster 2
#     all                   all sites
#   Per library and set, from the step 05 normalized signal (DiffBind reads / background size factor, steps 00 and 05):
#     total   sum over the set's sites (depends on how many sites a set has)
#     median  median signal per site (does not depend on the number of sites)
#   Outcome per library: log2(set A / set B) for each comparison and metric.
#   Model: lm(outcome ~ batch + day), 9 libraries (batch 3 vs earlier; 5 residual df).
#     level_<day>   the outcome on that day at the libraries' batch mix (1 of 3 libraries per day is batch 3):
#                   whether set A binds more than set B on that day. Between different site sets, levels also
#                   reflect site number and intrinsic site strength, so the changes are the cleaner read-out.
#     change_<a>_vs_<b>  difference between days: whether A gains or loses relative to B over time.
# usage: Rscript 06_cluster_contrast.R <step 27 isl1_diff dir> <integration results dir> <markers.tsv> <out_dir>
# Output in <out_dir>: library_values.tsv, tests.tsv, set_sizes.tsv

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
cnr <- args[1]; res <- args[2]; mfile <- args[3]; od <- args[4]
dir.create(od, showWarnings = FALSE, recursive = TRUE)

smp <- fread(file.path(cnr, "samples.tsv"))
sig <- fread(file.path(res, "05_binding_signal", "site_signal.tsv"))
stopifnot(all(smp$sample %in% names(sig)))
sg <- fread(file.path(res, "03_binding_change", "site_groups.tsv"))
sgl <- fread(file.path(res, "01_great", "site_gene.tsv"), sep = "\t", quote = "")
gt <- fread(file.path(res, "02_direct_targets", "gene_table.tsv"), sep = "\t", quote = "", na.strings = "NA")
mk <- fread(mfile, skip = "cell_type")
mk[gt, group := i.group, on = c(gene = "gene_name")]
msites <- function(genes) unique(sgl[gene_name %in% genes, site])

sets <- list(
  cluster4 = sg[group == "cluster4", unique(site)], cluster2 = sg[group == "cluster2", unique(site)],
  alpha = msites(mk[cell_type == "alpha", gene]), beta = msites(mk[cell_type == "beta", gene]),
  alpha_in_c4 = msites(mk[cell_type == "alpha" & group == "cluster4", gene]),
  beta_in_c2 = msites(mk[cell_type == "beta" & group == "cluster2", gene]),
  all = sig$site)
fwrite(data.table(set = names(sets), n_sites = lengths(sets)), file.path(od, "set_sizes.tsv"), sep = "\t")
comparisons <- list(cluster4_vs_cluster2 = c("cluster4", "cluster2"), cluster4_vs_all = c("cluster4", "all"),
                    cluster2_vs_all = c("cluster2", "all"), alpha_vs_beta = c("alpha", "beta"),
                    alpha_in_c4_vs_beta_in_c2 = c("alpha_in_c4", "beta_in_c2"))

m <- as.matrix(sig[, smp$sample, with = FALSE]); rownames(m) <- sig$site
summ <- function(s, metric) if (metric == "total") colSums(m[s, , drop = FALSE]) else apply(m[s, , drop = FALSE], 2, median)
lv <- rbindlist(lapply(names(comparisons), function(cn) {
  ab <- comparisons[[cn]]
  rbindlist(lapply(c("total", "median"), function(metric) {
    a <- summ(sets[[ab[1]]], metric); b <- summ(sets[[ab[2]]], metric)
    data.table(comparison = cn, metric = metric, sample = smp$sample, day = smp$day, batch = smp$batch,
               value_a = round(a, 3), value_b = round(b, 3), log2_ratio = round(log2(a / b), 4))
  }))
}))
fwrite(lv, file.path(od, "library_values.tsv"), sep = "\t", quote = FALSE)

# weights on (Intercept, batchbatch3, dayD18, dayD22); batch at the libraries' mix (1/3 batch 3)
W <- list(level_D15 = c(1, 1 / 3, 0, 0), level_D18 = c(1, 1 / 3, 1, 0), level_D22 = c(1, 1 / 3, 0, 1),
          change_D18_vs_D15 = c(0, 0, 1, 0), change_D22_vs_D18 = c(0, 0, -1, 1), change_D22_vs_D15 = c(0, 0, 0, 1))
tests <- lv[, {
  d <- data.frame(y = log2_ratio, day = factor(day, levels = c("D15", "D18", "D22")),
                  batch = factor(batch, levels = c("earlier", "batch3")))
  fit <- lm(y ~ batch + day, data = d)
  stopifnot(identical(names(coef(fit)), c("(Intercept)", "batchbatch3", "dayD18", "dayD22")))
  rbindlist(lapply(names(W), function(tn) {
    w <- W[[tn]]; est <- sum(w * coef(fit)); se <- sqrt(drop(t(w) %*% vcov(fit) %*% w))
    data.table(term = tn, estimate = round(est, 4), se = round(se, 4), df = fit$df.residual,
               p = signif(2 * pt(-abs(est / se), fit$df.residual), 3))
  }))
}, by = .(comparison, metric)]
fwrite(tests, file.path(od, "tests.tsv"), sep = "\t", quote = FALSE)
print(data.table(set = names(sets), n_sites = lengths(sets)))
print(tests)
