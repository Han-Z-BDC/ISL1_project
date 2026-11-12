#!/usr/bin/env Rscript
# Step 06: ISL1 binding over time at the sites of cluster 4 (alpha-leaning: GCG, ARX, IRX2, TTR, LOXL4, ...) against
# cluster 2 (beta-leaning: INS, MAFA, PCSK1, ABCC8, ...), with the libraries as replicates.
#   Per library: total normalized signal (step 05 site_signal.tsv: reads / step 27 background size factor) over the
#   sites linked (GREAT) to each cluster's genes (step 03 site_groups.tsv) and over all sites.
#   Outcomes per library: log2(cluster 4 / cluster 2), log2(cluster 4 / all sites), log2(cluster 2 / all sites).
#   Model: lm(outcome ~ batch + day) on the 9 libraries (batch 3 vs earlier; day D15 / D18 / D22; 5 residual df);
#   contrasts D18 vs D15, D22 vs D18, D22 vs D15 (t test from the coefficient covariance).
#   Sites linked to genes of both clusters count in both totals; their number is reported.
# usage: Rscript 06_cluster_contrast.R <step 27 isl1_diff dir> <integration results dir> <out_dir>
# Output in <out_dir>: library_totals.tsv, contrast_tests.tsv

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
cnr <- args[1]; res <- args[2]; od <- args[3]
dir.create(od, showWarnings = FALSE, recursive = TRUE)

smp <- fread(file.path(cnr, "samples.tsv"))
sig <- fread(file.path(res, "05_binding_signal", "site_signal.tsv"))
stopifnot(all(smp$sample %in% names(sig)))
sg <- fread(file.path(res, "03_binding_change", "site_groups.tsv"))
s4 <- sg[group == "cluster4", site]; s2 <- sg[group == "cluster2", site]
n_shared <- length(intersect(s4, s2))

tot <- function(sites) colSums(as.matrix(sig[site %in% sites, smp$sample, with = FALSE]))
lt <- data.table(sample = smp$sample, day = factor(smp$day, levels = c("D15", "D18", "D22")),
                 batch = factor(smp$batch, levels = c("earlier", "batch3")),
                 total_cluster4 = round(tot(s4), 1), total_cluster2 = round(tot(s2), 1),
                 total_all = round(tot(sig$site), 1))
lt[, `:=`(log2_c4_vs_c2 = log2(total_cluster4 / total_cluster2), log2_c4_vs_all = log2(total_cluster4 / total_all),
          log2_c2_vs_all = log2(total_cluster2 / total_all))]
fwrite(lt[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], file.path(od, "library_totals.tsv"),
       sep = "\t", quote = FALSE)

cts <- list(D18_vs_D15 = c(dayD18 = 1, dayD22 = 0), D22_vs_D18 = c(dayD18 = -1, dayD22 = 1),
            D22_vs_D15 = c(dayD18 = 0, dayD22 = 1))
out <- rbindlist(lapply(c("log2_c4_vs_c2", "log2_c4_vs_all", "log2_c2_vs_all"), function(y) {
  fit <- lm(reformulate(c("batch", "day"), y), data = lt)
  b <- coef(fit)[c("dayD18", "dayD22")]; V <- vcov(fit)[c("dayD18", "dayD22"), c("dayD18", "dayD22")]
  rbindlist(lapply(names(cts), function(ct) {
    w <- cts[[ct]]; est <- sum(w * b); se <- sqrt(drop(t(w) %*% V %*% w))
    tt <- est / se; p <- 2 * pt(-abs(tt), df = fit$df.residual)
    data.table(outcome = y, contrast = ct, estimate = round(est, 4), se = round(se, 4), t = round(tt, 3),
               df = fit$df.residual, p = signif(p, 3))
  }))
}))
out[, n_sites_cluster4 := length(s4)][, n_sites_cluster2 := length(s2)][, n_sites_shared := n_shared]
fwrite(out, file.path(od, "contrast_tests.tsv"), sep = "\t", quote = FALSE)
print(lt)
print(out)
