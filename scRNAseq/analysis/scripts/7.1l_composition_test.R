# Step 7.1l: does the share of each cell label differ between ISL1 KO and control? (after 7.1h for each day)
#
# Request (user, 2026-11-22): per day (D22, D40), compare each label's proportion between genotypes with statistics,
# using the 2 control and 2 KO samples as replicates. Tested within a day only (each day = one differentiation and
# one 10x OCM reaction, so the replicates are wells of the same differentiation, not separate differentiations).
#   propeller (main; Phipson et al. 2022, as speckle::propeller, written out here with limma): per sample the
#     proportion of its cells in each label, arcsine square-root transformed; limma lmFit ~ genotype, eBayes(robust =
#     TRUE); BH FDR across the labels of the day.
#   quasi-binomial GLM (check): per label, cbind(label cells, other cells) ~ genotype over the 4 samples; BH FDR per
#     day. With a genotype at 0 cells the coefficient is not estimable (reported NA).
#   Effect: log2((mean KO proportion + 0.001) / (mean control proportion + 0.001)) (0.1% offset for labels near 0).
# Writes ${OUT_DIR}/day_analysis/composition_genotype/ (reruns overwrite):
#   proportions.csv   day, label, sample, genotype, cells, total, proportion
#   tests.csv         day, label, mean proportion control / KO, log2 ratio, propeller t / p / FDR, GLM log2 OR / p / FDR
# D15_endocrine (user, 2026-11-22): the 7.1d D15 endocrine subset (5,994 cells) with its r0.6 names
#   (day_analysis/D15/subcluster_endocrine/subclusters.csv + study_config/annotation/D15_endocrine_r0.6.tsv); proportions
#   are of the subset's cells per sample, not of all D15 cells.
# Usage: Rscript 7.1l_composition_test.R D15 D22 D40 D15_endocrine
suppressPackageStartupMessages({ library(limma); library(data.table) })
days <- commandArgs(TRUE); stopifnot(length(days) > 0)
da <- file.path(Sys.getenv("OUT_DIR"), "day_analysis")
od <- file.path(da, "composition_genotype"); dir.create(od, showWarnings = FALSE)
OFF <- 0.001

props <- list(); tests <- list()
for (day in days) {
  if (day == "D15_endocrine") {
    sc <- fread(file.path(da, "D15", "subcluster_endocrine", "subclusters.csv")); setnames(sc, 1, "barcode")
    nm <- fread(file.path(Sys.getenv("REPO_DIR"), "study_config", "annotation", "D15_endocrine_r0.6.tsv"))
    l <- sc[, .(barcode, sample = sub("#.*$", "", sub("^[^-]*-", "", barcode)),
                label = nm$label[match(sub_r0p6, nm$subcluster)])]
    l[, genotype := sub("_rep.*$", "", sample)]
    stopifnot(nrow(l) == 5994, !anyNA(l$label), uniqueN(l$sample) == 4)
  } else {
    l <- fread(file.path(da, day, "labels", "labels.csv"))
  }
  l[, geno := ifelse(toupper(genotype) %like% "^KO", "KO", "ctr")]
  cnt <- dcast(l[, .N, by = .(label, sample)], label ~ sample, value.var = "N", fill = 0)
  m <- as.matrix(cnt[, -1]); rownames(m) <- cnt$label
  smp <- unique(l[, .(sample, geno)])[match(colnames(m), sample)]
  tot <- colSums(m); p <- sweep(m, 2, tot, "/")
  props[[day]] <- data.table(day = day, label = rep(rownames(m), ncol(m)), sample = rep(colnames(m), each = nrow(m)),
                             genotype = rep(smp$geno, each = nrow(m)), cells = as.vector(m),
                             total = rep(tot, each = nrow(m)), proportion = as.vector(p))
  # propeller (arcsine square root)
  g <- factor(smp$geno, levels = c("ctr", "KO"))
  des <- model.matrix(~ g)
  fit <- eBayes(lmFit(asin(sqrt(p)), des), robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  mc <- rowMeans(p[, g == "ctr", drop = FALSE]); mk <- rowMeans(p[, g == "KO", drop = FALSE])
  # quasi-binomial GLM per label
  glm_r <- t(sapply(rownames(m), function(k) {
    y <- cbind(m[k, ], tot - m[k, ])
    if (any(tapply(m[k, ], g, sum) == 0)) return(c(NA, NA))
    f <- summary(glm(y ~ g, family = quasibinomial))$coefficients
    c(f[2, 1] / log(2), f[2, 4])
  }))
  tests[[day]] <- data.table(day = day, label = rownames(m), cells = rowSums(m), ctr_mean = mc, KO_mean = mk,
                             log2_ratio = log2((mk + OFF) / (mc + OFF)),
                             direction = ifelse(mk > mc, "higher in KO", "lower in KO"),
                             propeller_t = tt$t, propeller_p = tt$P.Value, propeller_FDR = p.adjust(tt$P.Value, "BH"),
                             glm_log2OR = glm_r[, 1], glm_p = glm_r[, 2],
                             glm_FDR = p.adjust(glm_r[, 2], "BH"))
}
props <- rbindlist(props); tests <- rbindlist(tests)
fwrite(props, file.path(od, "proportions.csv"))
fwrite(tests, file.path(od, "tests.csv"))
print(tests[, .(day, label, cells, ctr = sprintf("%.2f%%", 100 * ctr_mean), KO = sprintf("%.2f%%", 100 * KO_mean),
                log2_ratio = round(log2_ratio, 2), prop_FDR = signif(propeller_FDR, 2), glm_FDR = signif(glm_FDR, 2))],
      row.names = FALSE)
cat("wrote", od, "\n")
