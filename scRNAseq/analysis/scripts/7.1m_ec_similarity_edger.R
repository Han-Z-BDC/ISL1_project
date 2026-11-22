# Step 7.1m (pseudobulk): KO vs control SC_EC1 / SC_EC2 similarity on replicate pseudobulks (after
# 7.1m_ec_similarity.py). Per day (D22, D40):
#   correlation   pseudobulks (sample x label, >= MIN_CELLS cells), TMM logCPM (prior count 2) on the day's 3,000 SCT B
#                 variable genes, each gene centred on its mean over the day's pseudobulks, Pearson r between all
#                 pseudobulks -> corr_{day}.csv. Cross-day: the EC pseudobulks of D22 vs D40 on the variable genes
#                 shared by both days (each day centred on its own pseudobulks) -> corr_crossday.csv.
#   DE            edgeR QL per test, samples of the test only: filterByExpr, TMM, ~ group, estimateDisp, glmQLFit
#                 (robust), call = FDR < 0.05 and |log2FC| >= 0.5 (as 7.3l). Tests: KO vs ctr within SC_EC1 and within
#                 SC_EC2 (the question); references on the same 2 vs 2 design: SC_EC2 vs SC_EC1 (control; KO),
#                 SC_EC1 vs SC_beta (control) = what a difference between cell types looks like.
# Writes ${OUT_DIR}/day_analysis/ec_similarity/: corr_{day}.csv, corr_crossday.csv, de_{day}_{test}.csv (all genes),
#   de_summary.csv
suppressPackageStartupMessages({ library(edgeR); library(data.table) })
days <- commandArgs(TRUE)
da <- file.path(Sys.getenv("OUT_DIR"), "day_analysis"); od <- file.path(da, "ec_similarity")
MIN_CELLS <- 20; LFC <- 0.5
EC <- c("SC_EC1", "SC_EC2")
summ <- list(); cen <- list()

for (day in days) {
  cnt <- fread(file.path(od, sprintf("pb_counts_%s.csv.gz", day)))
  genes <- cnt[[1]]; m <- as.matrix(cnt[, -1]); rownames(m) <- genes
  s <- fread(file.path(od, sprintf("pb_samples_%s.csv", day)))[match(colnames(m), pb)]
  keep <- s$cells >= MIN_CELLS
  cat(day, ": pseudobulks kept", sum(keep), "of", length(keep), "( dropped:", paste(s$pb[!keep], collapse = "; "), ")\n")
  m <- m[, keep]; s <- s[keep]
  s[, id := sprintf("%s %s %s", label, genotype, sub(".*_", "", sample))]
  colnames(m) <- s$id

  # correlation
  y <- calcNormFactors(DGEList(m))
  lc <- cpm(y, log = TRUE, prior.count = 2)
  hv <- intersect(fread(file.path(da, day, "sct", "SCT_B_variable_features.csv"))$gene, rownames(lc))
  z <- lc[hv, ] - rowMeans(lc[hv, ])
  cen[[day]] <- z[, s$label %in% EC, drop = FALSE]
  fwrite(data.table(pb = colnames(z), cor(z)), file.path(od, sprintf("corr_%s.csv", day)))

  # DE
  pick <- function(lab, gt) s$id[s$label == lab & s$genotype == gt]
  tests <- list(
    "KO_vs_ctr_SC_EC1"         = list(pick("SC_EC1", "ctr"), pick("SC_EC1", "KO")),
    "KO_vs_ctr_SC_EC2"         = list(pick("SC_EC2", "ctr"), pick("SC_EC2", "KO")),
    "ref_EC2_vs_EC1_ctr"       = list(pick("SC_EC1", "ctr"), pick("SC_EC2", "ctr")),
    "ref_EC2_vs_EC1_KO"        = list(pick("SC_EC1", "KO"), pick("SC_EC2", "KO")),
    "ref_EC1_vs_beta_ctr"      = list(pick("SC_beta", "ctr"), pick("SC_EC1", "ctr")))
  for (tn in names(tests)) {
    a <- tests[[tn]][[1]]; b <- tests[[tn]][[2]]
    if (length(a) < 2 || length(b) < 2) { cat(day, tn, ": < 2 pseudobulks per group, skipped\n"); next }
    g <- factor(c(rep("A", length(a)), rep("B", length(b))))
    x <- DGEList(m[, c(a, b)], group = g)
    x <- x[filterByExpr(x, group = g), , keep.lib.sizes = FALSE]
    x <- calcNormFactors(x)
    des <- model.matrix(~ g)
    x <- estimateDisp(x, des)
    f <- glmQLFit(x, des, robust = TRUE)
    r <- topTags(glmQLFTest(f, coef = 2), n = Inf)$table
    r <- data.table(gene = rownames(r), r)
    r[, call := FDR < 0.05 & abs(logFC) >= LFC]
    fwrite(r, file.path(od, sprintf("de_%s_%s.csv", day, tn)))
    summ[[length(summ) + 1]] <- data.table(day = day, test = tn, group_A = paste(a, collapse = "; "),
      group_B = paste(b, collapse = "; "), genes_tested = nrow(r), DE = r[call == TRUE, .N],
      up_in_B = r[call == TRUE & logFC > 0, .N], down_in_B = r[call == TRUE & logFC < 0, .N],
      top_up = paste(head(r[call == TRUE & logFC > 0][order(FDR), gene], 15), collapse = ","),
      top_down = paste(head(r[call == TRUE & logFC < 0][order(FDR), gene], 15), collapse = ","))
    cat(sprintf("%s %-22s tested %5d  DE %5d (up %d / down %d)\n", day, tn, nrow(r), r[call == TRUE, .N],
                r[call == TRUE & logFC > 0, .N], r[call == TRUE & logFC < 0, .N]))
  }
}
fwrite(rbindlist(summ), file.path(od, "de_summary.csv"))

# cross-day EC correlation
if (length(cen) == 2) {
  g <- intersect(rownames(cen[[1]]), rownames(cen[[2]]))
  cc <- cor(cen[[1]][g, ], cen[[2]][g, ])
  rownames(cc) <- paste(names(cen)[1], rownames(cc)); colnames(cc) <- paste(names(cen)[2], colnames(cc))
  fwrite(data.table(pb = rownames(cc), cc), file.path(od, "corr_crossday.csv"))
  cat("cross-day EC correlation on", length(g), "shared variable genes\n"); print(round(cc, 2))
}
cat("wrote", od, "\n")
