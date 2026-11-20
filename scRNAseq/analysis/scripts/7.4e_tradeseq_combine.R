# Step 7.4e (combine): join the tradeSeq chunks of the D15 + D22 trajectory (after all 7.4e_tradeseq_fit.R tasks).
#
# FDR (BH) over all genes per lineage / pair, as the D15 run (7.3e); then for the genes with FDR < 0.05 in any
# condition test: smoothed curves (50 points per lineage x genotype, from each chunk's models) and modules by the shape
# of the control curve per lineage (scaled log1p, k-means k = 6, seed 1234). predictSmooth (tradeSeq
# .getPredictRangeDf) fills the stage covariate with the first cell's value; the fit step checks that cell is a D15
# cell, so the curves are drawn at the D15 level (stage D22 indicator = 0). Curve shapes, not levels, define modules.
# Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/trajectory/ (reruns overwrite):
#   tradeseq_association.csv, tradeseq_condition.csv, tradeseq_earlyDE.csv  (gene, test, which, stat, df, p, FDR)
#   tradeseq_smooth.csv, tradeseq_modules.csv
# With INT_SUBSET=<name>: the 7.2 subset folder (see 7.2a).
# Usage: Rscript 7.4e_tradeseq_combine.R <N>
suppressPackageStartupMessages({ library(tradeSeq); library(SingleCellExperiment) })
N <- as.integer(commandArgs(TRUE)[1]); stopifnot(!is.na(N))
SEED <- 1234
SUBSET <- Sys.getenv("INT_SUBSET")   # 7.2 subset run (see 7.2a); empty = the full integration
int_dir <- file.path(Sys.getenv("OUT_DIR"), "day_analysis",
                     paste0("integration_D15_D22_endocrine", if (nzchar(SUBSET)) paste0("_", SUBSET) else ""))
tr <- file.path(int_dir, "trajectory")
ch <- file.path(tr, "tradeseq_chunks")
unlink(list.files(tr, pattern = "^tradeseq_.*\\.csv$", full.names = TRUE))

join <- function(what) {
  fs <- file.path(ch, sprintf("chunk%d_%s.csv", seq_len(N) - 1, what))
  if (what != "earlyDE") stopifnot(all(file.exists(fs)))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) return(NULL)
  d <- do.call(rbind, lapply(fs, read.csv, check.names = FALSE))
  d$FDR <- ave(d$pvalue, d$which, FUN = function(p) p.adjust(p, "fdr"))
  d <- d[order(d$which, d$pvalue), ]
  write.csv(d, file.path(tr, sprintf("tradeseq_%s.csv", what)), row.names = FALSE)
  print(aggregate(FDR ~ test + which, d, function(q) sum(q < 0.05, na.rm = TRUE)))
  d
}
as_res <- join("association"); co_res <- join("condition"); ed_res <- join("earlyDE")
cat("genes tested:", length(unique(as_res$gene)), "\n")

sig <- unique(co_res$gene[!is.na(co_res$FDR) & co_res$FDR < 0.05])
cat("genes with FDR < 0.05 in a condition test:", length(sig), "\n")
if (length(sig)) {
  sm <- do.call(rbind, lapply(seq_len(N) - 1, function(k) {
    fit <- readRDS(file.path(ch, sprintf("chunk%d_models.rds", k)))
    g <- intersect(sig, rownames(fit))
    if (!length(g)) return(NULL)
    tryCatch(predictSmooth(fit, gene = g, nPoints = 50, tidy = TRUE),
             error = function(e) { cat("predictSmooth chunk", k, "failed:", conditionMessage(e), "\n"); NULL })
  }))
  if (!is.null(sm) && nrow(sm)) {
    write.csv(sm, file.path(tr, "tradeseq_smooth.csv"), row.names = FALSE)
    mods <- list()
    for (k in unique(co_res$which)) {
      if (k == "all") next
      g <- co_res$gene[co_res$which == k & !is.na(co_res$FDR) & co_res$FDR < 0.05]
      li <- as.integer(sub("^lineage", "", k))
      if (length(g) < 12 || is.na(li)) next
      s <- sm[sm$gene %in% g & sm$lineage == li & sm$condition == "ctr", ]
      y <- tapply(s$yhat, list(s$gene, s$time), mean)
      y <- y[, order(as.numeric(colnames(y)))]
      z <- t(scale(t(log1p(y)))); z[is.na(z)] <- 0
      set.seed(SEED)
      km <- kmeans(z, centers = min(6, nrow(z) - 1), nstart = 20)
      mods[[k]] <- data.frame(gene = rownames(z), lineage = k, module = km$cluster)
    }
    if (length(mods)) write.csv(do.call(rbind, mods), file.path(tr, "tradeseq_modules.csv"), row.names = FALSE)
  }
}
cat("wrote", tr, "\n")
