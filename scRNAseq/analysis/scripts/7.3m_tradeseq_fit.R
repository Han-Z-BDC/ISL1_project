# Step 7.3m (fit): full D15 tradeSeq, one gene chunk per array task, one core, serial (after 7.3m_tradeseq_prep.R).
#
# Per chunk: fitGAM (6 knots, conditions = genotype, offset, SoupX counts) on the cells of the beta / alpha / EC
# lineages, then gene by gene (FDR in the combine step):
#   association     associationTest per lineage and genotype (gene changes along the lineage)
#   condition       conditionTest per lineage, l2fc = 0 (no fold-change threshold; effect sizes in the combine step):
#                     shared  only the knots inside the KO range of that lineage (knot 1 to the last knot at or below
#                             the KO cells' 95th percentile of pseudotime; KO curves past the stop are extrapolated, so
#                             they are left out); a lineage KO cells cover to the end uses all knots
#                     whole   every knot (as 7.3e / 7.3k, for comparison)
#   earlyDE         earlyDETest knots 1-3, pairwise between lineages
# Reuses a chunk's saved models if they exist (a rerun after a failed test step does not refit).
# Usage: Rscript 7.3m_tradeseq_fit.R <chunk>
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/tradeseq_full/: chunk{k}_models.rds, chunk{k}_association.csv,
#   chunk{k}_condition.csv (which = lineage, range = shared / whole), chunk{k}_earlyDE.csv
suppressPackageStartupMessages({
  library(Matrix); library(tradeSeq); library(SingleCellExperiment); library(BiocParallel)
})
K <- as.integer(commandArgs(TRUE)[1]); stopifnot(!is.na(K))
SEED <- 1234; set.seed(SEED)
# eigenvalue threshold of the Wald test: the default 0.01 left 1,000-2,000 genes per lineage NA in the shared-range
# test (e.g. FEV); 1e-4 gives every gene a result with the full df (checked 2026-11-21 on FEV's chunk)
EIG <- 1e-4
ts <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
fo <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "D15", "trajectory", "tradeseq_full")
inp <- readRDS(file.path(fo, "inputs.rds"))
m <- readRDS(file.path(fo, sprintf("chunk%d_counts.rds", K)))
kr <- read.csv(file.path(fo, "ko_range.csv"))
mf <- file.path(fo, sprintf("chunk%d_models.rds", K))
unlink(file.path(fo, sprintf("chunk%d_%s.csv", K, c("association", "condition", "earlyDE"))))
ts("chunk", K, ":", nrow(m), "genes,", ncol(m), "cells")

t0 <- Sys.time()
if (file.exists(mf)) {
  fit <- readRDS(mf); stopifnot(identical(rownames(fit), rownames(m)))
  ts("reused saved models")
} else {
  fit <- fitGAM(counts = m, pseudotime = inp$pt, cellWeights = inp$w, conditions = inp$cond, nknots = 6,
                offset = inp$offset, parallel = FALSE, verbose = FALSE)
  saveRDS(fit, mf)
  ts("fitGAM:", format(Sys.time() - t0))
}
ok <- !is.na(rowData(fit)$tradeSeq$beta[, 1])
ts("converged:", sum(ok), "of", length(ok))

long <- function(res, what) {   # wide tradeSeq output -> long (gene, which, stat, df, p)
  st <- grep("^waldStat_", names(res), value = TRUE)
  do.call(rbind, lapply(st, function(s) {
    k <- sub("^waldStat_", "", s)
    data.frame(gene = rownames(res), test = what, which = k, waldStat = res[[s]], df = res[[paste0("df_", k)]],
               pvalue = res[[paste0("pvalue_", k)]])
  }))
}
write.csv(long(associationTest(fit, lineages = TRUE), "association"),
          file.path(fo, sprintf("chunk%d_association.csv", K)), row.names = FALSE)

kn <- metadata(fit)$tradeSeq$knots
co <- list(cbind(long(conditionTest(fit, global = FALSE, lineages = TRUE, l2fc = 0, eigenThresh = EIG), "condition"), range = "whole",
                 knots = paste0("1-", length(kn))))
for (i in seq_len(nrow(kr))) {             # lineage i of the fit = row i of ko_range (beta, alpha, EC)
  last <- max(which(kn <= kr$ko_q95[i]))
  if (last < 2) next
  r <- long(conditionTest(fit, global = FALSE, lineages = TRUE, knots = c(1, last), l2fc = 0,
                          eigenThresh = EIG), "condition")
  r <- r[r$which == paste0("lineage", i), ]
  co[[length(co) + 1]] <- cbind(r, range = "shared", knots = paste0("1-", last))
}
write.csv(do.call(rbind, co), file.path(fo, sprintf("chunk%d_condition.csv", K)), row.names = FALSE)
ed <- tryCatch(long(earlyDETest(fit, knots = c(1, 3), pairwise = TRUE), "earlyDE"),
               error = function(e) { cat("earlyDETest failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(ed)) write.csv(ed, file.path(fo, sprintf("chunk%d_earlyDE.csv", K)), row.names = FALSE)
ts("chunk", K, "done:", format(Sys.time() - t0), "; knots", paste(round(kn, 1), collapse = " "))
