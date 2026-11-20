# Step 7.4e (fit): tradeSeq on the D15 + D22 endocrine trajectory, one chunk of genes per array task (after 7.4a).
#
# Decisions (user, 2026-11-19): the 7.3e design on the 7.2 integration (scVI), with
#   genes       all genes detected (count > 0) in >= 10 of the trajectory cells (as D15)
#   conditions  genotype (control vs KO: one smoother per lineage x genotype)
#   covariate   stage (D15 / D22) as a fixed effect: U = indicator of D22 (tradeSeq's U, no intercept)
# Same model otherwise: SoupX counts (7.2a merged_counts.h5ad), offset = log of each cell's total over all genes,
# 6 knots, Slingshot pseudotime + curve weights (7.4a, scVI latent, fixed D22 ends).
# Genes are split into N chunks (gene i goes to chunk i mod N) so the fits run side by side; each chunk runs the tests
# gene by gene (no FDR here) and keeps its fitted models. 7.4e_tradeseq_combine.R adds the FDR over all genes,
# smoothed curves and modules.
# With INT_SUBSET=<name>: the 7.2 subset folder (see 7.2a).
# Usage: Rscript 7.4e_tradeseq_fit.R <chunk 0..N-1> <N>
# Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/trajectory/tradeseq_chunks/ (reruns overwrite):
#   chunk{k}_models.rds, chunk{k}_association.csv, chunk{k}_condition.csv, chunk{k}_earlyDE.csv (gene, which, stat,
#   df, p)
suppressPackageStartupMessages({
  library(rhdf5); library(Matrix); library(tradeSeq); library(BiocParallel)
})
args <- commandArgs(TRUE); K <- as.integer(args[1]); N <- as.integer(args[2])
stopifnot(!is.na(K), !is.na(N), K >= 0, K < N)
SEED <- 1234; set.seed(SEED)
od <- Sys.getenv("OUT_DIR")
SUBSET <- Sys.getenv("INT_SUBSET")   # 7.2 subset run (see 7.2a); empty = the full integration
int_dir <- file.path(Sys.getenv("OUT_DIR"), "day_analysis",
                     paste0("integration_D15_D22_endocrine", if (nzchar(SUBSET)) paste0("_", SUBSET) else ""))
int <- int_dir
tr <- file.path(int, "trajectory"); ch <- file.path(tr, "tradeseq_chunks")
dir.create(ch, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(ch, pattern = sprintf("^chunk%d_", K), full.names = TRUE))
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

f <- file.path(int, "merged_counts.h5ad")
rd <- function(p) h5read(f, p, bit64conversion = "double")
genes <- as.character(rd("var/_index")); cells <- as.character(rd("obs/_index"))
# X is CSR cells x genes == CSC genes x cells
m <- new("dgCMatrix", i = as.integer(rd("X/indices")), p = as.integer(rd("X/indptr")), x = as.numeric(rd("X/data")),
         Dim = c(length(genes), length(cells)))
h5closeAll()
dimnames(m) <- list(make.unique(genes), cells)
sl <- read.csv(file.path(tr, "slingshot_scvi_cells.csv"), row.names = 1, check.names = FALSE)
stopifnot(all(rownames(sl) %in% cells))
m <- m[, rownames(sl)]
stopifnot(all(m@x == round(m@x)))
offset <- log(Matrix::colSums(m))
keep <- which(Matrix::rowSums(m > 0) >= 10)
mine <- keep[(seq_along(keep) - 1) %% N == K]
m <- m[mine, ]
cat("chunk", K, "of", N, ":", nrow(m), "genes of", length(keep), "fitted in total;", ncol(m), "cells;", NCORE, "cores\n")

pt <- as.matrix(sl[, grep("^pt_", names(sl))]); w <- as.matrix(sl[, grep("^w_", names(sl))])
colnames(pt) <- colnames(w) <- sub("^pt_", "", colnames(pt))
pt[is.na(pt)] <- 0
cond <- factor(sl$genotype, levels = c("ctr", "KO"))
stage <- factor(sl$stage, levels = c("D15", "D22"))
U <- model.matrix(~ stage)[, "stageD22", drop = FALSE]     # no intercept (tradeSeq fitGAM, argument U)
# predictSmooth fills U from the first cell (combine step); keep that a D15 cell so curves are at the D15 level
stopifnot(stage[1] == "D15")
print(table(stage, cond))

bp <- MulticoreParam(NCORE, RNGseed = SEED)
t0 <- Sys.time()
fit <- fitGAM(counts = m, pseudotime = pt, cellWeights = w, conditions = cond, U = U, nknots = 6, offset = offset,
              parallel = TRUE, BPPARAM = bp, verbose = TRUE)
cat("fitGAM:", format(Sys.time() - t0), "\n")
saveRDS(fit, file.path(ch, sprintf("chunk%d_models.rds", K)))
ok <- !is.na(rowData(fit)$tradeSeq$beta[, 1])
cat("converged:", sum(ok), "of", length(ok), "\n")

long <- function(res, what) {   # wide tradeSeq output -> long (gene, which, stat, df, p); FDR in the combine step
  st <- grep("^waldStat_", names(res), value = TRUE)
  do.call(rbind, lapply(st, function(s) {
    k <- sub("^waldStat_", "", s)
    data.frame(gene = rownames(res), test = what, which = k, waldStat = res[[s]], df = res[[paste0("df_", k)]],
               pvalue = res[[paste0("pvalue_", k)]])
  }))
}
write.csv(long(associationTest(fit, lineages = TRUE), "association"),
          file.path(ch, sprintf("chunk%d_association.csv", K)), row.names = FALSE)
write.csv(long(conditionTest(fit, lineages = TRUE, l2fc = log2(1.5)), "condition"),
          file.path(ch, sprintf("chunk%d_condition.csv", K)), row.names = FALSE)
ed <- tryCatch(long(earlyDETest(fit, knots = c(1, 3), pairwise = TRUE), "earlyDE"),
               error = function(e) { cat("earlyDETest failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(ed)) write.csv(ed, file.path(ch, sprintf("chunk%d_earlyDE.csv", K)), row.names = FALSE)
cat("chunk", K, "done:", format(Sys.time() - t0), "\n")
