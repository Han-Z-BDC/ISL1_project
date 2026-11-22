# TEST (user, 2026-11-21): how long does tradeSeq fitGAM take per gene on the D15 trajectory (same inputs as 7.3e)?
#   2 genes (a 1-gene matrix crashes fitGAM: "invalid row.names length"), 10 genes, 50 genes serially; then 200 genes on NCORE cores (is the parallel run really parallel?).
# Writes nothing; prints timings. Run: sbatch analysis/scripts/7.3e_tradeseq_timing.sb
suppressPackageStartupMessages({
  library(rhdf5); library(Matrix); library(tradeSeq); library(BiocParallel)
})
SEED <- 1234; set.seed(SEED)
ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
od <- Sys.getenv("OUT_DIR")
dd <- file.path(od, "day_analysis", "D15"); tr <- file.path(dd, "trajectory")
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

f <- file.path(dd, "D15_merged.h5ad")
rd <- function(p) h5read(f, p, bit64conversion = "double")
genes <- as.character(rd("var/_index")); cells <- as.character(rd("obs/_index"))
m <- new("dgCMatrix", i = as.integer(rd("layers/soupx_counts/indices")), p = as.integer(rd("layers/soupx_counts/indptr")),
         x = as.numeric(rd("layers/soupx_counts/data")), Dim = c(length(genes), length(cells)))
h5closeAll(); dimnames(m) <- list(genes, cells)
sl <- read.csv(file.path(tr, "slingshot_pc30_cells.csv"), row.names = 1, check.names = FALSE)
m <- m[, rownames(sl)]
offset <- log(Matrix::colSums(m))
m <- m[Matrix::rowSums(m > 0) >= 10, ]
pt <- as.matrix(sl[, grep("^pt_", names(sl))]); w <- as.matrix(sl[, grep("^w_", names(sl))])
colnames(pt) <- colnames(w) <- sub("^pt_", "", colnames(pt)); pt[is.na(pt)] <- 0
cond <- factor(sl$genotype, levels = c("ctr", "KO"))
ts("data ready:", ncol(m), "cells,", nrow(m), "genes,", ncol(pt), "lineages;", NCORE, "cores")

fitn <- function(g, par = FALSE) {
  t0 <- Sys.time()
  fit <- fitGAM(counts = m[g, , drop = FALSE], pseudotime = pt, cellWeights = w, conditions = cond, nknots = 6,
                offset = offset, parallel = par, BPPARAM = if (par) MulticoreParam(NCORE, RNGseed = SEED) else bpparam(),
                verbose = FALSE)
  s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ts(sprintf("%4d genes, parallel=%s: %.1f s total, %.2f s per gene", length(g), par, s, s / length(g)))
  invisible(s)
}
set.seed(SEED)
g1 <- c("ISL1", "PDX1")
g10 <- c(g1, sample(setdiff(rownames(m), g1), 8))
g50 <- c(g10, sample(setdiff(rownames(m), g10), 40))
fitn(g1); fitn(g10); fitn(g50)
s200 <- fitn(sample(rownames(m), 200), par = TRUE)
ts(sprintf("extrapolated full run (%d genes) at the parallel rate: %.1f h", nrow(m), s200 / 200 * nrow(m) / 3600))
