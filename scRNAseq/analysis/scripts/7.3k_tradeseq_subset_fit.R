# Step 7.3k (fit): tradeSeq on the D15 beta / alpha / EC lineages, control vs KO, on a gene subset, one chunk of genes
# per array task (after 7.3a, 7.3j).
#
# Decisions (user, 2026-11-21): the full 7.3e run (23,279 genes x 6 lineages) runs past 2 days; instead
#   genes       7.3j tradeseq_subset_genes.csv: alpha / beta / EC CellRank drivers (control-only runs, either kernel)
#               + the 2,656 KO vs control trajectory genes of 7.3h + endocrine markers
#   lineages    the three of interest only (7.3a PCs 1-30): Lineage1 beta (-> pre_beta2), Lineage2 alpha
#               (-> pre_alpha), Lineage3 EC (-> pre_EC); cells with no weight on any of the three are dropped
#   model       as 7.3e: SoupX counts, offset = log of each cell's total over all genes (all 5,994 cells' totals,
#               computed before any subsetting), 6 knots, conditions = genotype
#   tests       per chunk, gene by gene (FDR in the combine step): associationTest per lineage, conditionTest per
#               lineage (fold-change threshold 1.5), earlyDETest knots 1-3 pairwise
# Genes are split into N chunks (gene i goes to chunk i mod N). Timing test (7.3e_tradeseq_timing, 2026-11-21):
# about 1.5-3 s per gene on one core with 6 lineages.
# Usage: Rscript 7.3k_tradeseq_subset_fit.R <chunk 0..N-1> <N> [reuse]
#   reuse: load that chunk's saved models (chunk{k}_models.rds) instead of fitting again, then run the tests
#   (2026-11-21: the first run saved its models, then stopped at rowData without SingleCellExperiment loaded)
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/tradeseq_subset_chunks/ (reruns overwrite that chunk):
#   chunk{k}_models.rds, chunk{k}_association.csv, chunk{k}_condition.csv, chunk{k}_earlyDE.csv
suppressPackageStartupMessages({
  library(rhdf5); library(Matrix); library(tradeSeq); library(BiocParallel); library(SingleCellExperiment)
})
args <- commandArgs(TRUE); K <- as.integer(args[1]); N <- as.integer(args[2]); REUSE <- isTRUE(args[3] == "reuse")
stopifnot(!is.na(K), !is.na(N), K >= 0, K < N)
SEED <- 1234; set.seed(SEED)
ts <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
dd <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "D15"); tr <- file.path(dd, "trajectory")
ch <- file.path(tr, "tradeseq_subset_chunks")
dir.create(ch, showWarnings = FALSE, recursive = TRUE)
mf <- file.path(ch, sprintf("chunk%d_models.rds", K))
REUSE <- REUSE && file.exists(mf)
old <- list.files(ch, pattern = sprintf("^chunk%d_", K), full.names = TRUE)
unlink(if (REUSE) setdiff(old, mf) else old)
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
LINS <- c(beta = "Lineage1", alpha = "Lineage2", EC = "Lineage3")

f <- file.path(dd, "D15_merged.h5ad")
rd <- function(p) h5read(f, p, bit64conversion = "double")
genes <- as.character(rd("var/_index")); cells <- as.character(rd("obs/_index"))
m <- new("dgCMatrix", i = as.integer(rd("layers/soupx_counts/indices")), p = as.integer(rd("layers/soupx_counts/indptr")),
         x = as.numeric(rd("layers/soupx_counts/data")), Dim = c(length(genes), length(cells)))
h5closeAll(); dimnames(m) <- list(genes, cells)
sl <- read.csv(file.path(tr, "slingshot_pc30_cells.csv"), row.names = 1, check.names = FALSE)
lin <- read.csv(file.path(tr, "slingshot_pc30_lineages.csv"))
stopifnot(grepl("pre_beta2$", lin$clusters[lin$lineage == "Lineage1"]),
          grepl("pre_alpha$", lin$clusters[lin$lineage == "Lineage2"]),
          grepl("pre_EC$", lin$clusters[lin$lineage == "Lineage3"]))
stopifnot(all(rownames(sl) %in% cells))
m <- m[, rownames(sl)]
offset <- log(Matrix::colSums(m))          # all 5,994 cells, all genes (as 7.3e)

pt <- as.matrix(sl[, paste0("pt_", LINS)]); w <- as.matrix(sl[, paste0("w_", LINS)])
colnames(pt) <- colnames(w) <- names(LINS)
pt[is.na(pt)] <- 0; w[is.na(w)] <- 0
on <- rowSums(w) > 0
m <- m[, on]; pt <- pt[on, ]; w <- w[on, ]; offset <- offset[on]
cond <- factor(sl$genotype[on], levels = c("ctr", "KO"))
gs <- read.csv(file.path(tr, "tradeseq_subset_genes.csv"))$gene
miss <- setdiff(gs, rownames(m)); if (length(miss)) cat("genes not in the counts:", paste(miss, collapse = " "), "\n")
gs <- intersect(gs, rownames(m))
gs <- gs[Matrix::rowSums(m[gs, ] > 0) >= 10]
mine <- gs[(seq_along(gs) - 1) %% N == K]
m <- m[mine, , drop = FALSE]
ts("chunk", K, "of", N, ":", nrow(m), "genes of", length(gs), ";", ncol(m), "cells (of", length(on), ") on beta / alpha / EC;",
   NCORE, "cores")
print(table(cond))

t0 <- Sys.time()
if (REUSE) {
  fit <- readRDS(mf); stopifnot(identical(sort(rownames(fit)), sort(rownames(m))))
  ts("reused saved models:", mf)
} else {
  # serial even when Slurm gives several cores: 5-core MulticoreParam took 28.5 min per 212-gene chunk, slower than
  # one core (timing test: 1.5-4.6 s per gene)
  fit <- fitGAM(counts = m, pseudotime = pt, cellWeights = w, conditions = cond, nknots = 6, offset = offset,
                parallel = FALSE, verbose = FALSE)
  ts("fitGAM:", format(Sys.time() - t0))
  saveRDS(fit, mf)
}
ok <- !is.na(rowData(fit)$tradeSeq$beta[, 1])
ts("converged:", sum(ok), "of", length(ok))

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
ts("chunk", K, "done:", format(Sys.time() - t0))
