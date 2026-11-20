# Step 7.3e: genes along the D15 endocrine trajectory, control vs KO (tradeSeq; after 7.3a).
#
# Request (user, 2026-11-19): find genes linked to where the KO trajectory stops (regression + gene modules), as
# candidates for ISL1-controlled steps toward alpha / beta; to be crossed later with CellRank drivers, ISL1 CUT&RUN
# and bulk KO.
#   trajectory  7.3a Slingshot, PCs 1-30, fixed ends (pseudotime + curve weights per lineage)
#   counts      SoupX counts (layers/soupx_counts of D15_merged.h5ad, 7.1) of the 5,994 cells; offset = log of each
#               cell's total counts over all genes
#   genes       all genes detected (count > 0) in at least 10 of the subset's cells (user, 2026-11-19: all genes, not
#               only the variable genes); genes seen in fewer cells cannot be fitted
#   model       fitGAM, 6 knots, conditions = genotype (one smoother per lineage x genotype)
#   tests       associationTest (gene changes along a lineage, per lineage, control and KO together);
#               conditionTest per lineage (control and KO curves differ; fold-change threshold 1.5);
#               earlyDETest knots 1-3 (genes that separate lineages right after their split), pairwise
#   modules     genes significant in conditionTest (FDR < 0.05) on a lineage: smoothed control curve (50 points),
#               scaled, k-means (k = 6, seed 1234) -> modules by curve shape, per lineage
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
#   tradeseq_association.csv, tradeseq_condition.csv, tradeseq_earlyDE.csv  (gene, lineage / pair, stat, p, FDR)
#   tradeseq_smooth.csv     smoothed curves (gene, lineage, genotype, point, pseudotime, fitted counts) of the genes
#                           significant in any condition test
#   tradeseq_modules.csv    gene, lineage, module
#   tradeseq_models.rds     the fitted model (for later plots / windowed tests)
# Run: sbatch analysis/scripts/7.3e_tradeseq.sb
suppressPackageStartupMessages({
  library(rhdf5); library(Matrix); library(tradeSeq); library(BiocParallel)
})
SEED <- 1234; set.seed(SEED)
od <- Sys.getenv("OUT_DIR")
dd <- file.path(od, "day_analysis", "D15")
sub <- file.path(dd, "subcluster_endocrine"); tr <- file.path(dd, "trajectory")
unlink(list.files(tr, pattern = "^tradeseq_", full.names = TRUE))
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
EXTRA <- c("CHGA", "NEUROG3", "NKX2-2", "NEUROD1", "INSM1", "PAX4", "ARX", "IRX2", "PAX6", "ISL1", "MAFB", "MAFA",
           "INS", "GCG", "SST", "GHRL", "PDX1", "NKX6-1", "FEV", "LMX1A", "TPH1", "PHOX2A", "HHEX", "MNX1", "TOX")

# counts of the subset's cells
f <- file.path(dd, "D15_merged.h5ad")
rd <- function(p) h5read(f, p, bit64conversion = "double")
genes <- as.character(rd("var/_index")); cells <- as.character(rd("obs/_index"))
m <- new("dgCMatrix", i = as.integer(rd("layers/soupx_counts/indices")), p = as.integer(rd("layers/soupx_counts/indptr")),
         x = as.numeric(rd("layers/soupx_counts/data")), Dim = c(length(genes), length(cells)))
h5closeAll()
dimnames(m) <- list(genes, cells)
sl <- read.csv(file.path(tr, "slingshot_pc30_cells.csv"), row.names = 1, check.names = FALSE)
stopifnot(all(rownames(sl) %in% cells))
m <- m[, rownames(sl)]
offset <- log(Matrix::colSums(m))
m <- m[Matrix::rowSums(m > 0) >= 10, ]
stopifnot(all(intersect(EXTRA, genes) %in% rownames(m)) || TRUE)
cat("endocrine genes not fitted (< 10 cells):", paste(setdiff(EXTRA, rownames(m)), collapse = " "), "\n")
cat(ncol(m), "cells,", nrow(m), "genes;", NCORE, "cores\n")

pt <- as.matrix(sl[, grep("^pt_", names(sl))]); w <- as.matrix(sl[, grep("^w_", names(sl))])
colnames(pt) <- colnames(w) <- sub("^pt_", "", colnames(pt))
pt[is.na(pt)] <- 0
cond <- factor(sl$genotype, levels = c("ctr", "KO"))
lin <- read.csv(file.path(tr, "slingshot_pc30_lineages.csv"))
print(lin)

bp <- MulticoreParam(NCORE, RNGseed = SEED)
fit <- fitGAM(counts = m, pseudotime = pt, cellWeights = w, conditions = cond, nknots = 6, offset = offset,
              parallel = TRUE, BPPARAM = bp, verbose = TRUE)
saveRDS(fit, file.path(tr, "tradeseq_models.rds"))
ok <- !is.na(rowData(fit)$tradeSeq$beta[, 1])
cat("converged:", sum(ok), "of", length(ok), "\n")

tidy <- function(res, what) {   # wide tradeSeq output -> long (gene, lineage/pair, stat, p, FDR)
  st <- grep("^waldStat_", names(res), value = TRUE)
  out <- do.call(rbind, lapply(st, function(s) {
    k <- sub("^waldStat_", "", s)
    data.frame(gene = rownames(res), test = what, which = k, waldStat = res[[s]], df = res[[paste0("df_", k)]],
               pvalue = res[[paste0("pvalue_", k)]])
  }))
  out$FDR <- ave(out$pvalue, out$which, FUN = function(p) p.adjust(p, "fdr"))
  out[order(out$which, out$pvalue), ]
}
as_res <- tidy(associationTest(fit, lineages = TRUE), "association")
write.csv(as_res, file.path(tr, "tradeseq_association.csv"), row.names = FALSE)
co_res <- tidy(conditionTest(fit, lineages = TRUE, l2fc = log2(1.5)), "condition")
write.csv(co_res, file.path(tr, "tradeseq_condition.csv"), row.names = FALSE)
ed_res <- tryCatch(tidy(earlyDETest(fit, knots = c(1, 3), pairwise = TRUE), "earlyDE"),
                   error = function(e) { cat("earlyDETest failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(ed_res)) write.csv(ed_res, file.path(tr, "tradeseq_earlyDE.csv"), row.names = FALSE)
for (r in list(as_res, co_res, ed_res)) if (!is.null(r))
  print(aggregate(FDR ~ test + which, r, function(q) sum(q < 0.05, na.rm = TRUE)))

# smoothed curves of the condition-significant genes, and modules by control-curve shape per lineage
sig <- unique(co_res$gene[!is.na(co_res$FDR) & co_res$FDR < 0.05])
cat("genes with FDR < 0.05 in a condition test:", length(sig), "\n")
if (length(sig)) {
  sm <- predictSmooth(fit, gene = sig, nPoints = 50, tidy = TRUE)
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
cat("wrote", tr, "\n")
