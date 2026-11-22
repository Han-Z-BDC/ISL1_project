# Step 7.3m (prep): inputs for the full D15 tradeSeq on the beta / alpha / EC lineages, split into gene chunks small
# enough for one acpu core (3.7 GB) per fit task (after 7.3a).
#
# Decisions (user, 2026-11-21): rerun the full tradeSeq (7.3e was cancelled: parallel fitting was the slow part, and it
# would have stopped at rowData before the tests) as
#   genes       all genes detected (count > 0) in >= 10 of the cells used (no pre-selection)
#   lineages    beta (Lineage1 -> pre_beta2), alpha (Lineage2 -> pre_alpha), EC (Lineage3 -> pre_EC) of 7.3a (PCs
#               1-30); cells with no curve weight on any of the three dropped (as 7.3k)
#   model       as 7.3e / 7.3k: SoupX counts, offset = log of each cell's total over all genes (all 5,994 cells, before
#               any subsetting), 6 knots, conditions = genotype
#   chunks      CHUNK genes per fit task, one core each, fitted serially (timing test: 1.5-4.6 s per gene on one core)
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/tradeseq_full/ (reruns overwrite):
#   inputs.rds        pseudotime, cell weights, conditions, offset (cells on the three lineages)
#   ko_range.csv      per lineage: KO cells' 95th percentile of pseudotime on the lineage (curve weight >= 0.5); the
#                     KO vs control test is limited to the knots inside it (combine / fit step)
#   chunk{k}_counts.rds  genes x cells counts of chunk k (dgCMatrix)
#   chunks.txt        number of chunks (read by the submit script)
# Usage: Rscript 7.3m_tradeseq_prep.R
suppressPackageStartupMessages({ library(rhdf5); library(Matrix) })
CHUNK <- 200
dd <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "D15"); tr <- file.path(dd, "trajectory")
fo <- file.path(tr, "tradeseq_full"); dir.create(fo, showWarnings = FALSE)
unlink(list.files(fo, pattern = "^(chunk\\d+_counts\\.rds|inputs\\.rds|chunks\\.txt|ko_range\\.csv)$", full.names = TRUE))
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
m <- m[, rownames(sl)]
offset <- log(Matrix::colSums(m))
pt <- as.matrix(sl[, paste0("pt_", LINS)]); w <- as.matrix(sl[, paste0("w_", LINS)])
colnames(pt) <- colnames(w) <- names(LINS)
pt[is.na(pt)] <- 0; w[is.na(w)] <- 0
on <- rowSums(w) > 0
m <- m[, on]; pt <- pt[on, ]; w <- w[on, ]; offset <- offset[on]
cond <- factor(sl$genotype[on], levels = c("ctr", "KO"))
m <- m[Matrix::rowSums(m > 0) >= 10, ]
kr <- data.frame(lineage = names(LINS), slingshot = LINS,
                 ko_q95 = sapply(names(LINS), function(l) quantile(pt[cond == "KO" & w[, l] >= 0.5, l], 0.95)),
                 ctr_max = sapply(names(LINS), function(l) max(pt[cond == "ctr" & w[, l] >= 0.5, l])))
write.csv(kr, file.path(fo, "ko_range.csv"), row.names = FALSE)
print(kr)
saveRDS(list(pt = pt, w = w, cond = cond, offset = offset), file.path(fo, "inputs.rds"))
N <- ceiling(nrow(m) / CHUNK)
for (k in seq_len(N) - 1) {
  idx <- seq(k * CHUNK + 1, min((k + 1) * CHUNK, nrow(m)))
  saveRDS(m[idx, , drop = FALSE], file.path(fo, sprintf("chunk%d_counts.rds", k)))
}
writeLines(as.character(N), file.path(fo, "chunks.txt"))
cat(ncol(m), "cells,", nrow(m), "genes,", N, "chunks of", CHUNK, "\n")
print(table(cond))
