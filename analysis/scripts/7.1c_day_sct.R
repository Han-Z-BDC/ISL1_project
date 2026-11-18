# Step 7.1c (part 1): SCTransform v2 on one day's merged object, for comparison with the 7.1 log1p route.
#
# Decisions (user, 2026-11-18): same cells and same S - G2M scores as 7.1 (damaged NF > 0.6 already removed);
# SCTransform on the merged object as one sample (like the log1p route; no per-sample models, which would act
# as a batch correction); 3,000 variable features; PCA 50. Three variants:
#   SCT    regress S - G2M
#   SCT_A  regress S - G2M and % MT (pct_counts_mt)
#   SCT_B  regress S - G2M; mitochondrial genes (MT-) dropped from the variable features, which are topped
#          back up to 3,000 from the next-ranked genes (by residual variance)
# Neighbours and UMAP are done in 7.1c_day_sct_compare.py with the same settings as 7.1, so only the
# normalization differs. No clustering.
#
# Input:  ${OUT_DIR}/day_analysis/{day}/{day}_merged.h5ad (layers/soupx_counts = SoupX counts, the same counts
#         7.1 normalised; obs S_minus_G2M, pct_counts_mt).
# Writes: ${OUT_DIR}/day_analysis/{day}/sct/ (reruns overwrite this variant's files)
#   {variant}_pca.csv.gz                cells x 50 PCs
#   {variant}_variable_features.csv     gene (original name), rank, residual_variance
#   {variant}_pca_variance.csv          PC, stdev, variance ratio (of the 50)
# Usage: Rscript 7.1c_day_sct.R <D15|D22|D40> <SCT|SCT_A|SCT_B>   (via 7.1c_day_sct.sb)
suppressPackageStartupMessages({
  library(rhdf5); library(Matrix); library(Seurat); library(glmGamPoi)
})
args <- commandArgs(TRUE); day <- args[1]; variant <- args[2]
stopifnot(variant %in% c("SCT", "SCT_A", "SCT_B"))
res <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", day)
sct_dir <- file.path(res, "sct")
dir.create(sct_dir, showWarnings = FALSE)
unlink(list.files(sct_dir, pattern = paste0("^", variant, "_(pca|variable|pca_variance)"), full.names = TRUE))
f <- file.path(res, paste0(day, "_merged.h5ad"))
SEED <- 1234; N_VF <- 3000; N_PC <- 50

rd <- function(p) h5read(f, p, bit64conversion = "double")
genes <- as.character(rd("var/_index")); cells <- as.character(rd("obs/_index"))
L <- "layers/soupx_counts"
# h5ad CSR (cells x genes) == CSC (genes x cells)
m <- new("dgCMatrix", i = as.integer(rd(paste0(L, "/indices"))), p = as.integer(rd(paste0(L, "/indptr"))),
         x = as.numeric(rd(paste0(L, "/data"))), Dim = c(length(genes), length(cells)))
stopifnot(all(m@x == round(m@x)))
seurat_names <- make.unique(gsub("_", "-", genes))      # Seurat's renaming, mapped back below
dimnames(m) <- list(seurat_names, cells)
smg <- as.numeric(rd("obs/S_minus_G2M")); pmt <- as.numeric(rd("obs/pct_counts_mt"))
h5closeAll()
cat(day, variant, ":", ncol(m), "cells,", nrow(m), "genes, nnz", length(m@x), "\n")

so <- CreateSeuratObject(m, min.cells = 0, min.features = 0)
rm(m); invisible(gc())
so$S_minus_G2M <- smg; so$pct_mt <- pmt
regress <- if (variant == "SCT_A") c("S_minus_G2M", "pct_mt") else "S_minus_G2M"
n_fit <- if (variant == "SCT_B") N_VF + 100 else N_VF   # spare ranked genes to replace MT- genes
so <- SCTransform(so, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = regress,
                  variable.features.n = n_fit, seed.use = SEED, verbose = TRUE)
vf <- VariableFeatures(so)
if (variant == "SCT_B") {
  mt <- grepl("^MT-", vf)
  cat("MT- genes dropped from the top", N_VF, ":", sum(mt[seq_len(N_VF)]), "-",
      paste(vf[seq_len(N_VF)][mt[seq_len(N_VF)]], collapse = " "), "\n")
  vf <- head(vf[!mt], N_VF)
  VariableFeatures(so) <- vf
}
stopifnot(length(vf) == N_VF)
so <- RunPCA(so, features = vf, npcs = N_PC, seed.use = SEED, verbose = FALSE)

pca <- Embeddings(so, "pca")
stopifnot(identical(rownames(pca), cells))
gz <- gzfile(file.path(sct_dir, paste0(variant, "_pca.csv.gz")), "w")
write.csv(pca, gz); close(gz)
sd <- Stdev(so, "pca")
write.csv(data.frame(PC = seq_along(sd), stdev = sd, variance_ratio = sd^2 / sum(sd^2)),
          file.path(sct_dir, paste0(variant, "_pca_variance.csv")), row.names = FALSE)
fa <- SCTResults(so[["SCT"]], slot = "feature.attributes")
write.csv(data.frame(gene = genes[match(vf, seurat_names)], rank = seq_along(vf),
                     residual_variance = fa[vf, "residual_variance"]),
          file.path(sct_dir, paste0(variant, "_variable_features.csv")), row.names = FALSE)
cat("variable features:", length(vf), "; wrote", sct_dir, "\n")
