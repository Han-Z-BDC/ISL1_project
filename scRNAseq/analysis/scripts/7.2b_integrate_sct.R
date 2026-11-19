# Step 7.2b: SCT B on the merged D15 + D22 endocrine cells; Harmony and RPCA on stage (after 7.2a).
#
# Decisions (user, 2026-11-19): SCT B as in the per-day work (SCTransform v2, glmGamPoi, S - G2M regressed only,
# 3,000 variable features with MT- genes removed and topped up from the next-ranked genes); batch = stage
# (D15 vs D22); genotype not corrected; no clustering.
#   uncorrected  SCT B on the merged cells as one model -> PCA (50; 30 used downstream)
#   Harmony      harmony::RunHarmony on PCs 1-30 of that model, group.by.vars = stage (default theta 2)
#   RPCA         Seurat's standard SCT + RPCA workflow: SCT B model per stage (same regression), anchor features =
#                the merged model's 3,000 genes found in both stage models, PrepSCTIntegration, PCA per stage, FindIntegrationAnchors
#                (reduction rpca, dims 1-30), IntegrateData (dims 1-30), PCA on the integrated assay
# Input:  ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/merged_counts.h5ad + cells.tsv (7.2a)
# Writes (same folder; reruns overwrite):
#   SCT_B_variable_features.csv   gene (original name), rank, residual_variance (merged model)
#   emb_uncorrected.csv.gz, emb_harmony.csv.gz, emb_rpca.csv.gz   cells x 30 dimensions
# Usage: Rscript 7.2b_integrate_sct.R   (via 7.2b_integrate_sct.sb)
suppressPackageStartupMessages({
  library(rhdf5); library(Matrix); library(Seurat); library(glmGamPoi); library(harmony)
})
options(future.globals.maxSize = 20 * 1024^3)
out <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "integration_D15_D22_endocrine")
f <- file.path(out, "merged_counts.h5ad")
SEED <- 1234; N_VF <- 3000; N_PC <- 50; N_DIM <- 30

rd <- function(p) h5read(f, p, bit64conversion = "double")
genes <- as.character(rd("var/_index")); cells <- as.character(rd("obs/_index"))
m <- new("dgCMatrix", i = as.integer(rd("X/indices")), p = as.integer(rd("X/indptr")),
         x = as.numeric(rd("X/data")), Dim = c(length(genes), length(cells)))
h5closeAll()
stopifnot(all(m@x == round(m@x)))
g <- tabulate(m@i + 1L, nbins = nrow(m)) > 0          # genes detected in no cell carry no information
m <- m[g, ]; genes <- genes[g]
seurat_names <- make.unique(gsub("_", "-", genes))
dimnames(m) <- list(seurat_names, cells)
obs <- read.delim(file.path(out, "cells.tsv"), row.names = 1, check.names = FALSE)
stopifnot(identical(rownames(obs), cells))
cat("merged:", ncol(m), "cells,", nrow(m), "genes\n"); print(table(obs$stage, obs$genotype))

write_emb <- function(x, name) {
  x <- x[, seq_len(N_DIM)]
  stopifnot(identical(rownames(x), cells))
  gz <- gzfile(file.path(out, paste0("emb_", name, ".csv.gz")), "w"); write.csv(x, gz); close(gz)
  cat("wrote emb_", name, "\n", sep = "")
}
sct_b <- function(so) {   # SCT B: S - G2M regressed; MT- genes removed from the variable features
  so <- SCTransform(so, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "S_minus_G2M",
                    variable.features.n = N_VF + 100, seed.use = SEED, verbose = FALSE)
  vf <- VariableFeatures(so); vf <- head(vf[!grepl("^MT-", vf)], N_VF)
  VariableFeatures(so) <- vf
  so
}

# ---- uncorrected merged SCT B model + Harmony
so <- CreateSeuratObject(m, meta.data = obs, min.cells = 0, min.features = 0)
so <- sct_b(so)
vf <- VariableFeatures(so); stopifnot(length(vf) == N_VF, !any(grepl("^MT-", vf)))
fa <- SCTResults(so[["SCT"]], slot = "feature.attributes")
write.csv(data.frame(gene = genes[match(vf, seurat_names)], rank = seq_along(vf),
                     residual_variance = fa[vf, "residual_variance"]),
          file.path(out, "SCT_B_variable_features.csv"), row.names = FALSE)
so <- RunPCA(so, features = vf, npcs = N_PC, seed.use = SEED, verbose = FALSE)
write_emb(Embeddings(so, "pca"), "uncorrected")
set.seed(SEED)
so <- RunHarmony(so, group.by.vars = "stage", reduction.use = "pca", dims.use = seq_len(N_DIM),
                 reduction.save = "harmony", verbose = FALSE)
write_emb(Embeddings(so, "harmony"), "harmony")
rm(so); invisible(gc())

# ---- RPCA: one SCT B model per stage, anchors on the merged model's genes
so <- CreateSeuratObject(m, meta.data = obs, min.cells = 0, min.features = 0)
rm(m); invisible(gc())
lst <- SplitObject(so, split.by = "stage"); rm(so); invisible(gc())
lst <- lapply(lst, function(x) SCTransform(x, vst.flavor = "v2", method = "glmGamPoi",
                                           vars.to.regress = "S_minus_G2M", variable.features.n = N_VF + 100,
                                           seed.use = SEED, verbose = FALSE))
# a gene enters a per-stage model only if detected in >= 5 of that stage's cells (SCTransform min_cells);
# anchor features = the merged model's genes that are in both stage models
in_all <- Reduce(intersect, lapply(lst, function(x) rownames(SCTResults(x[["SCT"]], slot = "feature.attributes"))))
vf_a <- vf[vf %in% in_all]
cat("anchor features:", length(vf_a), "of", length(vf), "; not in both stage models:",
    paste(genes[match(setdiff(vf, vf_a), seurat_names)], collapse = " "), "\n")
lst <- PrepSCTIntegration(lst, anchor.features = vf_a, verbose = FALSE)
lst <- lapply(lst, function(x) RunPCA(x, features = vf_a, npcs = N_PC, seed.use = SEED, verbose = FALSE))
anc <- FindIntegrationAnchors(lst, normalization.method = "SCT", anchor.features = vf_a, reduction = "rpca",
                              dims = seq_len(N_DIM), verbose = FALSE)
cat("RPCA anchors:", nrow(anc@anchors), "\n")
ig <- IntegrateData(anc, normalization.method = "SCT", dims = seq_len(N_DIM), verbose = FALSE)
DefaultAssay(ig) <- "integrated"
ig <- RunPCA(ig, npcs = N_PC, seed.use = SEED, verbose = FALSE)
e <- Embeddings(ig, "pca")[cells, ]
write_emb(e, "rpca")
cat("done\n")
