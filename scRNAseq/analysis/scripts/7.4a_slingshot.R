# Step 7.4a: Slingshot on the D15 + D22 endocrine integration (after 7.4 prep), in the scVI latent space.
#
# Decisions (user, 2026-11-19): same trajectory set as D15 (7.3), on the 7.2 integration with scVI;
#   space     scVI latent, all 30 dimensions (the integration chosen in 7.2; never a UMAP)
#   clusters  the 7.2 labels: 11 D15 states ("D15 <name>") + 7 D22 labels
#   start     D15 ENP1
#   ends      fixed to the D22 mature labels SC_beta, SC_alpha, SC_delta, SC_EC1, SC_EC2; other clusters
#             (D15 ENP_PHOX2A+, D15 KO_only, D22 GHRL+, D22 ENP, ...) may end lineages freely
#   tests     condiments, conditions = genotype: progressionTest per lineage, differentiationTest (global, pairwise)
# Writes ${OUT_DIR}/day_analysis/integration_D15_D22_endocrine/trajectory/ (reruns overwrite):
#   slingshot_scvi_cells.csv      barcode, label, stage, genotype, pseudotime and weight per lineage
#   slingshot_scvi_lineages.csv   lineage, cluster order
#   slingshot_scvi_curves_2d.csv / _curves_3d.csv   curves embedded in the scVI 2D / 3D UMAP (embedCurves)
#   slingshot_scvi_tests.csv      condiments tests
# With INT_SUBSET=<name>: the 7.2 subset folder (see 7.2a), and the ends are the mature labels above that are in
# the subset (alpha_beta_EC: SC_beta, SC_alpha, SC_EC1, SC_EC2).
# Run: sbatch analysis/scripts/7.4a_slingshot.sb
suppressPackageStartupMessages({ library(slingshot) })
SEED <- 1234; set.seed(SEED)
START <- "D15 ENP1"
END <- c("SC_beta", "SC_alpha", "SC_delta", "SC_EC1", "SC_EC2")
od <- Sys.getenv("OUT_DIR")
SUBSET <- Sys.getenv("INT_SUBSET")   # 7.2 subset run (see 7.2a); empty = the full integration
int_dir <- file.path(Sys.getenv("OUT_DIR"), "day_analysis",
                     paste0("integration_D15_D22_endocrine", if (nzchar(SUBSET)) paste0("_", SUBSET) else ""))
tr <- file.path(int_dir, "trajectory")
unlink(list.files(tr, pattern = "^slingshot_", full.names = TRUE))

ce <- read.csv(file.path(tr, "cells.csv"), row.names = 1, check.names = FALSE)
lat <- as.matrix(read.csv(gzfile(file.path(tr, "scvi_latent.csv.gz")), row.names = 1, check.names = FALSE))
END <- intersect(END, ce$label)     # a subset run keeps only some of the mature labels
cat("ends:", paste(END, collapse = ", "), "\n")
stopifnot(identical(rownames(lat), rownames(ce)), ncol(lat) == 30, START %in% ce$label, length(END) >= 2)
lab <- ce$label; geno <- ce$genotype
u2 <- as.matrix(ce[, c("UMAP1", "UMAP2")]); u3 <- as.matrix(ce[, c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")])
has_cond <- requireNamespace("condiments", quietly = TRUE)
cat(nrow(lat), "cells;", length(unique(lab)), "labels; condiments:", has_cond, "\n")
print(table(ce$stage, geno))

curves_out <- function(emb, f, lineages) {
  cv <- slingCurves(emb, as.df = FALSE)
  stopifnot(length(cv) == length(lineages))
  if (is.null(names(cv))) names(cv) <- lineages     # embedCurves (slingshot 2.14) drops the names (see 7.3a)
  d <- do.call(rbind, lapply(names(cv), function(l) {
    s <- cv[[l]]$s[cv[[l]]$ord, , drop = FALSE]
    data.frame(lineage = l, point = seq_len(nrow(s)), s, check.names = FALSE)
  }))
  stopifnot(!is.null(d), nrow(d) > 0)
  write.csv(d, f, row.names = FALSE)
}

sds <- slingshot(lat, clusterLabels = lab, start.clus = START, end.clus = END, approx_points = 150)
ln <- slingLineages(sds)
cat("\nscVI latent: ", length(ln), " lineages\n", sep = "")
for (l in names(ln)) cat(" ", l, ":", paste(ln[[l]], collapse = " -> "), "\n")
write.csv(data.frame(lineage = names(ln), clusters = sapply(ln, paste, collapse = " -> ")),
          file.path(tr, "slingshot_scvi_lineages.csv"), row.names = FALSE)
pt <- slingPseudotime(sds); w <- slingCurveWeights(sds)
colnames(pt) <- paste0("pt_", colnames(pt)); colnames(w) <- paste0("w_", colnames(w))
write.csv(data.frame(label = lab, stage = ce$stage, genotype = geno, pt, w, check.names = FALSE),
          file.path(tr, "slingshot_scvi_cells.csv"))
curves_out(embedCurves(sds, u2), file.path(tr, "slingshot_scvi_curves_2d.csv"), names(ln))
curves_out(embedCurves(sds, u3), file.path(tr, "slingshot_scvi_curves_3d.csv"), names(ln))

if (has_cond) {
  tests <- list()
  pr <- tryCatch(condiments::progressionTest(sds, conditions = geno, lineages = TRUE),
                 error = function(e) { cat("progressionTest failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(pr)) tests$progression <- data.frame(test = "progression", as.data.frame(pr))
  if (length(ln) > 1) {
    df <- tryCatch(condiments::differentiationTest(sds, conditions = geno, global = TRUE, pairwise = TRUE),
                   error = function(e) { cat("differentiationTest failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(df)) tests$differentiation <- data.frame(test = "differentiation", as.data.frame(df))
  }
  if (length(tests)) {
    cols <- unique(unlist(lapply(tests, names)))
    tt <- do.call(rbind, lapply(tests, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
    print(tt); write.csv(tt, file.path(tr, "slingshot_scvi_tests.csv"), row.names = FALSE)
  }
}
cat("wrote", tr, "\n")
