# Step 7.3a: Slingshot trajectories of the D15 endocrine cells (after the 7.1d rerun and the r0.6 names).
#
# Decisions (user, 2026-11-19): trajectory on D15 first; Slingshot computed in PCA space (not on a UMAP), shown on the
# 2D and 3D UMAPs; compared with RNA velocity (7.3b-c). Control and KO share one trajectory (condiments approach).
#   input   day_analysis/D15/subcluster_endocrine/: SCT_B_pca.csv.gz (SCT B refit on the 5,994 cells, 50 PCs),
#           subclusters.csv (sub_r0p6, UMAP1/2, UMAP3D_1-3); names from study_config/annotation/D15_endocrine_r0.6.tsv
#   runs    PCs 1-30 (the clustering's PCs) and PCs 1-15 (less noise); clusters = the r0.6 names; start = ENP1;
#           ends not fixed (Slingshot finds them from the cluster tree)
#   tests   condiments (if installed): progressionTest (does KO progress less far along each lineage?) and
#           differentiationTest (do control and KO split between lineages differently?)
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
#   slingshot_pc{n}_cells.csv     barcode, label, genotype, pseudotime and weight per lineage
#   slingshot_pc{n}_lineages.csv  lineage, cluster order
#   slingshot_pc{n}_curves_2d.csv / _curves_3d.csv   curves embedded in the 2D / 3D UMAP (embedCurves)
#   slingshot_pc{n}_tests.csv     condiments tests
# Run: sbatch analysis/scripts/7.3a_slingshot.sb
suppressPackageStartupMessages({ library(slingshot) })
SEED <- 1234; set.seed(SEED)
od <- Sys.getenv("OUT_DIR"); rd <- Sys.getenv("REPO_DIR")
sub <- file.path(od, "day_analysis", "D15", "subcluster_endocrine")
out <- file.path(od, "day_analysis", "D15", "trajectory"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(out, pattern = "^slingshot_", full.names = TRUE))

pca <- as.matrix(read.csv(gzfile(file.path(sub, "SCT_B_pca.csv.gz")), row.names = 1, check.names = FALSE))
sc <- read.csv(file.path(sub, "subclusters.csv"), row.names = 1, check.names = FALSE)
stopifnot(setequal(rownames(pca), rownames(sc)))
sc <- sc[rownames(pca), ]
nm <- read.delim(file.path(rd, "study_config", "annotation", "D15_endocrine_r0.6.tsv"), colClasses = "character")
lab <- nm$label[match(as.character(sc$sub_r0p6), nm$subcluster)]
stopifnot(!anyNA(lab))
geno <- ifelse(grepl("-KO_", rownames(pca)), "KO", "ctr")
u2 <- as.matrix(sc[, c("UMAP1", "UMAP2")]); u3 <- as.matrix(sc[, c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3")])
has_cond <- requireNamespace("condiments", quietly = TRUE)
cat(nrow(pca), "cells;", table(geno), "; condiments:", has_cond, "\n")

curves_out <- function(emb, f) {
  cv <- slingCurves(emb, as.df = FALSE)
  d <- do.call(rbind, lapply(names(cv), function(l) {
    s <- cv[[l]]$s[cv[[l]]$ord, , drop = FALSE]
    data.frame(lineage = l, point = seq_len(nrow(s)), s, check.names = FALSE)
  }))
  write.csv(d, f, row.names = FALSE)
}

for (n in c(30, 15)) {
  set.seed(SEED)
  sds <- slingshot(pca[, seq_len(n)], clusterLabels = lab, start.clus = "ENP1", approx_points = 150)
  ln <- slingLineages(sds)
  cat("\nPCs 1-", n, ": ", length(ln), " lineages\n", sep = "")
  for (l in names(ln)) cat(" ", l, ":", paste(ln[[l]], collapse = " -> "), "\n")
  write.csv(data.frame(lineage = names(ln), clusters = sapply(ln, paste, collapse = " -> ")),
            file.path(out, sprintf("slingshot_pc%d_lineages.csv", n)), row.names = FALSE)
  pt <- slingPseudotime(sds); w <- slingCurveWeights(sds)
  colnames(pt) <- paste0("pt_", colnames(pt)); colnames(w) <- paste0("w_", colnames(w))
  write.csv(data.frame(label = lab, genotype = geno, pt, w, check.names = FALSE),
            file.path(out, sprintf("slingshot_pc%d_cells.csv", n)))
  curves_out(embedCurves(sds, u2), file.path(out, sprintf("slingshot_pc%d_curves_2d.csv", n)))
  curves_out(embedCurves(sds, u3), file.path(out, sprintf("slingshot_pc%d_curves_3d.csv", n)))
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
      print(tt); write.csv(tt, file.path(out, sprintf("slingshot_pc%d_tests.csv", n)), row.names = FALSE)
    }
  }
}
cat("wrote", out, "\n")
