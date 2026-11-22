# Step 7.3k (combine): join the D15 subset tradeSeq chunks (after all 7.3k_tradeseq_subset_fit.R tasks).
#
# FDR (BH) over all genes of the subset per lineage / pair (as 7.3e, but over the ~3,400 subset genes, not all
# genes); lineage names added (lineage1 = beta, lineage2 = alpha, lineage3 = EC); 7.3j driver calls and the 7.3h
# cluster joined to each row. For the genes with FDR < 0.05 in any condition test: smoothed curves (50 points per
# lineage x genotype) and modules by the shape of the control curve per lineage (scaled log1p, k-means k = 6,
# seed 1234), as 7.3e.
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/ (reruns overwrite):
#   tradeseq_subset_association.csv, tradeseq_subset_condition.csv, tradeseq_subset_earlyDE.csv
#       gene, test, which, lineage_name, stat, df, p, FDR, driver_alpha / driver_beta / driver_EC, ko_vs_ctr, marker
#   tradeseq_subset_condition_wide.csv   one row per gene: condition-test FDR per lineage (beta, alpha, EC) and the
#                                        genotype difference call (FDR < 0.05) per lineage
#   tradeseq_subset_smooth.csv, tradeseq_subset_modules.csv
# Usage: Rscript 7.3k_tradeseq_subset_combine.R <N>
suppressPackageStartupMessages({ library(tradeSeq); library(SingleCellExperiment) })
N <- as.integer(commandArgs(TRUE)[1]); stopifnot(!is.na(N))
SEED <- 1234
tr <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "D15", "trajectory")
ch <- file.path(tr, "tradeseq_subset_chunks")
old <- list.files(tr, pattern = "^tradeseq_subset_.*\\.csv$", full.names = TRUE)
unlink(old[basename(old) != "tradeseq_subset_genes.csv"])   # keep the 7.3j input
NAMES <- c(lineage1 = "beta", lineage2 = "alpha", lineage3 = "EC")
gs <- read.csv(file.path(tr, "tradeseq_subset_genes.csv"))

nm <- function(w) {   # lineage1 -> beta; 1vs2 or lineage1vslineage2 -> beta_vs_alpha; all -> all
  w <- sub("^lineage(\\d)$", "\\1", gsub("lineage", "", as.character(w)))
  pr <- regmatches(w, regexec("^(\\d)vs(\\d)$", w))
  vapply(seq_along(w), function(i) {
    if (length(pr[[i]]) == 3) paste0(NAMES[as.integer(pr[[i]][2])], "_vs_", NAMES[as.integer(pr[[i]][3])])
    else if (grepl("^\\d$", w[i])) NAMES[as.integer(w[i])] else w[i]
  }, "")
}
join <- function(what) {
  fs <- file.path(ch, sprintf("chunk%d_%s.csv", seq_len(N) - 1, what))
  if (what != "earlyDE") stopifnot(all(file.exists(fs)))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) return(NULL)
  d <- do.call(rbind, lapply(fs, read.csv, check.names = FALSE))
  d$FDR <- ave(d$pvalue, d$which, FUN = function(p) p.adjust(p, "fdr"))
  d$lineage_name <- nm(d$which)
  d <- merge(d, gs, by = "gene", all.x = TRUE)
  d <- d[order(d$which, d$pvalue), ]
  write.csv(d, file.path(tr, sprintf("tradeseq_subset_%s.csv", what)), row.names = FALSE)
  print(aggregate(FDR ~ test + lineage_name, d, function(q) sum(q < 0.05, na.rm = TRUE)))
  d
}
as_res <- join("association"); co_res <- join("condition"); ed_res <- join("earlyDE")
cat("genes tested:", length(unique(as_res$gene)), "\n")

per <- co_res[co_res$lineage_name %in% NAMES, ]
wide <- reshape(per[, c("gene", "lineage_name", "FDR")], idvar = "gene", timevar = "lineage_name", direction = "wide")
names(wide) <- sub("^FDR\\.", "FDR_", names(wide))
for (l in NAMES) wide[[paste0("diff_", l)]] <- !is.na(wide[[paste0("FDR_", l)]]) & wide[[paste0("FDR_", l)]] < 0.05
wide$diff_lineages <- apply(wide[, paste0("diff_", NAMES)], 1, function(r) paste(NAMES[r], collapse = ","))
wide <- merge(wide, gs, by = "gene", all.x = TRUE)
write.csv(wide, file.path(tr, "tradeseq_subset_condition_wide.csv"), row.names = FALSE)
cat("genotype difference (condition FDR < 0.05) by lineage combination:\n")
print(sort(table(ifelse(wide$diff_lineages == "", "none", wide$diff_lineages)), decreasing = TRUE))

sig <- unique(per$gene[!is.na(per$FDR) & per$FDR < 0.05])
cat("genes with FDR < 0.05 in a per-lineage condition test:", length(sig), "\n")
if (length(sig)) {
  sm <- do.call(rbind, lapply(seq_len(N) - 1, function(k) {
    fit <- readRDS(file.path(ch, sprintf("chunk%d_models.rds", k)))
    g <- intersect(sig, rownames(fit))
    if (!length(g)) return(NULL)
    tryCatch(predictSmooth(fit, gene = g, nPoints = 50, tidy = TRUE),
             error = function(e) { cat("predictSmooth chunk", k, "failed:", conditionMessage(e), "\n"); NULL })
  }))
  if (!is.null(sm) && nrow(sm)) {
    sm$lineage_name <- NAMES[paste0("lineage", sm$lineage)]
    write.csv(sm, file.path(tr, "tradeseq_subset_smooth.csv"), row.names = FALSE)
    mods <- list()
    for (k in names(NAMES)) {
      g <- per$gene[per$which == k & !is.na(per$FDR) & per$FDR < 0.05]
      li <- as.integer(sub("^lineage", "", k))
      if (length(g) < 12) next
      s <- sm[sm$gene %in% g & sm$lineage == li & sm$condition == "ctr", ]
      y <- tapply(s$yhat, list(s$gene, s$time), mean)
      y <- y[, order(as.numeric(colnames(y)))]
      z <- t(scale(t(log1p(y)))); z[is.na(z)] <- 0
      set.seed(SEED)
      km <- kmeans(z, centers = min(6, nrow(z) - 1), nstart = 20)
      mods[[k]] <- data.frame(gene = rownames(z), lineage = NAMES[[k]], module = km$cluster)
    }
    if (length(mods)) write.csv(do.call(rbind, mods), file.path(tr, "tradeseq_subset_modules.csv"), row.names = FALSE)
  }
}
cat("wrote", tr, "\n")
