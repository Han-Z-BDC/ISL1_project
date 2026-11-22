# Step 7.3m (combine): join the full D15 tradeSeq chunks (after all 7.3m_tradeseq_fit.R tasks).
#
# FDR (BH) over all genes per test x lineage (x range); lineage names (lineage1 = beta, lineage2 = alpha,
# lineage3 = EC). Smoothed curves (50 points per lineage x genotype) for: genes with a shared-range condition FDR
# < 0.05 on any lineage, the 7.3l pseudobulk genes with curves (pb_genes.csv in_curves), and the endocrine markers.
# Effect size from those curves: max |log2((KO + 1) / (control + 1))| over the well-populated KO range of the lineage:
# from the later of the two genotypes' 5th percentile of pseudotime on the lineage (curve weight >= 0.5) to the KO
# cells' 95th percentile (ko_range.csv). (user, 2026-11-21: the first version measured from the curve start, where
# few cells sit and the smoothers are unreliable; 255 of 362 beta calls peaked there.) Signed log2FC and pseudotime of
# that maximum are kept; tradeSeq call = shared FDR < 0.05 and effect >= 0.5 (as 7.3l).
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/tradeseq_full/ (reruns overwrite):
#   tradeseq_association.csv, tradeseq_condition.csv, tradeseq_earlyDE.csv   (gene, test, which, lineage_name, range,
#                                                                             stat, df, p, FDR)
#   tradeseq_condition_wide.csv   one row per gene: FDR shared / whole, effect, signed log2FC, direction and
#                                 pseudotime of the maximum per lineage, tradeSeq call per lineage, TF (Lambert et al.
#                                 2018 list), 7.3l pseudobulk calls (curve / past stop / KO state) next to it
#   tradeseq_effect_range.csv     per lineage: pseudotime range used for the effect
#   tradeseq_smooth.csv           gene, lineage, lineage_name, condition, time, yhat
# Usage: Rscript 7.3m_tradeseq_combine.R
suppressPackageStartupMessages({ library(tradeSeq); library(SingleCellExperiment) })
LFC <- 0.5
tr <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "D15", "trajectory")
fo <- file.path(tr, "tradeseq_full")
N <- as.integer(readLines(file.path(fo, "chunks.txt")))
unlink(list.files(fo, pattern = "^tradeseq_.*\\.csv$", full.names = TRUE))
NAMES <- c("beta", "alpha", "EC")
MARKERS <- c("ETV1", "LMO2", "ISL1", "PDX1", "NKX6-1", "ARX", "MAFB", "FEV", "LMX1A", "MNX1", "INS", "GCG", "SST",
             "TPH1", "NEUROG3", "NEUROD1", "PAX4", "PAX6", "IRX2", "CHGA", "HHEX", "PHOX2A", "TOX")
kr <- read.csv(file.path(fo, "ko_range.csv"))
inp <- readRDS(file.path(fo, "inputs.rds"))
kr$from <- sapply(seq_along(NAMES), function(i) max(sapply(c("ctr", "KO"), function(g)
  quantile(inp$pt[inp$cond == g & inp$w[, i] >= 0.5, i], 0.05))))
write.csv(kr, file.path(fo, "tradeseq_effect_range.csv"), row.names = FALSE)
print(kr)
TF <- readLines(file.path(Sys.getenv("REPO_DIR"), "study_config", "annotation", "TF_names_v_1.01.txt"))

nm <- function(w) {   # lineage1 -> beta; 1vs2 -> beta_vs_alpha; lineage1_conditionKO -> beta_KO
  w <- as.character(w)
  out <- w
  p <- regmatches(w, regexec("^(\\d)vs(\\d)$", w))
  for (i in seq_along(w)) {
    if (length(p[[i]]) == 3) out[i] <- paste0(NAMES[as.integer(p[[i]][2])], "_vs_", NAMES[as.integer(p[[i]][3])])
    else if (grepl("^lineage\\d$", w[i])) out[i] <- NAMES[as.integer(sub("lineage", "", w[i]))]
    else if (grepl("^lineage\\d_condition", w[i])) out[i] <- paste0(NAMES[as.integer(substr(w[i], 8, 8))], "_",
                                                                      sub("^lineage\\d_condition", "", w[i]))
    else if (grepl("^\\d_condition", w[i])) out[i] <- sub("^(\\d)_condition", "\\1_", w[i]) |>
      (\(z) paste0(NAMES[as.integer(substr(z, 1, 1))], substring(z, 2)))()
  }
  out
}
join <- function(what) {
  fs <- file.path(fo, sprintf("chunk%d_%s.csv", seq_len(N) - 1, what))
  miss <- fs[!file.exists(fs)]
  if (what != "earlyDE" && length(miss)) stop("missing chunk outputs: ", paste(basename(miss), collapse = " "))
  fs <- fs[file.exists(fs)]
  d <- do.call(rbind, lapply(fs, read.csv, check.names = FALSE))
  if (!"range" %in% names(d)) d$range <- ""
  d$lineage_name <- nm(d$which)
  d$FDR <- ave(d$pvalue, d$which, d$range, FUN = function(p) p.adjust(p, "fdr"))
  d <- d[order(d$which, d$range, d$pvalue), ]
  write.csv(d, file.path(fo, sprintf("tradeseq_%s.csv", what)), row.names = FALSE)
  print(aggregate(FDR ~ test + lineage_name + range, d, function(q) sum(q < 0.05, na.rm = TRUE)))
  d
}
as_res <- join("association"); co_res <- join("condition"); ed_res <- join("earlyDE")
genes <- unique(as_res$gene)
cat("genes tested:", length(genes), "\n")

wide <- data.frame(gene = genes)
for (l in NAMES) for (r in c("shared", "whole")) {
  x <- co_res[co_res$lineage_name == l & co_res$range == r, ]
  wide[[paste0("FDR_", r, "_", l)]] <- x$FDR[match(genes, x$gene)]
}
pbg <- read.csv(file.path(tr, "pseudobulk", "pb_genes.csv"), na.strings = character(0))
sig_shared <- unique(co_res$gene[co_res$range == "shared" & !is.na(co_res$FDR) & co_res$FDR < 0.05])
sm_genes <- unique(c(sig_shared, pbg$gene[pbg$in_curves %in% c(TRUE, "TRUE")], intersect(MARKERS, genes)))
cat("smoothed curves for", length(sm_genes), "genes (tradeSeq shared-range calls", length(sig_shared),
    "+ pseudobulk curve genes + markers)\n")
sm <- do.call(rbind, lapply(seq_len(N) - 1, function(k) {
  fit <- readRDS(file.path(fo, sprintf("chunk%d_models.rds", k)))
  g <- intersect(sm_genes, rownames(fit))
  if (!length(g)) return(NULL)
  tryCatch(predictSmooth(fit, gene = g, nPoints = 50, tidy = TRUE),
           error = function(e) { cat("predictSmooth chunk", k, "failed:", conditionMessage(e), "\n"); NULL })
}))
sm$lineage_name <- NAMES[sm$lineage]
write.csv(sm, file.path(fo, "tradeseq_smooth.csv"), row.names = FALSE)

for (i in seq_along(NAMES)) {
  # KO and control curves come on their own time grids: interpolate the control curve at the KO time points inside
  # both curves' ranges and the KO range of the lineage
  s <- sm[sm$lineage == i, ]
  # evaluated on a fixed grid of 30 points inside the well-populated range, both curves interpolated
  grid <- seq(kr$from[i], kr$ko_q95[i], length.out = 30)
  ef <- t(vapply(split(s, s$gene), function(d) {
    k_ <- d[d$condition == "KO", ]; c_ <- d[d$condition == "ctr", ]
    g <- grid[grid >= max(min(k_$time), min(c_$time)) & grid <= min(max(k_$time), max(c_$time))]
    if (length(g) < 2) return(c(NA_real_, NA_real_))
    r <- log2((approx(k_$time, k_$yhat, xout = g)$y + 1) / (approx(c_$time, c_$yhat, xout = g)$y + 1))
    j <- which.max(abs(r)); c(r[j], g[j])
  }, numeric(2)))
  ef0 <- matrix(NA_real_, nrow(wide), 2, dimnames = list(wide$gene, NULL))   # genes without curves stay NA
  ef0[intersect(rownames(ef), wide$gene), ] <- ef[intersect(rownames(ef), wide$gene), ]
  ef <- ef0
  wide[[paste0("effect_", NAMES[i])]] <- round(abs(ef[wide$gene, 1]), 3)
  wide[[paste0("log2FC_", NAMES[i])]] <- round(ef[wide$gene, 1], 3)
  wide[[paste0("direction_", NAMES[i])]] <- ifelse(is.na(ef[wide$gene, 1]), NA,
                                                   ifelse(ef[wide$gene, 1] > 0, "higher in KO", "lower in KO"))
  wide[[paste0("pt_max_", NAMES[i])]] <- round(ef[wide$gene, 2], 1)
  fs <- wide[[paste0("FDR_shared_", NAMES[i])]]
  wide[[paste0("call_", NAMES[i])]] <- !is.na(fs) & fs < 0.05 & !is.na(wide[[paste0("effect_", NAMES[i])]]) &
    wide[[paste0("effect_", NAMES[i])]] >= LFC
}
wide$TF <- wide$gene %in% TF
wide$tradeseq_lineages <- apply(wide[, paste0("call_", NAMES)], 1, function(r) paste(NAMES[r], collapse = ","))
wide <- merge(wide, pbg[, c("gene", "curve", "after_stop", "ko_state", "in_curves", "old_k10_cluster")], by = "gene",
              all.x = TRUE)
wide$pb_curve_genes <- wide$in_curves %in% c(TRUE, "TRUE")
write.csv(wide[order(wide$tradeseq_lineages == "", wide$gene), ], file.path(fo, "tradeseq_condition_wide.csv"),
          row.names = FALSE)

cat("\ntradeSeq calls (shared-range FDR < 0.05 and effect >= 0.5) by lineage combination:\n")
print(sort(table(ifelse(wide$tradeseq_lineages == "", "none", wide$tradeseq_lineages)), decreasing = TRUE))
ts_any <- wide$tradeseq_lineages != ""
cat("\nagreement with the 7.3l pseudobulk curve test (any lineage):\n")
pb_curve <- wide$curve != "" & !is.na(wide$curve)
print(table(tradeSeq = ts_any, pseudobulk_curve = pb_curve))
for (g in intersect(c("ETV1", "LMO2", "ISL1", "PDX1", "NKX6-1", "ARX", "FEV", "MNX1"), wide$gene)) {
  r <- wide[wide$gene == g, ]
  cat(sprintf("%-7s tradeSeq: %-14s | FDR shared b/a/EC %s | effect %s | pseudobulk curve: %s | past stop: %s\n", g,
              ifelse(r$tradeseq_lineages == "", "-", r$tradeseq_lineages),
              paste(signif(unlist(r[paste0("FDR_shared_", NAMES)]), 2), collapse = "/"),
              paste(unlist(r[paste0("effect_", NAMES)]), collapse = "/"), r$curve, r$after_stop))
}
cat("wrote", fo, "\n")
