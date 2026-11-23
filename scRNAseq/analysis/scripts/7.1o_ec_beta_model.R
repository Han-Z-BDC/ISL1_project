# Step 7.1o (edgeR): D22 KO EC vs control EC vs control beta, two layers (after 7.1m prep; user, 2026-11-23).
#
# Question (user): how different are KO EC from control EC and from control beta, and do some KO EC keep beta features
# (INS is spread over the KO EC population, so no INS cut-off: population model + per-cell scores).
# Mode "model" (layer 1): one edgeR QL model of the D22 replicate pseudobulks (7.1m, SoupX counts summed per
#   sample x label, >= 20 cells): groups ctr_EC1, ctr_EC2, ctr_beta, KO_EC1, KO_EC2; design ~ 0 + group + ctr_rep2 +
#   KO_rep2 (replicate effect within genotype: the control groups share samples, the KO groups share samples);
#   filterByExpr (group), TMM, estimateDisp, glmQLFit robust; call = FDR < 0.05 and |log2FC| >= 0.5 (as 7.1m).
#   Contrasts: Q1 KO_EC1 - ctr_EC1, KO_EC2 - ctr_EC2; Q2 KO_EC1 - ctr_beta, KO_EC2 - ctr_beta; reference axis
#   ctr_beta - ctr_EC1, ctr_beta - ctr_EC2.
#   Beta axis (per EC type k): genes with FDR < 0.05 and |log2FC| >= 1 in ctr_beta - ctr_ECk, INS left out. Shift of
#   KO ECk toward beta = least-squares slope (through 0) of log2FC(KO_ECk - ctr_ECk) on log2FC(ctr_beta - ctr_ECk)
#   over the axis genes (0 = like control EC, 1 = like control beta), Spearman r, and quadrant counts of the KO calls
#   (beta-high up / down in KO, EC-high up / down in KO). Control EC carry some beta ambient RNA (control samples have
#   beta cells), so the shift is, if anything, underestimated; genes "down in KO EC" that are beta-high are flagged.
#   Pseudobulk correlation (group means of TMM logCPM, axis genes, centred) between the 5 groups.
#   Score genes for layer 2: beta-high / EC-high in both EC types (log2FC >= 1 / <= -1 vs EC1 and vs EC2, FDR < 0.05),
#   INS left out, top 100 each by the smaller |log2FC|.
# Mode "subset" (layer 2, after 7.1o_ec_beta_cells.py): pseudobulks of beta-leaning KO EC vs the other KO EC of the
#   same EC type (paired within KO sample, ~ sample + group) and vs control beta (unpaired, ~ group), when every KO
#   sample has >= 20 beta-leaning cells of that EC type.
# Writes ${OUT_DIR}/day_analysis/ec_beta_D22/ (reruns overwrite): de_{contrast}.csv, de_wide.csv, beta_shift.csv,
#   quadrants.csv, pb_corr.csv, score_genes.csv, summary.csv; subset mode: de_subset_{test}.csv, subset_summary.csv
# Run: sbatch analysis/scripts/7.1o_ec_beta.sb
suppressPackageStartupMessages({ library(edgeR); library(data.table) })
mode <- commandArgs(TRUE)[1]; stopifnot(mode %in% c("model", "subset"))
da <- file.path(Sys.getenv("OUT_DIR"), "day_analysis"); od <- file.path(da, "ec_beta_D22")
dir.create(od, showWarnings = FALSE, recursive = TRUE)
MIN_CELLS <- 20; LFC <- 0.5; AXIS_LFC <- 1; N_SCORE <- 100

if (mode == "model") {
  cnt <- fread(file.path(da, "ec_similarity", "pb_counts_D22.csv.gz"))
  genes <- cnt[[1]]; m <- as.matrix(cnt[, -1]); rownames(m) <- genes
  s <- fread(file.path(da, "ec_similarity", "pb_samples_D22.csv"))[match(colnames(m), pb)]
  s[, grp := paste(genotype, sub("^SC_", "", label), sep = "_")]
  keep <- s$cells >= MIN_CELLS & s$grp %in% c("ctr_EC1", "ctr_EC2", "ctr_beta", "KO_EC1", "KO_EC2")
  s <- s[keep]; m <- m[, keep]
  print(s[, .(pb, cells, grp)])
  grp <- factor(s$grp, levels = c("ctr_EC1", "ctr_EC2", "ctr_beta", "KO_EC1", "KO_EC2"))
  rep2 <- grepl("rep2$", s$sample)
  des <- cbind(model.matrix(~ 0 + grp), ctr_rep2 = as.numeric(rep2 & s$genotype == "ctr"),
               KO_rep2 = as.numeric(rep2 & s$genotype == "KO"))
  colnames(des) <- sub("^grp", "", colnames(des))
  stopifnot(qr(des)$rank == ncol(des))
  y <- DGEList(m, group = grp)
  y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y); y <- estimateDisp(y, des)
  fit <- glmQLFit(y, des, robust = TRUE)
  cat("genes tested:", nrow(y), "\n")
  CON <- c(KO_EC1_vs_ctr_EC1 = "KO_EC1 - ctr_EC1", KO_EC2_vs_ctr_EC2 = "KO_EC2 - ctr_EC2",
           KO_EC1_vs_ctr_beta = "KO_EC1 - ctr_beta", KO_EC2_vs_ctr_beta = "KO_EC2 - ctr_beta",
           ctr_beta_vs_ctr_EC1 = "ctr_beta - ctr_EC1", ctr_beta_vs_ctr_EC2 = "ctr_beta - ctr_EC2")
  wide <- data.table(gene = rownames(y)); summ <- list()
  for (cn in names(CON)) {
    r <- topTags(glmQLFTest(fit, contrast = makeContrasts(contrasts = CON[[cn]], levels = des)), n = Inf, sort.by = "none")$table
    r <- data.table(gene = rownames(r), r); r[, call := FDR < 0.05 & abs(logFC) >= LFC]
    fwrite(r[order(PValue)], file.path(od, sprintf("de_%s.csv", cn)))
    wide[, (paste0("lfc_", cn)) := r$logFC][, (paste0("FDR_", cn)) := r$FDR]
    summ[[cn]] <- data.table(contrast = cn, up = sum(r$call & r$logFC > 0), down = sum(r$call & r$logFC < 0))
  }
  wide[, logCPM := aveLogCPM(y)]
  fwrite(wide, file.path(od, "de_wide.csv"))
  summ <- rbindlist(summ); print(summ); fwrite(summ, file.path(od, "summary.csv"))

  # beta axis, shift, quadrants per EC type
  sh <- list(); qd <- list()
  for (k in c("EC1", "EC2")) {
    ax <- wide[get(paste0("FDR_ctr_beta_vs_ctr_", k)) < 0.05 & abs(get(paste0("lfc_ctr_beta_vs_ctr_", k))) >= AXIS_LFC &
               gene != "INS"]
    xa <- ax[[paste0("lfc_ctr_beta_vs_ctr_", k)]]; ya <- ax[[sprintf("lfc_KO_%s_vs_ctr_%s", k, k)]]
    fa <- ax[[sprintf("FDR_KO_%s_vs_ctr_%s", k, k)]]
    called <- fa < 0.05 & abs(ya) >= LFC
    sh[[k]] <- data.table(ec = k, axis_genes = nrow(ax), beta_high = sum(xa > 0), EC_high = sum(xa < 0),
                          slope = sum(xa * ya) / sum(xa^2), spearman = cor(xa, ya, method = "spearman"),
                          slope_called_only = if (any(called)) sum(xa[called] * ya[called]) / sum(xa[called]^2) else NA)
    qd[[k]] <- ax[, .(ec = k, gene, beta_vs_EC = xa, KO_vs_ctr_EC = ya, FDR_KO = fa, called,
                      quadrant = fifelse(!called, "not changed in KO",
                        fifelse(xa > 0 & ya > 0, "beta-high, up in KO EC",
                        fifelse(xa > 0 & ya < 0, "beta-high, down in KO EC (check ambient)",
                        fifelse(xa < 0 & ya < 0, "EC-high, down in KO EC", "EC-high, up in KO EC")))))]
  }
  sh <- rbindlist(sh); qd <- rbindlist(qd)
  print(sh); print(qd[, .N, by = .(ec, quadrant)][order(ec, quadrant)])
  fwrite(sh, file.path(od, "beta_shift.csv")); fwrite(qd, file.path(od, "quadrants.csv"))

  # group-mean correlation on the axis genes (union of both EC types), centred over the 5 groups
  axg <- unique(qd$gene)
  lc <- cpm(y, log = TRUE, prior.count = 2)[axg, ]
  gm <- sapply(levels(grp), function(g) rowMeans(lc[, grp == g, drop = FALSE]))
  z <- gm - rowMeans(gm)
  fwrite(data.table(group = colnames(z), round(cor(z), 4)), file.path(od, "pb_corr.csv"))
  print(round(cor(z), 3))

  # score genes for layer 2
  b <- wide[gene != "INS" & FDR_ctr_beta_vs_ctr_EC1 < 0.05 & FDR_ctr_beta_vs_ctr_EC2 < 0.05]
  b[, m1 := pmin(lfc_ctr_beta_vs_ctr_EC1, lfc_ctr_beta_vs_ctr_EC2)][, m2 := pmax(lfc_ctr_beta_vs_ctr_EC1, lfc_ctr_beta_vs_ctr_EC2)]
  bh <- b[m1 >= AXIS_LFC][order(-m1)][seq_len(min(.N, N_SCORE))]
  eh <- b[m2 <= -AXIS_LFC][order(m2)][seq_len(min(.N, N_SCORE))]
  sg <- rbind(bh[, .(set = "beta", gene, min_abs_lfc = m1)], eh[, .(set = "EC", gene, min_abs_lfc = -m2)])
  fwrite(sg, file.path(od, "score_genes.csv"))
  cat("score genes: beta", nrow(bh), "| EC", nrow(eh), "\n")
} else {
  f <- file.path(od, "pb_subset_counts.csv.gz")
  cnt <- fread(f); genes <- cnt[[1]]; m <- as.matrix(cnt[, -1]); rownames(m) <- genes
  s <- fread(file.path(od, "pb_subset_samples.csv"))[match(colnames(m), pb)]
  out <- list()
  run <- function(cols, g, des, coef, name) {
    y <- DGEList(m[, cols], group = g)
    y <- y[filterByExpr(y, group = g), , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y); y <- estimateDisp(y, des)
    fit <- glmQLFit(y, des, robust = TRUE)
    r <- topTags(glmQLFTest(fit, coef = coef), n = Inf)$table
    r <- data.table(gene = rownames(r), r); r[, call := FDR < 0.05 & abs(logFC) >= LFC]
    fwrite(r, file.path(od, sprintf("de_subset_%s.csv", name)))
    out[[name]] <<- data.table(test = name, pseudobulks = length(cols), genes = nrow(y),
                               up = sum(r$call & r$logFC > 0), down = sum(r$call & r$logFC < 0))
  }
  for (k in c("EC1", "EC2")) {
    lean <- s[group == sprintf("KO_%s_beta_leaning", k)]; rest <- s[group == sprintf("KO_%s_other", k)]
    if (nrow(lean) < 2 || any(lean$cells < MIN_CELLS)) { cat(k, ": too few beta-leaning cells per sample, skipped\n"); next }
    cols <- c(lean$pb, rest$pb)
    g <- factor(ifelse(cols %in% lean$pb, "lean", "other"), levels = c("other", "lean"))
    smp <- factor(s$sample[match(cols, s$pb)])
    run(cols, g, model.matrix(~ smp + g), "glean", sprintf("KO_%s_beta_leaning_vs_other", k))
    beta <- s[group == "ctr_beta"]
    cols <- c(beta$pb, lean$pb); g <- factor(ifelse(cols %in% lean$pb, "lean", "beta"), levels = c("beta", "lean"))
    run(cols, g, model.matrix(~ g), "glean", sprintf("KO_%s_beta_leaning_vs_ctr_beta", k))
  }
  if (length(out)) { o <- rbindlist(out); print(o); fwrite(o, file.path(od, "subset_summary.csv")) }
}
