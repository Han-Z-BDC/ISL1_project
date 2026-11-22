# Step 7.3l (tests): KO vs control on replicate pseudobulk, all genes, edgeR quasi-likelihood (after
# 7.3l_pseudobulk_prep.py).
#
# Decisions (user, 2026-11-21): one simple, comprehensive, replicate-aware test in three comparisons (see the prep
# header for the samples). Samples with < MIN_REP_CELLS cells are dropped. Per test: filterByExpr on its samples, TMM,
# ~0 + group, estimateDisp, glmQLFit (robust). A call needs FDR < 0.05 and an effect of |log2FC| >= LFC (0.5, as
# 7.3f / 7.3g) somewhere.
#   1 curve       per lineage, windows where all four replicates have >= MIN_REP_CELLS cells ("shared" windows;
#                 windows without both genotypes are listed as untested, e.g. the late beta / alpha path KO never
#                 reaches). d_w = KO - control in window w.
#                   any   joint F-test of all d_w = 0: KO changes the curve in any way (level, shape, timing)
#                   level mean of d_w = 0 (same shape, shifted up / down)
#                   shape d_w - d_1 = 0 for all w (the difference changes along the lineage: timing / shape)
#                   per window: d_w with its own FDR
#                 call = (any FDR < 0.05 and |d_w| >= LFC in at least one window) OR (user, 2026-11-21) at least
#                 one window with its own FDR < 0.05 and |d_w| >= LFC (a narrow difference the joint test dilutes);
#                 pattern = level / shape / level + shape / any only / window only; direction from the windows with FDR < 0.05 and |d_w| >= LFC (higher in KO
#                 / lower in KO / mixed)
#   2 past stop   beta, alpha: rise = control past the stop - control at the stop; gap = KO at the stop - control past
#                 the stop; stop = KO at the stop - control at the stop (side column)
#                   missing in KO      rise >= LFC and gap <= -LFC (both FDR < 0.05): control switches the gene on
#                                      right after the point where KO stops, KO cells at the stop do not have it
#                   kept on in KO      rise <= -LFC and gap >= LFC (both FDR < 0.05): control switches it off right
#                                      after the stop, KO cells at the stop still have it
#   3 KO states   KO state - nearest control state (KO_only vs ENP_FEV_high, pre_alpha_KO vs pre_alpha,
#                 ENP_PHOX2A+ vs ENP_PHOX2A+): higher / lower in KO state
# Writes ${OUT_DIR}/day_analysis/D15/trajectory/pseudobulk/ (reruns overwrite):
#   pb_curve.csv        gene x lineage: F / p / FDR for any, level, shape; log2FC and FDR per window; call, pattern,
#                       direction, windows (significant), max_abs_log2FC
#   pb_curve_windows.csv  lineage x window: tested, pseudotime range, cells per replicate
#   pb_after_stop.csv   gene x lineage: rise, gap, stop (log2FC, FDR), call
#   pb_ko_states.csv    gene x state pair: log2FC, logCPM, F, p, FDR, call
#   pb_genes.csv        one row per gene called anywhere: calls per comparison, old 7.3h k = 10 cluster and source
#                       categories (comparison with the 2,656-gene set)
# Usage: Rscript 7.3l_pseudobulk_edger.R
suppressPackageStartupMessages({ library(edgeR) })
MIN_REP_CELLS <- 10; LFC <- 0.5; FDR_MAX <- 0.05
pb <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "D15", "trajectory", "pseudobulk")
tr <- dirname(pb)
st <- read.csv(file.path(pb, "pb_samples.csv"), check.names = FALSE)
cnt <- as.matrix(read.csv(file.path(pb, "pb_counts.csv.gz"), row.names = 1, check.names = FALSE))
stopifnot(identical(colnames(cnt), st$sample))
cat(nrow(cnt), "genes,", ncol(cnt), "samples\n")

fit_group <- function(s) {   # samples (rows of st) -> QL fit on ~0 + group
  y <- DGEList(cnt[, s$sample, drop = FALSE], group = factor(s$group))
  y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  des <- model.matrix(~ 0 + group, data = data.frame(group = factor(s$group)))
  colnames(des) <- levels(factor(s$group))
  y <- estimateDisp(y, des)
  list(y = y, des = des, fit = glmQLFit(y, des, robust = TRUE))
}
qtest <- function(f, con) {   # contrast (vector or matrix over the group columns) -> table
  con <- as.matrix(con); rownames(con) <- colnames(f$des)
  topTags(glmQLFTest(f$fit, contrast = con), n = Inf, sort.by = "none")$table
}
cvec <- function(f, plus, minus) { v <- setNames(numeric(ncol(f$des)), colnames(f$des)); v[plus] <- 1; v[minus] <- -1; v }

# ---- 1 curve
curve <- list(); wtab <- list()
for (lin in c("beta", "alpha", "EC")) {
  s <- st[st$comparison == "curve" & st$set == lin, ]
  ok_w <- sort(unique(s$window[s$cells >= MIN_REP_CELLS]))
  ok_w <- ok_w[vapply(ok_w, function(w) sum(s$window == w & s$cells >= MIN_REP_CELLS) == 4, TRUE)]
  for (w in sort(unique(s$window))) {
    sw <- s[s$window == w, ]
    wtab[[length(wtab) + 1]] <- data.frame(lineage = lin, window = w, pt_from = sw$pt_from[1], pt_to = sw$pt_to[1],
      tested = w %in% ok_w, t(tapply(sw$cells, sw$replicate, sum)), check.names = FALSE)
  }
  cat(sprintf("\n== curve %s: windows tested %s of %d\n", lin, paste(ok_w, collapse = ","), length(unique(s$window))))
  s <- s[s$window %in% ok_w & s$cells >= MIN_REP_CELLS, ]   # the prep writes every genotype x replicate (empty ones too)
  stopifnot(all(table(s$group) == 2))
  f <- fit_group(s)
  D <- sapply(ok_w, function(w) cvec(f, paste0("KO_w", w), paste0("ctr_w", w)))
  colnames(D) <- paste0("w", ok_w)
  any_ <- qtest(f, D)
  lev <- qtest(f, rowMeans(D))
  shp <- if (length(ok_w) > 1) qtest(f, D[, -1, drop = FALSE] - D[, 1]) else NULL
  out <- data.frame(gene = rownames(any_), lineage = lin, logCPM = round(any_$logCPM, 3),
                    F_any = any_$F, p_any = any_$PValue, FDR_any = any_$FDR,
                    log2FC_level = round(lev$logFC, 3), p_level = lev$PValue, FDR_level = lev$FDR,
                    p_shape = if (is.null(shp)) NA else shp$PValue, FDR_shape = if (is.null(shp)) NA else shp$FDR)
  L <- matrix(NA, nrow(out), length(ok_w)); Q <- L
  for (j in seq_along(ok_w)) {
    t1 <- qtest(f, D[, j])
    L[, j] <- t1$logFC; Q[, j] <- t1$FDR
    out[[paste0("log2FC_w", ok_w[j])]] <- round(t1$logFC, 3)
    out[[paste0("FDR_w", ok_w[j])]] <- signif(t1$FDR, 3)
  }
  hit <- Q < FDR_MAX & abs(L) >= LFC
  out$max_abs_log2FC <- round(apply(abs(L), 1, max), 3)
  out$call_joint <- out$FDR_any < FDR_MAX & out$max_abs_log2FC >= LFC
  out$call_window <- rowSums(hit) > 0
  out$call <- out$call_joint | out$call_window
  sl <- out$FDR_level < FDR_MAX; ss <- !is.na(out$FDR_shape) & out$FDR_shape < FDR_MAX
  out$pattern <- ifelse(!out$call, "-", ifelse(!out$call_joint, "window only", ifelse(sl & ss, "level + shape",
                 ifelse(sl, "level", ifelse(ss, "shape", "any only")))))
  up <- rowSums(hit & L > 0); dn <- rowSums(hit & L < 0)
  out$direction <- ifelse(!out$call, "-", ifelse(up > 0 & dn > 0, "mixed", ifelse(up > 0, "higher in KO",
                   ifelse(dn > 0, "lower in KO", "no single window"))))
  out$windows <- apply(hit, 1, function(h) paste0("w", ok_w[h], collapse = ","))
  out$windows[!out$call] <- ""
  curve[[lin]] <- out
  cat(sprintf("  genes tested %d; called %d (joint test %d; single window only %d)\n", nrow(out), sum(out$call),
              sum(out$call_joint), sum(out$call & !out$call_joint)))
  print(table(pattern = out$pattern[out$call], direction = out$direction[out$call]))
}
cv <- do.call(rbind, lapply(curve, function(d) d[, c("gene", "lineage", "logCPM", "F_any", "p_any", "FDR_any",
  "log2FC_level", "p_level", "FDR_level", "p_shape", "FDR_shape", "max_abs_log2FC", "call", "call_joint",
  "call_window", "pattern", "direction",
  "windows")]))
# per-window columns (tested windows differ between lineages; untested = empty)
allw <- unlist(lapply(1:8, function(w) paste0(c("log2FC_w", "FDR_w"), w)))
cvw <- do.call(rbind, lapply(curve, function(d) {
  d[setdiff(allw, names(d))] <- NA
  d[, c("gene", "lineage", allw)]
}))
cv <- merge(cv, cvw, by = c("gene", "lineage"), all.x = TRUE)
cv <- cv[order(cv$lineage, !cv$call, cv$p_any), ]
write.csv(cv, file.path(pb, "pb_curve.csv"), row.names = FALSE)
wt <- do.call(rbind, lapply(wtab, function(d) { d[setdiff(c("ctr_rep1", "ctr_rep2", "KO_rep1", "KO_rep2"), names(d))] <- NA; d }))
write.csv(wt, file.path(pb, "pb_curve_windows.csv"), row.names = FALSE)
print(wt)

# ---- 2 past the stop
after <- list()
for (lin in c("beta", "alpha")) {
  s <- st[st$comparison == "after_stop" & st$set == lin & st$cells >= MIN_REP_CELLS, ]
  print(table(s$group, s$replicate))
  stopifnot(all(table(s$group) >= 2))
  f <- fit_group(s)
  rise <- qtest(f, cvec(f, "ctr_after", "ctr_stop"))
  gap <- qtest(f, cvec(f, "KO_stop", "ctr_after"))
  stp <- qtest(f, cvec(f, "KO_stop", "ctr_stop"))
  d <- data.frame(gene = rownames(rise), lineage = lin, logCPM = round(rise$logCPM, 3),
                  log2FC_rise_ctr = round(rise$logFC, 3), FDR_rise_ctr = rise$FDR,
                  log2FC_KO_stop_vs_ctr_after = round(gap$logFC, 3), FDR_KO_stop_vs_ctr_after = gap$FDR,
                  log2FC_KO_vs_ctr_at_stop = round(stp$logFC, 3), FDR_KO_vs_ctr_at_stop = stp$FDR)
  miss <- d$log2FC_rise_ctr >= LFC & d$FDR_rise_ctr < FDR_MAX & d$log2FC_KO_stop_vs_ctr_after <= -LFC &
    d$FDR_KO_stop_vs_ctr_after < FDR_MAX
  kept <- d$log2FC_rise_ctr <= -LFC & d$FDR_rise_ctr < FDR_MAX & d$log2FC_KO_stop_vs_ctr_after >= LFC &
    d$FDR_KO_stop_vs_ctr_after < FDR_MAX
  d$call <- ifelse(miss, "missing in KO", ifelse(kept, "kept on in KO", "-"))
  after[[lin]] <- d
  cat(sprintf("\n== past the stop %s: %d genes; missing in KO %d, kept on in KO %d\n", lin, nrow(d), sum(miss), sum(kept)))
  for (cl in c("missing in KO", "kept on in KO")) {
    x <- d[d$call == cl, ]; x <- x[order(-abs(x$log2FC_rise_ctr)), ]
    cat("  ", cl, ":", paste(head(x$gene, 40), collapse = ", "), "\n")
  }
}
af <- do.call(rbind, after); af <- af[order(af$lineage, af$call == "-", -abs(af$log2FC_rise_ctr)), ]
write.csv(af, file.path(pb, "pb_after_stop.csv"), row.names = FALSE)

# ---- 3 KO-diverted states
kos <- list()
for (set in unique(st$set[st$comparison == "ko_state"])) {
  s <- st[st$comparison == "ko_state" & st$set == set & st$cells >= MIN_REP_CELLS, ]
  print(table(s$group, s$replicate))
  if (any(table(s$group) < 2) || length(unique(s$group)) < 2) { cat("  ", set, ": fewer than 2 replicates per side, skipped\n"); next }
  f <- fit_group(s)
  kg <- grep("^KO_", colnames(f$des), value = TRUE); cg <- grep("^ctr_", colnames(f$des), value = TRUE)
  t1 <- qtest(f, cvec(f, kg, cg))
  d <- data.frame(gene = rownames(t1), state = set, comparison = paste(kg, "vs", cg), log2FC = round(t1$logFC, 3),
                  logCPM = round(t1$logCPM, 3), F = t1$F, p = t1$PValue, FDR = t1$FDR)
  d$call <- ifelse(d$FDR < FDR_MAX & d$log2FC >= LFC, "higher in KO state",
                   ifelse(d$FDR < FDR_MAX & d$log2FC <= -LFC, "lower in KO state", "-"))
  kos[[set]] <- d
  cat(sprintf("\n== KO state %s (%s): %d genes; higher %d, lower %d\n", set, d$comparison[1], nrow(d),
              sum(d$call == "higher in KO state"), sum(d$call == "lower in KO state")))
  for (cl in c("higher in KO state", "lower in KO state")) {
    x <- d[d$call == cl, ]; x <- x[order(-abs(x$log2FC)), ]
    cat("  ", cl, ":", paste(head(x$gene, 40), collapse = ", "), "\n")
  }
}
ks <- do.call(rbind, kos); ks <- ks[order(ks$state, ks$call == "-", ks$p), ]
write.csv(ks, file.path(pb, "pb_ko_states.csv"), row.names = FALSE)

# ---- gene summary and comparison with the 7.3h set
cc <- cv[cv$call, ]
g1 <- tapply(paste0(cc$lineage, ": ", cc$direction, " (", cc$pattern, ")"), cc$gene, paste, collapse = "; ")
aa <- af[af$call != "-", ]; g2 <- tapply(paste0(aa$lineage, ": ", aa$call), aa$gene, paste, collapse = "; ")
kk <- ks[ks$call != "-", ]; g3 <- tapply(paste0(kk$state, ": ", kk$call), kk$gene, paste, collapse = "; ")
all_g <- sort(unique(c(names(g1), names(g2), names(g3))))
gs <- data.frame(gene = all_g, curve = unname(g1[all_g]), after_stop = unname(g2[all_g]), ko_state = unname(g3[all_g]))
gs[is.na(gs)] <- ""
gs$in_curves <- gs$curve != "" | gs$after_stop != ""   # genes for the re-clustering (curves on the three lineages)
old <- read.csv(file.path(tr, "gene_patterns_k10.csv"))[, c("gene", "cluster", "sources")]
names(old) <- c("gene", "old_k10_cluster", "old_sources")
gs <- merge(gs, old, by = "gene", all.x = TRUE)
dr <- read.csv(file.path(tr, "tradeseq_subset_genes.csv"))[, c("gene", "driver_alpha", "driver_beta", "driver_EC")]
gs <- merge(gs, dr, by = "gene", all.x = TRUE)
for (c in c("driver_alpha", "driver_beta", "driver_EC")) gs[[c]][is.na(gs[[c]])] <- "-"
write.csv(gs, file.path(pb, "pb_genes.csv"), row.names = FALSE)

cat("\n===== summary\n")
cat("curve calls (any lineage):", length(g1), "| past-the-stop calls:", length(g2), "| KO-state calls:", length(g3), "\n")
cat("genes with curves for re-clustering (curve or past stop):", sum(gs$in_curves), "\n")
cat("KO-state-only genes (separate category):", sum(!gs$in_curves), "\n")
src <- strsplit(old$old_sources, "; ")
cat_of <- function(s) { g <- any(startsWith(s, "7.3g")); st_ <- any(startsWith(s, "7.3f") & !grepl("on after stop", s))
  af_ <- any(grepl("on after stop", s)); trimws(paste(if (g) "windows", if (st_) "stop", if (af_) "after-stop")) }
old$cat <- vapply(src, cat_of, "")
old$recovered_curves <- old$gene %in% gs$gene[gs$in_curves]
old$recovered_any <- old$gene %in% gs$gene
cat("old 2,656 genes recovered (in the re-clustering set / anywhere), by old source:\n")
print(do.call(rbind, lapply(split(old, old$cat), function(d) data.frame(genes = nrow(d),
  in_curves = sum(d$recovered_curves), anywhere = sum(d$recovered_any)))))
cat("new genes (not in the old 2,656), in the re-clustering set:", sum(gs$in_curves & is.na(gs$old_k10_cluster)), "\n")
cat("wrote", pb, "\n")
