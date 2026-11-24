# Step 12: ISL1 binding of the D15 scRNA-seq trajectory gene clusters (user, 2026-11-24).
# Genes: the 3,118 KO vs control trajectory genes in k = 10 clusters (scRNAseq 7.3h, gene_patterns_k10.csv).
# Background: the other genes tested along the D15 lineages (scRNAseq 7.3l pb_curve.csv), i.e. expressed in the D15
#   endocrine cells but not called.
# Binding: CUT&RUN step 13 links (cutnrun_ISL1_timecourse/abc_targets/links.tsv): a merged ISL1 site (pooled
#   macs3_ge2 consensus of D15 / D18 / D22) within +/- 1 kb of any TSS of the gene ("TSS") or overlapping an ABC
#   element linked to it (ENP2, ENP3, early_SC_alpha; "enhancer only" when no TSS link).
#   1 per stage: bound at a day = >= 1 linked site called on that day; 2 any stage = >= 1 linked site.
#   Each cluster vs background: Fisher's exact test on bound / not bound, odds ratio, BH over all cluster x stage
#   tests (also for TSS-bound alone).
# 3 TSS signal: one TSS per gene (GTF gene record start, strand-aware), +/- 3 kb, 50 bp bins, the per-day mean ISL1
#   bigWigs (CUT&RUN step 08, background-normalised CPM) via deepTools computeMatrix (run by the .sb between the two
#   modes) for clusters C1-C10 and 1,000 background genes (random, seed 1234). Heatmap: clusters only, genes sorted
#   within a cluster by the mean signal over +/- 500 bp and the three days (highest at the top); sized to match the
#   scRNAseq gene-pattern heatmap (3.5 in tall, each day as wide as one beta / alpha curve). Profiles: clusters and
#   background.
# usage: Rscript 12_trajectory_binding.R prep <gene_patterns_k10.csv> <pb_curve.csv> <links.tsv> <genes.gtf> <out dir>
#        Rscript 12_trajectory_binding.R plot <out dir>
#        Rscript 12_trajectory_binding.R mode <out dir> <links.tsv> <pb_curve.csv> <ABC results dir> <groups, comma-separated>
# mode (user, 2026-11-24: which clusters are bound at the TSS or through enhancers; per stage only): per gene and
#   stage, from the sites called that day: TSS = >= 1 site within +/- 1 kb of a TSS of the gene (step 13 source TSS
#   or TSS+ABC), enhancer = >= 1 site linked only through an ABC element (source ABC); categories TSS only, TSS +
#   enhancer, enhancer only, and all (any of the three). Each category: Fisher's exact test cluster vs background
#   (category vs not, all genes), BH within a category over cluster x stage. Adjusted check: logistic regression per stage and
#   category, category ~ group (background = reference) + log2(1 + ABC elements linked to the gene) + mean logCPM
#   along the D15 lineages (7.3l), BH within a category over the cluster x stage terms. ABC elements = distinct non-self-promoter elements
#   linked to the gene (ABC >= 0.02) in any of the ABC groups (the opportunity for an enhancer site).
# Output in <out dir>: genes.tsv (gene, cluster, bound per stage, TSS / enhancer only), binding_by_cluster.tsv,
#   tss.bed (computeMatrix regions), tss_matrix.gz (computeMatrix), tss_signal.tsv (per gene and day, +/- 500 bp),
#   isl1_binding_by_cluster.pdf, isl1_tss_heatmap.pdf, isl1_tss_profile.pdf
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE); mode <- a[1]
DAYS <- c("D15", "D18", "D22"); SEED <- 1234; N_BG <- 1000
FLANK <- 3000; BIN <- 50; CORE <- 500

if (mode == "prep") {
  stopifnot(length(a) == 6); od <- a[6]; dir.create(od, showWarnings = FALSE, recursive = TRUE)
  cl <- fread(a[2], select = c("gene", "cluster"))
  uni <- union(unique(fread(a[3], select = "gene")$gene), cl$gene)
  g <- data.table(gene = uni, cluster = cl$cluster[match(uni, cl$gene)])
  g[, group := fifelse(is.na(cluster), "background", sprintf("C%d", cluster))]
  cat("genes:", nrow(g), "| clusters:", sum(!is.na(g$cluster)), "| background:", sum(is.na(g$cluster)), "\n")

  # gene TSS (GTF gene records); names as in the scRNA-seq (same GTF); duplicated names: first record
  gt <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' %s", a[5]), header = FALSE, sep = "\t", quote = "")
  gt <- gt[, .(chr = V1, start = V4, end = V5, strand = V7,
               gene = sub('.*gene_name "([^"]+)".*', "\\1", V9), gene_id = sub('.*gene_id "([^"]+)".*', "\\1", V9))]
  gt <- gt[grepl("^chr([0-9]+|X|Y)$", chr)][!duplicated(gene)]
  g[gt, `:=`(gene_id = i.gene_id, chr = i.chr, tss = fifelse(i.strand == "+", i.start, i.end), strand = i.strand),
    on = "gene"]
  cat("genes without a GTF TSS (dropped):", sum(is.na(g$tss)), paste(head(g[is.na(tss), gene], 10), collapse = ", "), "\n")
  g <- g[!is.na(tss)]

  # binding from the step 13 links
  ln <- fread(a[4], select = c("site", "days", "gene_name", "gene_id", "source"))
  ln <- ln[gene_id %in% g$gene_id]
  for (d in c("any", DAYS)) {
    s <- if (d == "any") ln else ln[grepl(d, days)]
    b <- s[, .(bound = TRUE, tss = any(grepl("TSS", source))), by = gene_id]
    g[, (paste0("bound_", d)) := gene_id %in% b$gene_id]
    g[, (paste0("tss_", d)) := gene_id %in% b[(tss), gene_id]]
  }
  g[, n_sites := ln[, uniqueN(site), by = gene_id][match(g$gene_id, gene_id), V1]][is.na(n_sites), n_sites := 0L]
  fwrite(g, file.path(od, "genes.tsv"), sep = "\t")

  # cluster vs background, per stage
  st <- rbindlist(lapply(c("any", DAYS), function(d) rbindlist(lapply(setdiff(sort(unique(g$group)), "background"),
    function(k) rbindlist(lapply(c("bound", "tss"), function(w) {
      v <- g[[paste0(w, "_", d)]]; inn <- g$group == k; bg <- g$group == "background"
      ft <- fisher.test(matrix(c(sum(v & inn), sum(!v & inn), sum(v & bg), sum(!v & bg)), 2))
      data.table(stage = d, group = k, measure = fifelse(w == "bound", "TSS or enhancer", "TSS"), n = sum(inn),
                 n_bound = sum(v & inn), frac = mean(v[inn]), frac_background = mean(v[bg]),
                 odds_ratio = unname(ft$estimate), p = ft$p.value)
    }))))))
  st[, padj := p.adjust(p, "BH"), by = measure]
  fwrite(st, file.path(od, "binding_by_cluster.tsv"), sep = "\t")
  print(dcast(st[measure == "TSS or enhancer"], group ~ stage, value.var = "frac", fun.aggregate = function(x) round(x, 3)))
  print(unique(st[, .(stage, measure, frac_background = round(frac_background, 3))]))
  print(st[padj < 0.05, .(stage, group, measure, frac = round(frac, 3), OR = round(odds_ratio, 2), padj = signif(padj, 2))])

  # computeMatrix regions: clusters + N_BG background genes
  set.seed(SEED)
  bgk <- sample(g[group == "background", gene], N_BG)
  r <- g[group != "background" | gene %in% bgk]
  fwrite(r[, .(chr, start = tss - 1L, end = tss, name = gene, score = 0, strand)], file.path(od, "tss.bed"),
         sep = "\t", col.names = FALSE)
  cat("tss.bed:", nrow(r), "regions\n")
}

if (mode == "plot") {
  od <- a[2]
  g <- fread(file.path(od, "genes.tsv")); st <- fread(file.path(od, "binding_by_cluster.tsv"))
  GROUPS <- c(sprintf("C%d", 1:10), "background")

  # ---- 1 + 2: fraction bound per cluster and stage, TSS vs enhancer only
  fr <- rbindlist(lapply(c("any", DAYS), function(d) g[, .(TSS = mean(get(paste0("tss_", d))),
    enh = mean(get(paste0("bound_", d)) & !get(paste0("tss_", d))), n = .N), by = group][, stage := d]))
  fr <- melt(fr, id.vars = c("group", "stage", "n"), variable.name = "link", value.name = "frac")
  fr[, link := factor(fifelse(link == "TSS", "TSS (+/- 1 kb)", "enhancer only (ABC)"),
                      levels = c("enhancer only (ABC)", "TSS (+/- 1 kb)"))]
  SLAB <- c(any = "any stage", D15 = "D15", D18 = "D18", D22 = "D22")
  fr[, `:=`(group = factor(group, levels = rev(GROUPS)), stage = factor(SLAB[stage], levels = SLAB))]
  tot <- st[measure == "TSS or enhancer", .(group, stage = factor(SLAB[stage], levels = SLAB), frac, odds_ratio, padj)]
  tot[, lab := sprintf("%.0f%%%s", 100 * frac, fifelse(padj < 0.001, "***", fifelse(padj < 0.01, "**",
                                                                                     fifelse(padj < 0.05, "*", ""))))]
  tot[, group := factor(group, levels = rev(GROUPS))]
  bgt <- g[group == "background", .(stage = factor(SLAB, levels = SLAB), frac = c(mean(bound_any), mean(bound_D15),
                                                                                  mean(bound_D18), mean(bound_D22)))]
  bgl <- fr[group == "background", .(frac = sum(frac)), by = stage]
  ns <- g[, .N, by = group]; GL <- setNames(sprintf("%s (%s)", ns$group, format(ns$N, big.mark = ",")), ns$group)
  p1 <- ggplot(fr, aes(frac, group, fill = link)) +
    geom_col(width = 0.75) +
    geom_vline(data = bgl, aes(xintercept = frac), linetype = "dashed", linewidth = 0.3, colour = "grey30") +
    geom_text(data = tot, aes(x = frac + 0.02, y = group, label = lab), inherit.aes = FALSE, hjust = 0, size = 1.9) +
    facet_grid(. ~ stage) +
    scale_fill_manual(values = c("enhancer only (ABC)" = "#9ecae1", "TSS (+/- 1 kb)" = "#08519c"), name = NULL) +
    scale_x_continuous(limits = c(0, 1.18), breaks = c(0, 0.5, 1), labels = c("0", "50%", "100%"), expand = c(0, 0)) +
    scale_y_discrete(labels = GL) +
    labs(x = "genes with an ISL1 site", y = NULL, title = "ISL1 binding of the D15 trajectory gene clusters",
         subtitle = paste0("site at the TSS (+/- 1 kb) or in an ABC-linked enhancer; dashed = background (genes tested,\n",
                           "not called); * BH FDR < 0.05, ** < 0.01, *** < 0.001 vs background (Fisher, TSS or enhancer)")) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom",
          plot.title = element_text(face = "bold"), plot.title.position = "plot", legend.key.size = unit(0.3, "cm"),
          panel.spacing.x = unit(0.25, "cm"))
  ggsave(file.path(od, "isl1_binding_by_cluster.pdf"), p1, width = 6, height = 3.2, useDingbats = FALSE)

  # ---- 3: TSS signal heatmap
  m <- fread(cmd = sprintf("zcat %s | tail -n +2", file.path(od, "tss_matrix.gz")), header = FALSE)
  nb <- 2 * FLANK / BIN; stopifnot(ncol(m) == 6 + 3 * nb)
  genes <- m$V4; v <- as.matrix(m[, 7:ncol(m)]); v[is.na(v)] <- 0
  core <- which(abs(seq(-FLANK + BIN / 2, FLANK - BIN / 2, by = BIN)) <= CORE)
  sig <- sapply(seq_along(DAYS), function(i) rowMeans(v[, (i - 1) * nb + core, drop = FALSE]))
  colnames(sig) <- DAYS
  s <- data.table(gene = genes, sig, mean = rowMeans(sig))
  s[g, group := i.group, on = "gene"]
  fwrite(s, file.path(od, "tss_signal.tsv"), sep = "\t")
  print(s[, .(n = .N, D15 = round(median(D15), 2), D18 = round(median(D18), 2), D22 = round(median(D22), 2),
              mean = round(mean(mean), 2)), by = group][order(match(group, GROUPS))])
  setorderv(s, c("mean", "gene"), c(-1, 1))
  s[, group := factor(group, levels = GROUPS)]
  s <- s[order(group)][, row := rev(seq_len(.N))]   # C1 at the top, highest signal first within a block
  pos <- seq(-FLANK + BIN / 2, FLANK - BIN / 2, by = BIN)
  hm <- rbindlist(lapply(seq_along(DAYS), function(i) {
    x <- v[, (i - 1) * nb + seq_len(nb), drop = FALSE]; colnames(x) <- seq_len(nb)
    melt(data.table(gene = genes, day = DAYS[i], x), id.vars = c("gene", "day"), variable.name = "bin",
         value.name = "signal")
  }))
  hm[, pos := pos[as.integer(as.character(bin))]][, bin := NULL]
  cat("heatmap cells:", nrow(hm), "| NA position:", sum(is.na(hm$pos)), "| NA signal:", sum(is.na(hm$signal)),
      "| min signal:", min(hm$signal), "\n")
  hm[s, `:=`(row = i.row, group = i.group), on = "gene"]
  CAP <- quantile(hm$signal, 0.99); cat("colour cap (99th percentile):", signif(CAP, 3), "\n")
  hm[, sc := pmin(pmax(signal, 0), CAP)]
  # heatmap: clusters only (user, 2026-11-24: no background), sized to sit next to the scRNAseq gene-pattern heatmap
  # (fig_d15_gene_patterns.R: heatmap area 2.5 x 3.5 in, 52 bins over 6 curves with 0.06 cm gaps, 0.04 cm gaps between
  # clusters): each day as wide as one beta / alpha curve there (9 bins), panel area 3.5 in tall, cluster heights in
  # proportion to their genes, same gaps
  CL <- setdiff(GROUPS, "background")
  hmc <- hm[group %in% CL]; hmc[, group := factor(as.character(group), levels = CL)]
  GAP_X <- 0.06; GAP_Y <- 0.04                                   # cm
  DAY_W <- (2.5 - 5 * GAP_X / 2.54) / 52 * 9; HM_H <- 3.5        # inches
  ns <- s[group %in% CL, .N, by = group]
  BL <- setNames(sprintf("%s\n%s", ns$group, format(ns$N, big.mark = ",")), as.character(ns$group))
  hmc[, glab := factor(BL[as.character(group)], levels = BL[CL])]
  p3 <- ggplot(hmc, aes(pos, row, fill = sc)) + geom_raster() +
    facet_grid(glab ~ day, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradientn(colours = c("white", "#fee391", "#fe9929", "#cc4c02", "#662506"), limits = c(0, CAP),
                         name = "ISL1\nsignal") +
    scale_x_continuous(breaks = c(-2000, 0, 2000), labels = c("-2", "0", "+2"), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    labs(x = "kb from TSS", y = NULL, title = "ISL1 CUT&RUN at the TSS",
         subtitle = sprintf(paste0("mean of 3 libraries per day,\nbackground-normalised CPM;\nsorted by mean signal ",
                                   "+/- %d bp;\ncapped at the 99th percentile"), CORE)) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(), panel.spacing.x = unit(GAP_X, "cm"), panel.spacing.y = unit(GAP_Y, "cm"),
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 5.5, lineheight = 0.9),
          strip.text.x = element_text(face = "bold", size = 6), strip.placement = "outside", axis.text.x = element_text(size = 5),
          plot.title = element_text(face = "bold"), plot.title.position = "plot", legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.25, "cm"))
  g3 <- ggplotGrob(p3); pp <- g3$layout[grepl("^panel", g3$layout$name), ]
  hts <- ns[match(CL, as.character(group)), N]; PH <- HM_H - (length(CL) - 1) * GAP_Y / 2.54
  for (j in seq_along(unique(pp$t))) g3$heights[sort(unique(pp$t))[j]] <- unit(PH * hts[j] / sum(hts), "in")
  for (i in unique(pp$l)) g3$widths[i] <- unit(DAY_W, "in")
  cat(sprintf("TSS heatmap: day width %.3f in, panel area %.2f in tall\n", DAY_W, HM_H))
  w <- grid::convertWidth(sum(g3$widths), "in", valueOnly = TRUE); h <- grid::convertHeight(sum(g3$heights), "in", valueOnly = TRUE)
  pdf(file.path(od, "isl1_tss_heatmap.pdf"), width = w, height = h, useDingbats = FALSE, bg = "white"); grid::grid.draw(g3); invisible(dev.off())
  cat(sprintf("isl1_tss_heatmap.pdf: %.2f x %.2f in\n", w, h))

  # mean profile per block and day
  pr <- hm[, .(signal = mean(signal)), by = .(group, day, pos)]
  GCOL <- setNames(c(scales::hue_pal()(10), "grey40"), GROUPS)
  p4 <- ggplot(pr, aes(pos, signal, colour = group, linetype = group == "background")) +
    geom_line(linewidth = 0.4) + facet_grid(. ~ day) +
    scale_colour_manual(values = GCOL, name = NULL) + scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"),
                                                                            guide = "none") +
    scale_x_continuous(breaks = c(-2000, 0, 2000), labels = c("-2 kb", "TSS", "+2 kb")) +
    labs(x = NULL, y = "mean ISL1 signal", title = "ISL1 CUT&RUN at the TSS: mean profile per cluster",
         subtitle = "background (dashed) = 1,000 random genes tested, not called") +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
          plot.title.position = "plot", legend.key.size = unit(0.3, "cm"))
  ggsave(file.path(od, "isl1_tss_profile.pdf"), p4, width = 5, height = 2.2, useDingbats = FALSE)
  cat("wrote figures in", od, "\n")
}

if (mode == "mode") {
  stopifnot(length(a) == 6); od <- a[2]
  groups <- strsplit(a[6], ",")[[1]]
  g <- fread(file.path(od, "genes.tsv"))
  GROUPS <- c(sprintf("C%d", 1:10), "background"); CL <- GROUPS[1:10]
  ln <- fread(a[3], select = c("site", "days", "gene_id", "source"))[gene_id %in% g$gene_id]
  # opportunity: ABC elements per gene (not the gene's own promoter), union of the ABC groups
  ab <- unique(rbindlist(lapply(groups, function(k) {
    d <- fread(file.path(a[5], k, "Predictions", "EnhancerPredictionsFull_threshold0.02_self_promoter.tsv"),
               select = c("chr", "start", "end", "TargetGeneEnsembl_ID", "isSelfPromoter"))
    d[isSelfPromoter == FALSE, .(chr, start, end, gene_id = sub("\\..*$", "", TargetGeneEnsembl_ID))]
  })))
  na <- ab[, .(n_abc = .N), by = gene_id]
  g[na, n_abc := i.n_abc, on = "gene_id"][is.na(n_abc), n_abc := 0L]
  pb <- fread(a[4], select = c("gene", "logCPM"))[, .(logCPM = mean(logCPM)), by = gene]
  g[pb, logCPM := i.logCPM, on = "gene"]
  cat("ABC elements per gene, median by group:\n")
  print(g[, .(n = .N, n_abc = as.numeric(median(n_abc)), logCPM = round(median(logCPM, na.rm = TRUE), 2)),
          by = group][order(match(group, GROUPS))])
  CAT <- c("TSS only", "TSS + enhancer", "enhancer only", "all")
  for (d in DAYS) {
    s <- ln[grepl(d, days)]
    tg <- s[grepl("TSS", source), unique(gene_id)]; eg <- s[source == "ABC", unique(gene_id)]
    tt <- g$gene_id %in% tg; ee <- g$gene_id %in% eg
    g[, (paste(d, CAT, sep = "|")) := list(tt & !ee, tt & ee, ee & !tt, tt | ee)]
    stopifnot(all(g[[paste0("bound_", d)]] == (tt | ee)))   # same as the per-stage binding above
  }
  fwrite(g, file.path(od, "genes_mode.tsv"), sep = "\t")
  g[, group := factor(group, levels = rev(GROUPS))]
  # Fisher per cluster x stage x category
  fs <- rbindlist(lapply(DAYS, function(d) rbindlist(lapply(CAT, function(k) rbindlist(lapply(CL, function(cl) {
    v <- g[[paste(d, k, sep = "|")]]; inn <- g$group == cl; bg <- g$group == "background"
    ft <- fisher.test(matrix(c(sum(v & inn), sum(!v & inn), sum(v & bg), sum(!v & bg)), 2))
    data.table(stage = d, category = k, group = cl, n = sum(inn), n_cat = sum(v & inn), frac = mean(v[inn]),
               frac_background = mean(v[bg]), odds_ratio = unname(ft$estimate), p = ft$p.value)
  }))))))
  fs[, padj := p.adjust(p, "BH"), by = category]   # BH within a category, over cluster x stage
  # adjusted: logistic regression
  g[, `:=`(grp = relevel(factor(as.character(group)), ref = "background"), la = log2(1 + n_abc))]
  gl <- rbindlist(lapply(DAYS, function(d) rbindlist(lapply(CAT, function(k) {
    y <- g[[paste(d, k, sep = "|")]]
    m <- glm(y ~ grp + la + logCPM, family = binomial, data = g)
    co <- summary(m)$coefficients
    r <- co[grepl("^grp", rownames(co)), , drop = FALSE]
    data.table(stage = d, category = k, group = sub("^grp", "", rownames(r)), log_or_adj = r[, 1], se = r[, 2],
               p_adj_model = r[, 4], or_abc_per_doubling = exp(co["la", 1]), or_logCPM = exp(co["logCPM", 1]),
               n_model = nobs(m))
  }))))
  gl[, `:=`(odds_ratio_adj = exp(log_or_adj), padj_adj = p.adjust(p_adj_model, "BH")), by = category]
  out <- merge(fs, gl[, .(stage, category, group, odds_ratio_adj, p_adj_model, padj_adj, or_abc_per_doubling,
                          or_logCPM, n_model)], by = c("stage", "category", "group"))
  out[, `:=`(group = factor(group, levels = CL), category = factor(category, levels = CAT))]
  setorder(out, stage, category, group)
  fwrite(out, file.path(od, "binding_mode.tsv"), sep = "\t")
  cat("genes in the models:", out$n_model[1], "(genes without a logCPM dropped)\n")
  print(unique(out[, .(stage, category, OR_abc = round(or_abc_per_doubling, 2), OR_logCPM = round(or_logCPM, 2))]))
  print(out[padj < 0.05 | padj_adj < 0.05, .(stage, category, group, frac = round(frac, 3),
                                             bg = round(frac_background, 3), OR = round(odds_ratio, 2),
                                             padj = signif(padj, 2), OR_adj = round(odds_ratio_adj, 2),
                                             padj_adj = signif(padj_adj, 2))])

  # ---- figure: stacked fractions per stage, then the odds ratios (Fisher and adjusted)
  star <- function(p) fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", "")))
  fr <- rbindlist(lapply(DAYS, function(d) g[, lapply(.SD, mean), by = group,
                                             .SDcols = paste(d, CAT[1:3], sep = "|")][, stage := d]), use.names = FALSE)
  setnames(fr, 2:4, CAT[1:3])
  fr <- melt(fr, id.vars = c("group", "stage"), variable.name = "category", value.name = "frac")
  fr[, category := factor(category, levels = rev(CAT[1:3]))]
  tot <- out[category == "all", .(group = factor(as.character(group), levels = rev(GROUPS)), stage,
                                  frac = pmax(frac, frac_background), lab = sprintf("%.0f%%%s", 100 * frac, star(padj)))]
  bgl <- g[group == "background", .(frac = sapply(DAYS, function(d) mean(get(paste(d, "all", sep = "|")))))][,
    stage := DAYS]
  ns <- g[, .N, by = group]; GL <- setNames(sprintf("%s (%s)", ns$group, format(ns$N, big.mark = ",")), ns$group)
  MCOL <- c("TSS only" = "#08306b", "TSS + enhancer" = "#4292c6", "enhancer only" = "#c6dbef")
  pa <- ggplot(fr, aes(frac, group, fill = category)) + geom_col(width = 0.75) +
    geom_vline(data = bgl, aes(xintercept = frac), linetype = "dashed", linewidth = 0.3, colour = "grey30") +
    geom_text(data = tot, aes(x = frac + 0.02, y = group, label = lab), inherit.aes = FALSE, hjust = 0, size = 1.9) +
    facet_grid(. ~ stage) +
    scale_fill_manual(values = MCOL, breaks = CAT[1:3], name = NULL) +
    scale_x_continuous(limits = c(0, 1.12), breaks = c(0, 0.5, 1), labels = c("0", "50%", "100%"), expand = c(0, 0)) +
    scale_y_discrete(labels = GL) +
    labs(x = "genes with an ISL1 site (sites called that day)", y = NULL,
         title = "ISL1 binding of the D15 trajectory gene clusters: TSS and enhancers",
         subtitle = paste0("TSS = site within +/- 1 kb of a TSS; enhancer = site linked only through an ABC element; ",
                           "dashed = background (all);\n* BH FDR < 0.05, ** < 0.01, *** < 0.001, all bound vs background ",
                           "(Fisher)")) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom",
          plot.title = element_text(face = "bold"), plot.title.position = "plot", legend.key.size = unit(0.3, "cm"),
          panel.spacing.x = unit(0.25, "cm"))
  orl <- melt(out[, .(stage, category, group, Fisher = odds_ratio, padj_F = padj, adjusted = odds_ratio_adj,
                      padj_A = padj_adj)], id.vars = c("stage", "category", "group"),
              measure.vars = list(c("Fisher", "adjusted"), c("padj_F", "padj_A")), variable.name = "test",
              value.name = c("or", "padj"))
  orl[, test := factor(c("Fisher", "adjusted")[as.integer(test)], levels = c("Fisher", "adjusted"))]
  orl[, `:=`(l2 = pmax(pmin(log2(or), 1.5), -1.5), lab = star(padj), group = factor(as.character(group), levels = rev(CL)))]
  pb2 <- ggplot(orl, aes(category, group, fill = l2)) + geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 2, vjust = 0.75) +
    facet_grid(test ~ stage) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", limits = c(-1.5, 1.5),
                         name = "log2 odds ratio\nvs background") +
    labs(x = NULL, y = NULL, subtitle = paste0("odds ratio vs background per category. Fisher = Fisher's exact test; ",
                                               "adjusted = logistic regression with\nlog2(1 + ABC elements of the gene) and ",
                                               "mean logCPM as covariates. * BH FDR < 0.05, ** < 0.01, *** < 0.001\n",
                                               "(BH within a category, over cluster x stage)")) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1),
          strip.text = element_text(face = "bold"), legend.key.size = unit(0.3, "cm"))
  p <- pa / pb2 + plot_layout(heights = c(1, 1.5))
  ggsave(file.path(od, "isl1_binding_mode.pdf"), p, width = 6.2, height = 6.4, useDingbats = FALSE)
  cat("wrote", file.path(od, "isl1_binding_mode.pdf"), "\n")
}
