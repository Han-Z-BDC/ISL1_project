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
#        Rscript 12_trajectory_binding.R timing <out dir>
# timing (user, 2026-11-24: are the clusters enriched among the ISL1 direct targets, and is it stage-specific): a target
#   = gene with a TSS or ABC-enhancer site (as "bound"). 1 first bound: each gene's first stage with a site (D15 /
#   first at D18 / first at D22 / never); per cluster vs background: chi-square on the 2 x 4 table (class mix), BH
#   over clusters; per class Fisher's exact test (class vs not), BH within a class. 2 interaction: bound (gene x stage,
#   3 rows per gene) ~ group * stage (background, D15 = reference), logistic GEE with an independence working model
#   and gene-clustered robust (sandwich) variance (base R; lme4 / geepack are not installed); Wald tests: all
#   group x stage terms jointly, and per cluster its 2 terms (BH over clusters). A cluster x stage odds ratio > 1 =
#   the cluster's enrichment among targets is larger at that stage than at D15.
#        Rscript 12_trajectory_binding.R within <out dir>
# within (user, 2026-11-24): universe = the C1-C10 genes only (no background). Per stage (target = TSS or ABC-enhancer
#   site called that day) and cluster: is the cluster's share of the targets different from its share of the cluster
#   genes, i.e. 2 x 2 Fisher's exact test, two-sided, cluster vs the other nine clusters x target vs not; BH over
#   10 clusters x 3 stages. Per stage also a 10 x 2 chi-square (does the target fraction differ among the clusters).
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
         subtitle = paste0("TSS = site within +/- 1 kb of a TSS; enhancer = site linked only through an ABC element;\n",
                           "dashed = background (all); * BH FDR < 0.05, ** < 0.01, *** < 0.001, all bound vs background ",
                           "(Fisher)")) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom",
          plot.title = element_text(face = "bold"), plot.title.position = "plot", legend.key.size = unit(0.3, "cm"),
          panel.spacing.x = unit(0.25, "cm"))
  # bar plot only (user, 2026-11-24: "make the bar plot only, keep the current size"): stage panels fixed at the size
  # they had above the odds-ratio tiles (1.32 x 1.58 in); the odds ratios stay in binding_mode.tsv
  ga <- ggplotGrob(pa); pp <- ga$layout[grepl("^panel", ga$layout$name), ]
  for (i in unique(pp$l)) ga$widths[i] <- unit(1.32, "in")
  for (i in unique(pp$t)) ga$heights[i] <- unit(1.58, "in")
  w <- grid::convertWidth(sum(ga$widths), "in", valueOnly = TRUE); h <- grid::convertHeight(sum(ga$heights), "in", valueOnly = TRUE)
  pdf(file.path(od, "isl1_binding_mode.pdf"), width = w, height = h, useDingbats = FALSE, bg = "white")
  grid::grid.draw(ga); invisible(dev.off())
  cat(sprintf("isl1_binding_mode.pdf: %.2f x %.2f in\n", w, h))
  cat("wrote", file.path(od, "isl1_binding_mode.pdf"), "\n")
}

if (mode == "timing") {
  od <- a[2]
  g <- fread(file.path(od, "genes.tsv"))
  GROUPS <- c(sprintf("C%d", 1:10), "background"); CL <- GROUPS[1:10]
  FB <- c("D15", "first at D18", "first at D22", "never")
  g[, first := fifelse(bound_D15, FB[1], fifelse(bound_D18, FB[2], fifelse(bound_D22, FB[3], FB[4])))]
  g[, first := factor(first, levels = FB)]
  fwrite(g[, .(gene, gene_id, group, bound_D15, bound_D18, bound_D22, first)], file.path(od, "genes_timing.tsv"), sep = "\t")
  tab <- dcast(g[, .N, by = .(group, first)], group ~ first, value.var = "N", fill = 0)
  print(tab[order(match(group, GROUPS))])
  bgv <- unlist(tab[group == "background", ..FB])
  # 1a class mix per cluster vs background
  mix <- rbindlist(lapply(CL, function(cl) {
    m <- rbind(unlist(tab[group == cl, ..FB]), bgv)
    ct <- suppressWarnings(chisq.test(m))
    data.table(group = cl, n = sum(m[1, ]), chisq = unname(ct$statistic), df = unname(ct$parameter), p_mix = ct$p.value)
  }))
  mix[, padj_mix := p.adjust(p_mix, "BH")]
  # 1b per class
  cls <- rbindlist(lapply(FB, function(k) rbindlist(lapply(CL, function(cl) {
    v <- g$first == k; inn <- g$group == cl; bg <- g$group == "background"
    ft <- fisher.test(matrix(c(sum(v & inn), sum(!v & inn), sum(v & bg), sum(!v & bg)), 2))
    data.table(class = k, group = cl, n = sum(inn), n_class = sum(v & inn), frac = mean(v[inn]),
               frac_background = mean(v[bg]), odds_ratio = unname(ft$estimate), p = ft$p.value)
  }))))
  cls[, padj := p.adjust(p, "BH"), by = class]
  # 2 interaction: logistic GEE (independence), gene-clustered sandwich variance
  long <- melt(g[, .(gene, group, D15 = bound_D15, D18 = bound_D18, D22 = bound_D22)], id.vars = c("gene", "group"),
               variable.name = "stage", value.name = "bound")
  long[, `:=`(group = relevel(factor(group), ref = "background"), stage = factor(stage, levels = DAYS),
              bound = as.numeric(bound))]
  fit <- glm(bound ~ group * stage, family = binomial, data = long)
  X <- model.matrix(fit); r <- long$bound - fitted(fit)
  B <- solve(crossprod(X * sqrt(fit$weights)))                            # (X'WX)^-1
  U <- rowsum(X * r, long$gene, reorder = FALSE)                          # score per gene
  V <- B %*% crossprod(U) %*% B; b <- coef(fit)
  wald <- function(idx) { bb <- b[idx]; s <- as.numeric(t(bb) %*% solve(V[idx, idx, drop = FALSE]) %*% bb)
                          c(stat = s, df = length(idx), p = pchisq(s, length(idx), lower.tail = FALSE)) }
  ii <- grep(":", names(b))
  ov <- wald(ii)
  cat(sprintf("all cluster x stage terms: Wald chi2 = %.1f, df = %d, p = %.2g (robust, %d genes)\n", ov[1], ov[2], ov[3],
              uniqueN(long$gene)))
  inter <- rbindlist(lapply(CL, function(cl) {
    k <- grep(sprintf("^group%s:", cl), names(b)); w <- wald(k)
    data.table(group = cl, or_D18_vs_D15 = exp(b[sprintf("group%s:stageD18", cl)]),
               or_D22_vs_D15 = exp(b[sprintf("group%s:stageD22", cl)]),
               or_cluster_D15 = exp(b[sprintf("group%s", cl)]), wald_chi2 = w[1], p_interaction = w[3])
  }))
  inter[, padj_interaction := p.adjust(p_interaction, "BH")]
  res <- merge(mix, inter, by = "group")[order(match(group, CL))]
  fwrite(res, file.path(od, "timing_tests.tsv"), sep = "\t")
  fwrite(cls[order(class, match(group, CL))], file.path(od, "timing_classes.tsv"), sep = "\t")
  fwrite(data.table(test = "all cluster x stage terms", wald_chi2 = ov[1], df = ov[2], p = ov[3]),
         file.path(od, "timing_overall.tsv"), sep = "\t")
  print(res[, .(group, n, p_mix = signif(p_mix, 2), padj_mix = signif(padj_mix, 2), OR_D15 = round(or_cluster_D15, 2),
                OR_D18vD15 = round(or_D18_vs_D15, 2), OR_D22vD15 = round(or_D22_vs_D15, 2),
                padj_int = signif(padj_interaction, 2))])
  print(cls[padj < 0.05, .(class, group, frac = round(frac, 3), bg = round(frac_background, 3), OR = round(odds_ratio, 2),
                           padj = signif(padj, 2))])

  # ---- figure: class mix per cluster (stacked) + per-class odds ratios
  star <- function(p) fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", "")))
  fr <- g[, .N, by = .(group, first)][, frac := N / sum(N), by = group]
  fr[, `:=`(group = factor(group, levels = rev(GROUPS)), first = factor(first, levels = rev(FB)))]
  ns <- g[, .N, by = group]
  GL <- setNames(sprintf("%s (%s)", ns$group, format(ns$N, big.mark = ",", trim = TRUE)), ns$group)
  lab <- res[, .(group = factor(group, levels = rev(GROUPS)), lab = sprintf("%s | %s", fifelse(star(padj_mix) == "", "ns",
                 star(padj_mix)), fifelse(star(padj_interaction) == "", "ns", star(padj_interaction))))]
  TCOL <- c("D15" = "#00441b", "first at D18" = "#41ab5d", "first at D22" = "#c7e9c0", "never" = "#e0e0e0")
  p1 <- ggplot(fr, aes(frac, group, fill = first)) + geom_col(width = 0.75) +
    geom_text(data = lab, aes(x = 1.02, y = group, label = lab), inherit.aes = FALSE, hjust = 0, size = 1.9) +
    scale_fill_manual(values = TCOL, breaks = FB, name = "first bound") +
    scale_x_continuous(limits = c(0, 1.32), breaks = c(0, 0.5, 1), labels = c("0", "50%", "100%"), expand = c(0, 0)) +
    scale_y_discrete(labels = GL) +
    labs(x = "genes", y = NULL) +
    theme_classic(base_size = 7) +
    theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"), plot.title = element_text(face = "bold"),
          plot.title.position = "plot")
  cl2 <- cls[, .(class = factor(class, levels = FB), group = factor(group, levels = rev(GROUPS)),
                 l2 = pmax(pmin(log2(odds_ratio), 1.5), -1.5), lab = star(padj))]
  p2 <- ggplot(cl2, aes(class, group, fill = l2)) + geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 2, vjust = 0.75) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", limits = c(-1.5, 1.5),
                         name = "log2 OR\nvs background") +
    scale_y_discrete(drop = FALSE, limits = rev(GROUPS)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_blank(),
          legend.key.size = unit(0.3, "cm"))
  p <- p1 + p2 + plot_layout(widths = c(2.2, 1)) +
    plot_annotation(title = "When the trajectory genes become ISL1 targets",
                    subtitle = paste0("first stage with an ISL1 site (TSS or ABC enhancer). Bar labels: class mix vs ",
                                      "background (chi-square) | change of enrichment\nacross stages (cluster x stage, ",
                                      "robust logistic GEE). Tiles: odds ratio of each class vs background (Fisher, BH ",
                                      "within a class).\nBH FDR * < 0.05, ** < 0.01, *** < 0.001"),
                    theme = theme(plot.title = element_text(face = "bold", size = 8), plot.subtitle = element_text(size = 6.5)))
  ggsave(file.path(od, "isl1_target_timing.pdf"), p, width = 6.4, height = 3.2, useDingbats = FALSE)
  cat("wrote", file.path(od, "isl1_target_timing.pdf"), "\n")
}

if (mode == "within") {
  od <- a[2]
  g <- fread(file.path(od, "genes.tsv"))[group != "background"]
  CL <- sprintf("C%d", 1:10)
  res <- rbindlist(lapply(DAYS, function(d) rbindlist(lapply(CL, function(cl) {
    v <- g[[paste0("bound_", d)]]; inn <- g$group == cl
    ft <- fisher.test(matrix(c(sum(v & inn), sum(!v & inn), sum(v & !inn), sum(!v & !inn)), 2))
    data.table(stage = d, group = cl, n = sum(inn), n_target = sum(v & inn), frac = mean(v[inn]),
               frac_other = mean(v[!inn]), frac_all = mean(v), share_of_targets = sum(v & inn) / sum(v),
               share_of_genes = mean(inn), odds_ratio = unname(ft$estimate), p = ft$p.value)
  }))))
  res[, padj := p.adjust(p, "BH")]
  chi <- rbindlist(lapply(DAYS, function(d) {
    ct <- suppressWarnings(chisq.test(table(g$group, g[[paste0("bound_", d)]])))
    data.table(stage = d, chisq = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value)
  }))
  fwrite(res, file.path(od, "within_clusters.tsv"), sep = "\t"); fwrite(chi, file.path(od, "within_chisq.tsv"), sep = "\t")
  print(chi)
  print(res[, .(stage, group, frac = round(frac, 3), others = round(frac_other, 3), share_t = round(share_of_targets, 3),
                share_g = round(share_of_genes, 3), OR = round(odds_ratio, 2), padj = signif(padj, 2))])

  # figure (user, 2026-11-25: "statistics in this way, but the plot like the previous version with TSS only, TSS and
  # enhancer and enhancer only stacks"): bars stacked by binding mode as isl1_binding_mode.pdf (genes_mode.tsv), label =
  # all bound + stars of the within-cluster Fisher test, dashed = all C1-C10; no background row
  star <- function(p) fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", "")))
  CAT <- c("TSS only", "TSS + enhancer", "enhancer only")
  gm <- fread(file.path(od, "genes_mode.tsv"), check.names = FALSE)[group != "background"]
  fr <- rbindlist(lapply(DAYS, function(d) {
    x <- gm[, lapply(.SD, mean), by = group, .SDcols = paste(d, CAT, sep = "|")]
    setnames(x, c("group", CAT)); x[, stage := d]
  }))
  stopifnot(all.equal(fr[, rowSums(.SD), .SDcols = CAT],
                      res[match(paste(fr$stage, fr$group), paste(stage, group)), frac]))   # same totals as the test
  fr <- melt(fr, id.vars = c("group", "stage"), variable.name = "category", value.name = "frac")
  fr[, category := factor(category, levels = rev(CAT))]
  res[, `:=`(group = factor(group, levels = rev(CL)), lab = sprintf("%.0f%%%s", 100 * frac, star(padj)))]
  avg <- res[, .(frac_all = frac_all[1]), by = stage]
  SL <- setNames(sprintf("%s\nchi-square p = %s", chi$stage, format.pval(chi$p, digits = 2)), chi$stage)
  res[, stage := factor(SL[stage], levels = SL)]; avg[, stage := factor(SL[stage], levels = SL)]
  fr[, stage := factor(SL[stage], levels = SL)]
  fr[, group := factor(group, levels = rev(CL))]
  ns <- g[, .N, by = group]; GL <- setNames(sprintf("%s (%s)", ns$group, ns$N), ns$group)
  MCOL <- c("TSS only" = "#08306b", "TSS + enhancer" = "#4292c6", "enhancer only" = "#c6dbef")   # as isl1_binding_mode.pdf
  p <- ggplot(fr, aes(frac, group, fill = category)) + geom_col(width = 0.75) +
    geom_vline(data = avg, aes(xintercept = frac_all), linetype = "dashed", linewidth = 0.3, colour = "grey30") +
    geom_text(data = res, aes(x = pmax(frac, frac_all) + 0.02, y = group, label = lab), inherit.aes = FALSE,
              hjust = 0, size = 1.9) +
    facet_grid(. ~ stage) +
    scale_fill_manual(values = MCOL, breaks = CAT, name = NULL) +
    scale_x_continuous(limits = c(0, 1.12), breaks = c(0, 0.5, 1), labels = c("0", "50%", "100%"), expand = c(0, 0)) +
    scale_y_discrete(labels = GL) +
    labs(x = "genes with an ISL1 site (sites called that day)", y = NULL,
         title = "ISL1 targets: each cluster vs the other trajectory clusters",
         subtitle = paste0("TSS = site within +/- 1 kb of a TSS; enhancer = site linked only through an ABC element;\n",
                           "dashed = all 3,118 C1-C10 genes; all bound, cluster vs the other nine (Fisher, two-sided),\n",
                           "BH over 10 clusters x 3 stages: * FDR < 0.05, ** < 0.01, *** < 0.001; chi-square: does the\n",
                           "target fraction differ among the ten clusters")) +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom",
          plot.title = element_text(face = "bold"), plot.title.position = "plot", legend.key.size = unit(0.3, "cm"),
          panel.spacing.x = unit(0.25, "cm"))
  ga <- ggplotGrob(p); pp <- ga$layout[grepl("^panel", ga$layout$name), ]
  for (i in unique(pp$l)) ga$widths[i] <- unit(1.32, "in")          # as isl1_binding_mode.pdf
  for (i in unique(pp$t)) ga$heights[i] <- unit(1.58, "in")
  w <- grid::convertWidth(sum(ga$widths), "in", valueOnly = TRUE); h <- grid::convertHeight(sum(ga$heights), "in", valueOnly = TRUE)
  pdf(file.path(od, "isl1_target_within_clusters.pdf"), width = w, height = h, useDingbats = FALSE, bg = "white")
  grid::grid.draw(ga); invisible(dev.off())
  cat(sprintf("isl1_target_within_clusters.pdf: %.2f x %.2f in\n", w, h))
}
