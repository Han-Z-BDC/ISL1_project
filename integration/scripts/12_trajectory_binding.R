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
#   modes). Rows: clusters C1-C10, then 1,000 background genes (random, seed 1234) as reference; within a block,
#   genes sorted by the mean signal over +/- 500 bp and the three days (highest at the top).
# usage: Rscript 12_trajectory_binding.R prep <gene_patterns_k10.csv> <pb_curve.csv> <links.tsv> <genes.gtf> <out dir>
#        Rscript 12_trajectory_binding.R plot <out dir>
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
  ns <- s[, .N, by = group]
  BL <- setNames(sprintf("%s\n%s", ns$group, format(ns$N, big.mark = ",")), as.character(ns$group))
  BL["background"] <- sprintf("background\n%s random", format(ns[group == "background", N], big.mark = ","))
  hm[, glab := factor(BL[as.character(group)], levels = BL[GROUPS])]
  HM_W <- 2.4; HM_H <- 4.2
  p3 <- ggplot(hm, aes(pos, row, fill = sc)) + geom_raster() +
    facet_grid(glab ~ day, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradientn(colours = c("white", "#fee391", "#fe9929", "#cc4c02", "#662506"), limits = c(0, CAP),
                         name = "ISL1\nsignal") +
    scale_x_continuous(breaks = c(-2000, 0, 2000), labels = c("-2 kb", "TSS", "+2 kb"), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    labs(x = NULL, y = NULL, title = "ISL1 CUT&RUN at the TSS\nof the trajectory gene clusters",
         subtitle = sprintf(paste0("per-day mean of 3 libraries (background-normalised CPM,\n50 bp bins); genes sorted ",
                                   "by mean signal over +/- %d bp\nand the 3 days; colour capped at the 99th percentile"), CORE)) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(), panel.spacing.x = unit(0.08, "cm"), panel.spacing.y = unit(0.04, "cm"),
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 5.5, lineheight = 0.9),
          strip.text.x = element_text(face = "bold"), strip.placement = "outside", axis.text.x = element_text(size = 5),
          plot.title = element_text(face = "bold"), plot.title.position = "plot", legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.25, "cm"))
  g3 <- ggplotGrob(p3); pp <- g3$layout[grepl("^panel", g3$layout$name), ]
  hts <- ns[match(GROUPS, as.character(group)), N]
  for (j in seq_along(unique(pp$t))) g3$heights[sort(unique(pp$t))[j]] <- unit(HM_H * hts[j] / sum(hts), "in")
  for (i in unique(pp$l)) g3$widths[i] <- unit(HM_W / 3, "in")
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
