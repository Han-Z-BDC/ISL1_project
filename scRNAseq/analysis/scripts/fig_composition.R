# Figure: label proportions KO vs control per day, from step 7.1l (user, 2026-11-22). Runs locally (R, ggplot2).
#
# 1. plot/D22D40/composition_{day}.pdf: one small panel per label (label order and colours from label_colours.tsv); x =
#    genotype, bar = mean proportion of the genotype, dots = the two replicate samples; propeller FDR on top
#    (*** < 0.001, ** < 0.01, * < 0.05, ns). Free y axis per label (labels differ 1,000-fold in size).
# 2. plot/D22D40/composition_log2ratio.pdf: log2(KO / control) mean proportion per label (0.1% offset), D22 and D40 side by
#    side, filled = propeller FDR < 0.05.
# Inputs (copied from Alpine, git-ignored): reports/composition/{proportions,tests}.csv (7.1l)
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_composition.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
# output folder: plot/D22D40 (user, 2026-11-22: D22 / D40 figures in their own subfolder); override with PLOT_DIR
PLOT <- Sys.getenv("PLOT_DIR", "plot/D22D40"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
pr <- fread("reports/composition/proportions.csv"); te <- fread("reports/composition/tests.csv")
GENO <- c(ctr = "Control", KO = "ISL1 KO")
stars <- function(q) fifelse(q < 0.001, "***", fifelse(q < 0.01, "**", fifelse(q < 0.05, "*", "ns")))
th <- theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7),
        plot.title = element_text(face = "bold"), legend.position = "none")

for (d in unique(te$day)) {
  x <- pr[day == d]; t <- te[day == d]
  lv <- intersect(names(pal), unique(x$label))
  x[, label := factor(label, lv)]; t[, label := factor(label, lv)]
  x[, g := factor(GENO[genotype], GENO)]
  m <- x[, .(proportion = mean(proportion)), by = .(label, g)]
  top <- x[, .(ymax = max(proportion)), by = label]
  t <- merge(t, top, by = "label")[, lab := sprintf("%s\nFDR %s", stars(propeller_FDR), signif(propeller_FDR, 2))]
  p <- ggplot(x, aes(g, 100 * proportion)) +
    geom_col(data = m, aes(fill = label), width = 0.65, alpha = 0.85) +
    geom_point(size = 1.2, shape = 21, fill = "white", colour = "black", stroke = 0.4,
               position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_text(data = t, aes(x = 1.5, y = 100 * ymax * 1.28, label = lab), inherit.aes = FALSE, size = 2.1,
              lineheight = 0.9) +
    geom_blank(data = t, aes(x = 1.5, y = 100 * ymax * 1.5), inherit.aes = FALSE) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    facet_wrap(~label, scales = "free_y", nrow = 2) +
    labs(x = NULL, y = "% of cells in the sample",
         title = sprintf("%s: cell type proportions, control vs ISL1 KO (2 replicates each; propeller)", d)) +
    th + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  nc <- ceiling(length(lv) / 2)
  f <- file.path(PLOT, sprintf("composition_%s.pdf", d))
  ggsave(f, p, width = 1.15 * nc + 0.6, height = 4.2, useDingbats = FALSE)
  cat("wrote", f, "\n")
}

# log2 ratio summary, one PDF per day with only that day's labels (user, 2026-11-22): plot/D22D40/composition_log2ratio_{day}.pdf;
# dot size = propeller FDR in 4 bins (user, 2026-11-22), hollow = not significant
lv <- rev(intersect(names(pal), unique(te$label)))
FB <- c("ns", "< 0.05", "< 0.01", "< 0.001")
te[, label := factor(label, lv)]
te[, fdr_bin := factor(fifelse(propeller_FDR < 0.001, FB[4], fifelse(propeller_FDR < 0.01, FB[3],
                         fifelse(propeller_FDR < 0.05, FB[2], FB[1]))), levels = FB)]
xr <- range(c(0, te$log2_ratio)) + c(-0.5, 0.5)   # same x range on both days
for (d in unique(te$day)) {
  td <- te[day == d][, label := factor(as.character(label), intersect(lv, as.character(label)))]
  p <- ggplot(td, aes(log2_ratio, label)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = log2_ratio, yend = label), colour = "grey70", linewidth = 0.4) +
    geom_point(aes(colour = label, size = fdr_bin, shape = fdr_bin), stroke = 0.7) +
    # one invisible point per FDR bin so the legend always shows all 4 bins (D22 has no ns label)
    geom_point(data = data.table(log2_ratio = 0, label = levels(td$label)[1], fdr_bin = factor(FB, FB)),
               aes(size = fdr_bin, shape = fdr_bin), alpha = 0) +
    scale_size_manual(values = setNames(c(1.6, 1.6, 2.6, 3.8), FB), name = "FDR", breaks = FB) +
    scale_shape_manual(values = setNames(c(1, 16, 16, 16), FB), name = "FDR", breaks = FB) +
    scale_colour_manual(values = pal, guide = "none") +
    guides(size = guide_legend(override.aes = list(alpha = 1, colour = "black"))) +
    scale_x_continuous(limits = xr) +
    labs(x = "log2(KO / control proportion)", y = NULL, title = sprintf("%s: ISL1 KO vs control", d)) +
    th + theme(legend.position = "bottom")
  f <- sprintf(file.path(PLOT, "composition_log2ratio_%s.pdf"), d)
  ggsave(f, p, width = 3.4, height = 0.95 + 0.2 * nlevels(td$label), useDingbats = FALSE)
  cat("wrote", f, "\n")
}

# 3. plot/D22D40/composition_stacked.pdf (user, 2026-11-22): one stacked bar per sample (100% of its cells), coloured by
#    label in the label_colours.tsv order, control then KO replicates, D22 and D40 side by side (each day's own labels).
#    The legend is a table (user): colour key, label, and per day the propeller FDR stars with the KO direction
#    (up / down arrow; "-" = label not present that day).
suppressPackageStartupMessages(library(patchwork))
sm <- copy(pr)
lv_all <- intersect(names(pal), unique(sm$label))
sm[, label := factor(label, rev(lv_all))]
sm[, sample := factor(sprintf("%s %s", GENO[genotype], sub("^.*rep", "rep", sample)),
                      levels = c(sprintf("Control rep%d", 1:2), sprintf("ISL1 KO rep%d", 1:2)))]
bars <- ggplot(sm, aes(sample, 100 * proportion, fill = label)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.1) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~day, nrow = 1) +
  labs(x = NULL, y = "% of cells in the sample", title = "Cell type composition per sample") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1))

days <- sort(unique(te$day))
tab <- CJ(label = lv_all, day = days, sorted = FALSE)
tab <- merge(tab, te[, .(label = as.character(label), day, propeller_FDR, log2_ratio)], by = c("label", "day"),
             all.x = TRUE)
tab[, txt := fifelse(is.na(propeller_FDR), "\u2013",
                     paste0(stars(propeller_FDR), " ", fifelse(log2_ratio > 0, "\u2191", "\u2193")))]
tab[stars(propeller_FDR) == "ns", txt := "ns"]
tab[, y := match(label, rev(lv_all))]
key <- data.table(label = lv_all, y = match(lv_all, rev(lv_all)))
xd <- setNames(2.6 + 0.75 * seq_along(days), days)
ny <- length(lv_all)
star_tab <- ggplot() +
  geom_tile(data = key, aes(x = 0, y = y, fill = label), width = 0.35, height = 0.7) +
  geom_text(data = key, aes(x = 0.3, y = y, label = label), hjust = 0, size = 2.5) +
  geom_text(data = tab, aes(x = xd[day], y = y, label = txt), size = 2.5) +
  annotate("text", x = xd, y = ny + 1, label = names(xd), fontface = "bold", size = 2.6) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(-0.3, max(xd) + 0.45), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, ny + 1.5), expand = c(0, 0)) +
  labs(caption = "KO vs control, propeller FDR\n*** < 0.001, ** < 0.01, * < 0.05; \u2191 / \u2193 KO higher / lower;\n\u2013 not present that day") +
  theme_void(base_size = 8) + theme(plot.caption = element_text(size = 6, hjust = 0, colour = "grey30"))
p <- bars + star_tab + plot_layout(widths = c(1.5, 1))
ggsave(file.path(PLOT, "composition_stacked.pdf"), p, width = 6.2, height = 3.6, device = cairo_pdf)
cat("wrote plot/D22D40/composition_stacked.pdf\n")
