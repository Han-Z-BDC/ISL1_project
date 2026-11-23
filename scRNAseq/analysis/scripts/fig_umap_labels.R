# Figure: UMAPs of one day coloured by the final cell labels (user, 2026-11-22). Runs locally (R, ggplot2).
#
# Style (user): no labels on the UMAP (legend only); bigger dots; every plotting region is a square of fixed size,
# identical for every UMAP (PANEL_CM), with equal x / y scales (the data range is centred in the square); colours from
# study_config/annotation/label_colours.tsv (one colour per label name, shared by all days), legend in that order.
# Split by genotype (user, 2026-11-22): control and KO side by side on the same axes; on each panel the other
# genotype's cells are drawn light grey underneath.
# Inputs (copied from Alpine, git-ignored): reports/labels/{day}/labels.csv (7.1h) and reports/labels/{day}/
#   sct_umaps.csv.gz (7.1c, SCT B UMAP = the UMAP of the labels).
# Output (git-ignored): plot/D22D40/umap_{day}.pdf (all cells), plot/D22D40/umap_{day}_genotype.pdf (control | KO)
# High-resolution UMAP cores (user, 2026-11-22): the plotting square only (no axes, title, legend), PANEL_CM x PANEL_CM
#   at PNG_DPI, white background: plot/D22D40/umap_core/umap_{day}.png, umap_{day}_ctr.png, umap_{day}_KO.png
# D15 endocrine compartment (user, 2026-11-22): "day" D15_endocrine = the 7.1d subset (5,994 cells) on its own 2D
#   UMAP, labels = r0.6 names (study_config/annotation/D15_endocrine_r0.6.tsv) from reports/labels/D15/
#   endocrine_subclusters.csv; colours: pre_alpha / pre_alpha_KO, pre_beta1 / pre_beta2 and pre_EC are lighter shades
#   of SC_alpha, SC_beta and SC_EC1 (pre_delta: lighter SC_delta); the ENP states and KO_only get distinct hues.
#   Legend in ENDO_ORDER (lineage order). PLOT_DIR=plot/D15 Rscript analysis/scripts/fig_umap_labels.R D15_endocrine
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_umap_labels.R D22 D40
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
# output folder: plot/D22D40 (user, 2026-11-22: D22 / D40 figures in their own subfolder); override with PLOT_DIR
PLOT <- Sys.getenv("PLOT_DIR", "plot/D22D40"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
PANEL_CM <- 6          # side of the square plotting region, all UMAPs
POINT_SIZE <- 0.6
GREY <- "#e3e3e3"      # the other genotype on the split UMAPs
GENO <- c(ctr = "Control", KO = "ISL1 KO")
SEED <- 1234
PNG_DPI <- 600         # UMAP core PNGs: 6 cm at 600 dpi = 1417 x 1417 px
days <- commandArgs(TRUE); stopifnot(length(days) > 0)
pal <- fread("study_config/annotation/label_colours.tsv")
pal <- setNames(pal$colour, pal$label)
n_fmt <- function(n) trimws(format(as.integer(n), big.mark = ","))

ENDO_ORDER <- c("ENP1", "ENP2", "ENP_FEV_high", "ENP_PHOX2A+", "pre_EC", "pre_beta1", "pre_beta2", "pre_alpha",
                "pre_alpha_KO", "pre_delta", "KO_only")
load_endocrine <- function() {
  s <- fread("reports/labels/D15/endocrine_subclusters.csv"); setnames(s, 1, "barcode")
  nm <- fread("study_config/annotation/D15_endocrine_r0.6.tsv")
  x <- s[, .(barcode, UMAP1, UMAP2, label = nm$label[match(sub_r0p6, nm$subcluster)],
             genotype = fifelse(grepl("-KO_", barcode), "KO", "ctr"))]
  stopifnot(nrow(x) == 5994, !anyNA(x$label), setequal(x$label, ENDO_ORDER), all(ENDO_ORDER %in% names(pal)))
  x[, label := factor(label, levels = ENDO_ORDER)]
  set.seed(SEED); x[sample(.N)]
}
load_day <- function(day) {
  if (day == "D15_endocrine") return(load_endocrine())
  l <- fread(file.path("reports/labels", day, "labels.csv"))
  u <- fread(file.path("reports/labels", day, "sct_umaps.csv.gz")); setnames(u, 1, "barcode")
  x <- merge(l[, .(barcode, genotype, label)], u[, .(barcode, UMAP1 = SCT_B_UMAP1, UMAP2 = SCT_B_UMAP2)],
             by = "barcode")
  stopifnot(nrow(x) == nrow(l), all(x$label %in% names(pal)), all(x$genotype %in% names(GENO)))
  x[, label := factor(label, levels = intersect(names(pal), unique(label)))]
  set.seed(SEED); x[sample(.N)]
}

# shared look: square equal-scale limits (longer data range + 4% per side, the other range centred), legend only
umap_base <- function(x, p) {
  half <- max(diff(range(x$UMAP1)), diff(range(x$UMAP2))) / 2 * 1.08
  cx <- mean(range(x$UMAP1)); cy <- mean(range(x$UMAP2))
  n <- table(x$label)
  p + scale_colour_manual(values = pal, name = NULL, drop = TRUE,
                          labels = function(b) sprintf("%s (%s)", b, n_fmt(n[b]))) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    scale_x_continuous(limits = cx + c(-half, half), expand = c(0, 0)) +
    scale_y_continuous(limits = cy + c(-half, half), expand = c(0, 0)) +
    coord_fixed() +
    theme_classic(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(face = "bold"),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8, hjust = 0),
          legend.key.height = unit(0.32, "cm"), legend.key.width = unit(0.3, "cm"), plot.margin = margin(4, 4, 4, 4))
}

# UMAP core: the panel only (same limits, point size and physical size as the PDF panels)
core_png <- function(p, f) {
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  ggsave(f, p + labs(title = NULL) + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0)),
         width = PANEL_CM, height = PANEL_CM, units = "cm", dpi = PNG_DPI, bg = "white")
}

# fix every panel to PANEL_CM x PANEL_CM and size the page to the whole figure
save_fixed <- function(p, f) {
  g <- ggplotGrob(p)
  pp <- g$layout[grepl("^panel", g$layout$name), ]
  for (i in unique(pp$l)) g$widths[i] <- unit(PANEL_CM, "cm")
  for (i in unique(pp$t)) g$heights[i] <- unit(PANEL_CM, "cm")
  w <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); h <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  pdf(f, width = w, height = h, useDingbats = FALSE); grid.draw(g); invisible(dev.off())
  cat(sprintf("%s: page %.2f x %.2f in, %d panel(s) of %g cm\n", f, w, h, nrow(pp), PANEL_CM))
}

for (day in days) {
  x <- load_day(day)
  p <- ggplot(x, aes(UMAP1, UMAP2, colour = label)) + geom_point(size = POINT_SIZE, stroke = 0) +
    labs(title = sprintf("%s (%s cells)", sub("_", " ", day), n_fmt(nrow(x))))
  save_fixed(umap_base(x, p), file.path(PLOT, sprintf("umap_%s.pdf", day)))
  core_png(umap_base(x, p), file.path(PLOT, "umap_core", sprintf("umap_%s.png", day)))
  for (gt in names(GENO)) {   # genotype cores: this genotype coloured, the other grey underneath
    pc <- ggplot(x, aes(UMAP1, UMAP2)) +
      geom_point(data = x[genotype != gt], colour = GREY, size = POINT_SIZE, stroke = 0) +
      geom_point(data = x[genotype == gt], aes(colour = label), size = POINT_SIZE, stroke = 0)
    core_png(umap_base(x, pc), file.path(PLOT, "umap_core", sprintf("umap_%s_%s.png", day, gt)))
  }

  # split: one facet per genotype; the other genotype's cells as grey background in each facet
  ng <- x[, .N, by = genotype]
  fl <- setNames(sprintf("%s (%s cells)", GENO[ng$genotype], n_fmt(ng$N)), ng$genotype)
  x[, panel := factor(fl[genotype], levels = fl[names(GENO)])]
  bg <- x[, .(UMAP1, UMAP2, panel = factor(fl[setdiff(names(GENO), genotype)], levels = levels(x$panel))),
          by = genotype]
  p <- ggplot(x, aes(UMAP1, UMAP2)) +
    geom_point(data = bg, colour = GREY, size = POINT_SIZE, stroke = 0) +
    geom_point(aes(colour = label), size = POINT_SIZE, stroke = 0) +
    facet_wrap(~panel, nrow = 1) +
    labs(title = sub("_", " ", day))
  save_fixed(umap_base(x, p) + theme(panel.spacing = unit(0.4, "cm")),
             file.path(PLOT, sprintf("umap_%s_genotype.pdf", day)))
}
