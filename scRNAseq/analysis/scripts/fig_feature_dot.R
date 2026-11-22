# Figures: marker feature plots and dot plots, split by genotype, one day at a time (user, 2026-11-22). Runs locally.
#
# Feature plots (plot/feature_{day}/{gene}.pdf): control | ISL1 KO on the day's SCT B UMAP, each panel the same
#   fixed square as every UMAP (PANEL_CM, fig_umap_labels.R); the other genotype's cells in light grey underneath;
#   cells of the shown genotype from a darker grey (log1p 0) to red, one colour scale per gene for both panels, capped
#   at the 99th percentile of the expressing cells of both genotypes (values above drawn at the cap); expressing
#   cells drawn on top. High-resolution cores (no axes / title / legend): plot/feature_{day}/png/{gene}_{ctr,KO}.png.
# Dot plot (plot/dotplot_{day}.pdf): cell types (label_colours.tsv order) x genes, each gene with its control and KO
#   column side by side (user, 2026-11-22: no grid, no cell numbers); size = % cells with log1p > 0, colour = mean
#   log1p scaled to the gene's maximum over all cell type x genotype groups shown; a group with < MIN_CELLS cells
#   (e.g. KO SC_alpha) is drawn as a small grey cross.
# Inputs (copied from Alpine, git-ignored): reports/labels/{day}/labels.csv, reports/labels/{day}/sct_umaps.csv.gz,
#   reports/figure_data/{day}_genes.csv.gz (fig_export_genes.py)
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_feature_dot.R D22 GCG ARX INS ...
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
PANEL_CM <- 6; POINT_SIZE <- 0.6; MIN_CELLS <- 20; SEED <- 1234
PNG_DPI <- 600   # feature-plot cores: plot/feature_{day}/png/{gene}_{ctr,KO}.png, 6 cm at 600 dpi
OTHER <- "#ececec"                                   # the other genotype
SCALE <- c("#bdbdbd", "#fcbba1", "#ef3b2c", "#99000d")   # darker grey (0) -> red
GENO <- c(ctr = "Control", KO = "ISL1 KO")
a <- commandArgs(TRUE); day <- a[1]; genes <- a[-1]; stopifnot(length(genes) > 0)
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
n_fmt <- function(n) trimws(format(as.integer(n), big.mark = ","))

l <- fread(file.path("reports/labels", day, "labels.csv"))
u <- fread(file.path("reports/labels", day, "sct_umaps.csv.gz")); setnames(u, 1, "barcode")
e <- fread(cmd = sprintf("gzcat reports/figure_data/%s_genes.csv.gz", day))
genes <- intersect(genes, names(e))
x <- merge(merge(l[, .(barcode, genotype, label)], u[, .(barcode, UMAP1 = SCT_B_UMAP1, UMAP2 = SCT_B_UMAP2)]), e,
           by = "barcode")
stopifnot(nrow(x) == nrow(l))
set.seed(SEED); x <- x[sample(.N)]
ng <- x[, .N, by = genotype]
fl <- setNames(sprintf("%s (%s cells)", GENO[ng$genotype], n_fmt(ng$N)), ng$genotype)
x[, panel := factor(fl[genotype], fl[names(GENO)])]
half <- max(diff(range(x$UMAP1)), diff(range(x$UMAP2))) / 2 * 1.08
cx <- mean(range(x$UMAP1)); cy <- mean(range(x$UMAP2))

save_fixed <- function(p, f) {
  g <- ggplotGrob(p)
  pp <- g$layout[grepl("^panel", g$layout$name), ]
  for (i in unique(pp$l)) g$widths[i] <- unit(PANEL_CM, "cm")
  for (i in unique(pp$t)) g$heights[i] <- unit(PANEL_CM, "cm")
  w <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); h <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  pdf(f, width = w, height = h, useDingbats = FALSE); grid.draw(g); invisible(dev.off())
}

# ---- feature plots
fd <- file.path("plot", sprintf("feature_%s", day)); dir.create(fd, showWarnings = FALSE, recursive = TRUE)
for (gn in genes) {
  v <- x[[gn]]
  cap <- max(quantile(v[v > 0], 0.99), 0.1)
  d <- x[, .(UMAP1, UMAP2, genotype, panel, val = pmin(get(gn), cap))]
  bg <- d[, .(UMAP1, UMAP2, panel = factor(fl[setdiff(names(GENO), genotype)], levels(d$panel))), by = genotype]
  setorder(d, val)
  p <- ggplot(d, aes(UMAP1, UMAP2)) +
    geom_point(data = bg, colour = OTHER, size = POINT_SIZE, stroke = 0) +
    geom_point(aes(colour = val), size = POINT_SIZE, stroke = 0) +
    scale_colour_gradientn(colours = SCALE, limits = c(0, cap), name = "log1p",
                           breaks = function(lim) pretty(lim, 3)) +
    scale_x_continuous(limits = cx + c(-half, half), expand = c(0, 0)) +
    scale_y_continuous(limits = cy + c(-half, half), expand = c(0, 0)) +
    coord_fixed() + facet_wrap(~panel, nrow = 1) +
    labs(title = bquote(italic(.(gn)) ~ "(" * .(day) * ")")) +
    theme_classic(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(face = "bold"),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8, hjust = 0),
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.25, "cm"), panel.spacing = unit(0.4, "cm"),
          plot.margin = margin(4, 4, 4, 4))
  save_fixed(p, file.path(fd, sprintf("%s.pdf", gn)))
  # high-resolution cores (user, 2026-11-22): plotting square only, one PNG per genotype, same scale as the PDF;
  # the colour bar is in the PDF (log1p 0 to the cap above)
  for (gt in names(GENO)) {
    pc <- ggplot(d[genotype == gt], aes(UMAP1, UMAP2)) +
      geom_point(data = d[genotype != gt], colour = OTHER, size = POINT_SIZE, stroke = 0) +
      geom_point(aes(colour = val), size = POINT_SIZE, stroke = 0) +
      scale_colour_gradientn(colours = SCALE, limits = c(0, cap)) +
      scale_x_continuous(limits = cx + c(-half, half), expand = c(0, 0)) +
      scale_y_continuous(limits = cy + c(-half, half), expand = c(0, 0)) +
      coord_fixed() + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
    ggsave(file.path(fd, "png", sprintf("%s_%s.png", gn, gt)), pc, width = PANEL_CM, height = PANEL_CM, units = "cm",
           dpi = PNG_DPI, bg = "white", create.dir = TRUE)
  }
}
cat(sprintf("wrote %d feature plots to %s\n", length(genes), fd))

# ---- dot plot
lv <- intersect(names(pal), unique(x$label))
long <- melt(x[, c("label", "genotype", genes), with = FALSE], id.vars = c("label", "genotype"), variable.name = "gene")
dp <- long[, .(pct = 100 * mean(value > 0), mean = mean(value), cells = .N), by = .(label, genotype, gene)]
full <- CJ(label = lv, genotype = names(GENO), gene = genes)
dp <- merge(full, dp[, gene := as.character(gene)], by = c("label", "genotype", "gene"), all.x = TRUE)
dp[is.na(cells), cells := 0L]
dp[, shown := cells >= MIN_CELLS]
dp[, scaled := mean / max(mean[shown]), by = gene]
dp[, `:=`(label = factor(label, rev(lv)), gene = factor(gene, genes),
          g = factor(c(ctr = "Ctr", KO = "KO")[genotype], c("Ctr", "KO")))]
p <- ggplot(dp, aes(g, label)) +
  geom_point(data = dp[shown == TRUE], aes(size = pct, colour = scaled)) +
  geom_point(data = dp[shown == FALSE], shape = 4, size = 1.3, stroke = 0.6, colour = "grey25") +
  scale_size_area(max_size = 3.6, limits = c(0, 100), breaks = c(25, 50, 100), name = "% expressing") +
  scale_colour_gradientn(colours = SCALE, limits = c(0, 1), breaks = c(0, 0.5, 1), name = "mean log1p\n(scaled)") +
  guides(colour = guide_colourbar(barheight = unit(1.4, "cm"), barwidth = unit(0.25, "cm"))) +
  facet_grid(. ~ gene, switch = "x") +
  labs(x = NULL, y = NULL, title = sprintf("%s: marker genes, control (Ctr) and ISL1 KO", day),
       caption = sprintf("\u00d7 = fewer than %d cells of that cell type in that genotype", MIN_CELLS)) +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.placement = "outside", strip.text.x.bottom = element_text(face = "italic", size = 6.5, angle = 40, hjust = 1, vjust = 1), strip.clip = "off",
        axis.text.x = element_text(size = 5, colour = "grey35"), axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0.08, "cm"), plot.caption = element_text(size = 6, colour = "grey40"),
        legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
        legend.key.height = unit(0.3, "cm"), legend.spacing.y = unit(0.1, "cm"))
f <- file.path("plot", sprintf("dotplot_%s.pdf", day))
ggsave(f, p, width = 2.1 + 0.34 * length(genes), height = 1.55 + 0.115 * length(lv), device = cairo_pdf)
cat("wrote", f, "\n")
