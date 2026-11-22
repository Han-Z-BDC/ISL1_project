# Figure (exploratory): D15 whole-day UMAP for annotation, before the endocrine subset (user, 2026-11-22). Runs locally.
# Same UMAP style as fig_umap_labels.R (fixed PANEL_CM square, equal scales); cluster numbers printed at the
# cluster centres because these plots are for choosing labels.
# Inputs (copied from Alpine, git-ignored): reports/labels/D15/clusters.csv (7.1b Leiden on SCT B, all cells),
#   reports/labels/D15/sct_umaps.csv.gz (SCT B UMAP), reports/labels/D15/endocrine_subclusters.csv (7.1d endocrine
#   subset, sub_r0p6) with study_config/annotation/D15_endocrine_r0.6.tsv, reports/figure_data/D15_genes.csv.gz
#   (fig_export_genes.py), genotypes from cells_obs via the sample name in the barcode.
# Outputs in plot/D15/ (git-ignored):
#   umap_D15_r{res}.pdf           all cells by Leiden cluster (res = 0.2, 0.4)
#   umap_D15_r0.2_genotype.pdf    control | KO
#   umap_D15_endocrine_names.pdf  the 11 named endocrine subclusters on the whole-day UMAP, other cells grey
#   dotplot_D15_r{res}.pdf        marker genes by cluster, control and KO side by side
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_d15_overview.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
PLOT <- Sys.getenv("PLOT_DIR", "plot/D15"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
PANEL_CM <- 6; POINT_SIZE <- 0.6; SEED <- 1234; GREY <- "#e3e3e3"; MIN_CELLS <- 20
GENO <- c(ctr = "Control", KO = "ISL1 KO")
SCALE <- c("#bdbdbd", "#fcbba1", "#ef3b2c", "#99000d")
CL <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        "#393b79", "#e6ab02", "#66a61e", "#a6761d", "#e7298a")
GENES <- c("CHGA", "INSM1", "NEUROG3", "FEV", "ISL1", "INS", "GCG", "ARX", "SST", "GHRL", "TPH1", "PHOX2A", "PDX1",
           "NKX6-1", "ONECUT1", "SOX9", "SPINK1", "DLK1", "KRT19", "PTF1A", "AFP", "APOB", "HNF4A", "CDX2", "COL1A2",
           "COL3A1", "PDGFRA", "CRABP2", "MKI67")
n_fmt <- function(n) trimws(format(as.integer(n), big.mark = ","))

cl <- fread("reports/labels/D15/clusters.csv"); setnames(cl, 1, "barcode")
u <- fread("reports/labels/D15/sct_umaps.csv.gz"); setnames(u, 1, "barcode")
x <- merge(cl, u[, .(barcode, UMAP1 = SCT_B_UMAP1, UMAP2 = SCT_B_UMAP2)], by = "barcode")
x[, genotype := fifelse(grepl("-KO_", barcode), "KO", "ctr")]
stopifnot(nrow(x) == nrow(cl), x[, uniqueN(genotype)] == 2)
sub <- fread("reports/labels/D15/endocrine_subclusters.csv"); setnames(sub, 1, "barcode")
nm <- fread("study_config/annotation/D15_endocrine_r0.6.tsv")
x[sub, endo := nm$label[match(i.sub_r0p6, nm$subcluster)], on = "barcode"]
set.seed(SEED); x <- x[sample(.N)]
half <- max(diff(range(x$UMAP1)), diff(range(x$UMAP2))) / 2 * 1.08
cx <- mean(range(x$UMAP1)); cy <- mean(range(x$UMAP2))
base_th <- theme_classic(base_size = 8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(face = "bold"),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8, hjust = 0),
        legend.key.height = unit(0.32, "cm"), legend.key.width = unit(0.3, "cm"), plot.margin = margin(4, 4, 4, 4),
        panel.spacing = unit(0.4, "cm"))
lims <- list(scale_x_continuous(limits = cx + c(-half, half), expand = c(0, 0)),
             scale_y_continuous(limits = cy + c(-half, half), expand = c(0, 0)), coord_fixed())
save_fixed <- function(p, f) {
  g <- ggplotGrob(p)
  pp <- g$layout[grepl("^panel", g$layout$name), ]
  for (i in unique(pp$l)) g$widths[i] <- unit(PANEL_CM, "cm")
  for (i in unique(pp$t)) g$heights[i] <- unit(PANEL_CM, "cm")
  w <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); h <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  pdf(f, width = w, height = h, useDingbats = FALSE); grid.draw(g); invisible(dev.off())
}

for (res in c("0.2", "0.4")) {
  k <- paste0("leiden_r", sub(".", "p", res, fixed = TRUE))
  x[, cluster := factor(get(k), sort(unique(get(k))))]
  n <- x[, .N, by = cluster]
  cen <- x[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = cluster]
  pal <- setNames(rep(CL, length.out = nlevels(x$cluster)), levels(x$cluster))
  p <- ggplot(x, aes(UMAP1, UMAP2, colour = cluster)) + geom_point(size = POINT_SIZE, stroke = 0) +
    geom_label(data = cen, aes(label = cluster), size = 2.6, fontface = "bold", colour = "black", fill = "white",
               alpha = 0.85, label.size = 0, label.padding = unit(0.08, "lines")) +
    scale_colour_manual(values = pal, name = sprintf("cluster (r%s)", res),
                        labels = function(b) sprintf("%s (%s)", b, n_fmt(n$N[match(b, n$cluster)]))) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) + lims +
    labs(title = sprintf("D15, all cells (%s), Leiden r%s", n_fmt(nrow(x)), res)) + base_th
  save_fixed(p, file.path(PLOT, sprintf("umap_D15_r%s.pdf", res)))

  # dot plot by cluster, control and KO side by side
  e <- fread(cmd = "gzcat reports/figure_data/D15_genes.csv.gz")
  gs <- intersect(GENES, names(e))
  long <- melt(merge(x[, .(barcode, cluster, genotype)], e, by = "barcode")[, c("cluster", "genotype", gs), with = FALSE],
               id.vars = c("cluster", "genotype"), variable.name = "gene")
  dp <- long[, .(pct = 100 * mean(value > 0), mean = mean(value), cells = .N), by = .(cluster, genotype, gene)]
  dp <- merge(CJ(cluster = levels(x$cluster), genotype = names(GENO), gene = gs),
              dp[, `:=`(cluster = as.character(cluster), gene = as.character(gene))], all.x = TRUE,
              by = c("cluster", "genotype", "gene"))
  dp[is.na(cells), cells := 0L][, shown := cells >= MIN_CELLS][, scaled := mean / max(mean[shown]), by = gene]
  dp[, `:=`(cluster = factor(cluster, rev(levels(x$cluster))), gene = factor(gene, gs),
            g = factor(c(ctr = "Ctr", KO = "KO")[genotype], c("Ctr", "KO")))]
  p <- ggplot(dp, aes(g, cluster)) +
    geom_point(data = dp[shown == TRUE], aes(size = pct, colour = scaled)) +
    geom_point(data = dp[shown == FALSE], shape = 4, size = 1.3, stroke = 0.6, colour = "grey25") +
    scale_size_area(max_size = 3.6, limits = c(0, 100), breaks = c(25, 50, 100), name = "% expressing") +
    scale_colour_gradientn(colours = SCALE, limits = c(0, 1), breaks = c(0, 0.5, 1), name = "mean log1p\n(scaled)") +
    guides(colour = guide_colourbar(barheight = unit(1.4, "cm"), barwidth = unit(0.25, "cm"))) +
    facet_grid(. ~ gene, switch = "x") +
    labs(x = NULL, y = sprintf("cluster (r%s)", res), title = sprintf("D15 Leiden r%s: marker genes, control and KO", res),
         caption = sprintf("× = fewer than %d cells of that cluster in that genotype", MIN_CELLS)) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(face = "bold"), strip.background = element_blank(), strip.placement = "outside",
          strip.text.x.bottom = element_text(face = "italic", size = 6.5, angle = 40, hjust = 1, vjust = 1),
          strip.clip = "off", axis.text.x = element_text(size = 5, colour = "grey35"), axis.ticks.x = element_blank(),
          panel.spacing.x = unit(0.08, "cm"), plot.caption = element_text(size = 6, colour = "grey40"),
          legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
          legend.key.height = unit(0.3, "cm"))
  ggsave(file.path(PLOT, sprintf("dotplot_D15_r%s.pdf", res)), p, width = 2.3 + 0.34 * length(gs),
         height = 1.6 + 0.15 * nlevels(x$cluster), device = cairo_pdf)
}

# genotype split at r0.2
x[, cluster := factor(leiden_r0p2, sort(unique(leiden_r0p2)))]
pal <- setNames(rep(CL, length.out = nlevels(x$cluster)), levels(x$cluster))
ng <- x[, .N, by = genotype]; fl <- setNames(sprintf("%s (%s cells)", GENO[ng$genotype], n_fmt(ng$N)), ng$genotype)
x[, panel := factor(fl[genotype], fl[names(GENO)])]
bg <- x[, .(UMAP1, UMAP2, panel = factor(fl[setdiff(names(GENO), genotype)], levels(x$panel))), by = genotype]
p <- ggplot(x, aes(UMAP1, UMAP2)) + geom_point(data = bg, colour = GREY, size = POINT_SIZE, stroke = 0) +
  geom_point(aes(colour = cluster), size = POINT_SIZE, stroke = 0) + facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = pal, name = "cluster (r0.2)") +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) + lims + labs(title = "D15, Leiden r0.2") + base_th
save_fixed(p, file.path(PLOT, "umap_D15_r0.2_genotype.pdf"))

# endocrine subcluster names on the whole-day UMAP
lvn <- nm$label
x[, en := factor(endo, lvn)]
np <- setNames(c("#ff7f0e", "#bcbd22", "#8c6d31", "#17becf", "#d62728", "#1f77b4", "#8c564b", "#9467bd", "#fdae6b",
                 "#6baed6", "#000000"), lvn)
n <- x[!is.na(en), .N, by = en]
p <- ggplot(x, aes(UMAP1, UMAP2)) + geom_point(data = x[is.na(en)], colour = GREY, size = POINT_SIZE, stroke = 0) +
  geom_point(data = x[!is.na(en)], aes(colour = en), size = POINT_SIZE, stroke = 0) +
  scale_colour_manual(values = np, name = "endocrine subcluster (7.1d)",
                      labels = function(b) sprintf("%s (%s)", b, n_fmt(n$N[match(b, n$en)]))) +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) + lims +
  labs(title = sprintf("D15: named endocrine subclusters (%s cells); other cells grey", n_fmt(x[!is.na(en), .N]))) +
  base_th
save_fixed(p, file.path(PLOT, "umap_D15_endocrine_names.pdf"))
cat("wrote", PLOT, "\n")
