# Figure: one gene along the D15 beta, alpha and EC lineages, control vs ISL1 KO (user, 2026-11-24: "single gene, like
# the average curve plots, coloured as the D15 endocrine UMAP, with lines where the KO cells stop for alpha and beta").
# Curves as 7.3h / gene_dynamics (fig_d15_trajectory.R figure 8): Slingshot lineage cells with curve weight >= 0.5,
#   10 equal-width pseudotime bins per lineage, curve-weighted mean log1p per genotype, 3-bin moving average, no point
#   where a genotype has < 15 cells in the bin. Style as gene_patterns_curves.pdf (fig_d15_gene_patterns.R): colours =
#   the UMAP colours of pre_beta2 / pre_alpha / pre_EC, control solid / KO dashed, red dashed line = KO stop on beta /
#   alpha (KO 95th percentile of pseudotime, curve weight >= 0.5, as ko_stop_density.pdf), panels 0.55 x 0.72 in.
# Inputs (git-ignored): reports/trajectory/D15/slingshot_pc30_{cells,lineages}.csv (7.3a),
#   reports/figure_data/D15_genes.csv.gz (fig_export_genes.py D15; a gene must be exported there first).
# Output: plot/D15/trajectory/genes/{gene}_curve.pdf
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_d15_gene_curves.R ARX [more genes]
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
PLOT <- file.path(Sys.getenv("PLOT_DIR", "plot/D15/trajectory"), "genes"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
TR <- "reports/trajectory/D15"
GENES <- commandArgs(TRUE); stopifnot(length(GENES) > 0)
GENO <- c(ctr = "Control", KO = "ISL1 KO")
LIN <- c(pre_beta2 = "beta", pre_alpha = "alpha", pre_EC = "EC")
NBIN <- 10; MIN_CELLS <- 15; W_MIN <- 0.5; CURVE_W <- 0.55; CURVE_H <- 0.72
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
LCOL <- setNames(unname(pal[names(LIN)]), LIN)
KO_COL <- "#d6604d"

sl <- fread(file.path(TR, "slingshot_pc30_cells.csv")); setnames(sl, 1, "barcode")
ln <- fread(file.path(TR, "slingshot_pc30_lineages.csv"))
lend <- setNames(sapply(strsplit(ln$clusters, " -> "), function(p) p[length(p)]), ln$lineage)
e <- fread(cmd = "gzcat reports/figure_data/D15_genes.csv.gz")
miss <- setdiff(GENES, names(e)); if (length(miss)) cat("not in D15_genes.csv.gz (export first), skipped:", miss, "\n")
GENES <- intersect(GENES, names(e))
e <- e[match(sl$barcode, barcode)]; stopifnot(!anyNA(e$barcode))

stops <- rbindlist(lapply(c("pre_beta2", "pre_alpha"), function(end) {
  L <- names(lend)[lend == end]; w <- sl[[paste0("w_", L)]]; pt <- sl[[paste0("pt_", L)]]
  m <- !is.na(w) & w >= W_MIN & !is.na(pt) & sl$genotype == "KO"
  data.table(lineage = factor(LIN[[end]], levels = LIN), q = unname(quantile(pt[m], 0.95)))
}))

for (gn in GENES) {
  cur <- rbindlist(lapply(names(LIN), function(end) {
    L <- names(lend)[lend == end]
    pt <- sl[[paste0("pt_", L)]]; w <- sl[[paste0("w_", L)]]; w[is.na(w)] <- 0
    on <- w >= W_MIN & !is.na(pt)
    br <- seq(min(pt[on]), max(pt[on]), length.out = NBIN + 1)
    b <- pmin(pmax(findInterval(pt, br, rightmost.closed = TRUE), 1), NBIN)
    rbindlist(lapply(names(GENO), function(gt) {
      r <- rbindlist(lapply(seq_len(NBIN), function(i) {
        m <- on & sl$genotype == gt & b == i
        data.table(bin = i, mid = (br[i] + br[i + 1]) / 2, n = sum(m),
                   mean = if (any(m)) weighted.mean(e[[gn]][m], w[m]) else NA_real_)
      }))
      v <- ifelse(r$n >= MIN_CELLS, r$mean, NA)
      r[, smooth := sapply(seq_len(.N), function(i) { k <- max(1, i - 1):min(.N, i + 1)
                                                      if (is.na(v[i])) NA else mean(v[k], na.rm = TRUE) })]
      r[, `:=`(genotype = gt, lineage = LIN[[end]])]
    }))
  }))
  cur[, lineage := factor(lineage, levels = LIN)]
  p <- ggplot(cur, aes(mid, smooth, colour = lineage, linetype = genotype)) +
    geom_vline(data = stops, aes(xintercept = q), linetype = "dashed", colour = KO_COL, linewidth = 0.35) +
    geom_line(linewidth = 0.6, na.rm = TRUE) +
    facet_grid(. ~ lineage, scales = "free_x") +
    scale_colour_manual(values = LCOL, guide = "none") +
    scale_linetype_manual(values = c(ctr = "solid", KO = "22"), breaks = names(GENO), labels = GENO, name = NULL) +
    scale_x_continuous(breaks = c(0, 100, 200)) +
    labs(x = "Slingshot pseudotime", y = "mean log1p", title = gn,
         subtitle = "control solid, ISL1 KO dashed;\nred dashed line = KO stop") +
    theme_classic(base_size = 7) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold.italic"), plot.title.position = "plot", legend.position = "top",
          legend.key.width = unit(0.5, "cm"), panel.spacing.x = unit(0.1, "cm"))
  g <- ggplotGrob(p); pp <- g$layout[grepl("^panel", g$layout$name), ]
  for (i in unique(pp$l)) g$widths[i] <- unit(CURVE_W, "in")
  for (i in unique(pp$t)) g$heights[i] <- unit(CURVE_H, "in")
  w <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); h <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  f <- file.path(PLOT, sprintf("%s_curve.pdf", gn))
  pdf(f, width = w, height = h, useDingbats = FALSE, bg = "white"); grid.draw(g); invisible(dev.off())
  cat(sprintf("wrote %s (%.2f x %.2f in)\n", f, w, h))
}
