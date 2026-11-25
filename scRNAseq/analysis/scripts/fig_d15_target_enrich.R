# Figure: GO / pathway enrichment of the direct ISL1 targets per D15 trajectory cluster, paper version (user,
# 2026-11-24: "focus on top 5-8 terms across all categories in each gene cluster and plot them"). Runs locally (R).
# Input (copied from Alpine, git-ignored): ../integration/results/12_trajectory_binding/target_enrich.csv
#   (integration step 12, 12_target_tf_enrich.R: fgsea::fora, MSigDB 2026.1 GO BP / MF / CC, KEGG Medicus / Legacy,
#   Reactome, Hallmark; universe = the 14,490 genes tested along the D15 lineages; BH per cluster x collection).
# Term choice per cluster: terms with padj < 0.05 from all collections pooled, by padj; a term is skipped when its
#   target genes overlap an already chosen term's by Jaccard >= 0.5 (removes near-duplicate GO terms); up to N_TOP.
# Dot plot as fig_ec_deseq2.R / fig_ec_beta.R: columns = clusters, rows = chosen terms [database], dot where padj < 0.05
#   (size = target genes in the term, colour = -log10 padj capped at 10).
# Output (git-ignored): plot/D15/trajectory/isl1_target_enrich_top.pdf; the chosen terms in
#   ../integration/results/12_trajectory_binding/target_enrich_top.csv
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_d15_target_enrich.R [N_TOP, default 6]
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
IN <- "../integration/results/12_trajectory_binding"
PLOT <- Sys.getenv("PLOT_DIR", "plot/D15/trajectory"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
a <- commandArgs(TRUE); N_TOP <- if (length(a)) as.integer(a[1]) else 6L
JAC <- 0.5
COLL <- c("GO BP", "GO MF", "GO CC", "KEGG Medicus", "KEGG Legacy", "Reactome", "Hallmark")
CCOL <- c(`GO BP` = "#1b7837", `GO MF` = "#5aae61", `GO CC` = "#a6dba0", `KEGG Medicus` = "#762a83",
          `KEGG Legacy` = "#9970ab", Reactome = "#2166ac", Hallmark = "#b2182b")
CL <- paste0("C", 1:10)

r <- fread(file.path(IN, "target_enrich.csv"))[padj < 0.05]
pick <- rbindlist(lapply(CL, function(cl) {
  x <- r[cluster == cl][order(padj, -overlap)]
  keep <- integer(0); gl <- list()
  for (i in seq_len(nrow(x))) {
    g <- strsplit(x$genes[i], ",")[[1]]
    if (all(vapply(gl, function(h) length(intersect(g, h)) / length(union(g, h)) < JAC, logical(1)))) {
      keep <- c(keep, i); gl[[length(gl) + 1]] <- g
    }
    if (length(keep) == N_TOP) break
  }
  x[keep]
}))
short <- function(x) {
  x <- sub("^(GOBP|GOMF|GOCC|KEGG_MEDICUS|KEGG|REACTOME|HALLMARK)_", "", x)
  x <- tolower(gsub("_", " ", x))
  x <- gsub("\\bgpcr\\b", "GPCR", gsub("\\bigf\\b", "IGF", gsub("\\becm\\b", "ECM", gsub("\\bg2m\\b", "G2M", x))))
  x <- gsub("\\be2f\\b", "E2F", gsub("\\bkras\\b", "KRAS", gsub("\\ber\\b", "ER", x)))
  ifelse(nchar(x) > 52, paste0(substr(x, 1, 50), "..."), x)
}
pick[, label := short(term)]
pick[, cluster := factor(cluster, CL)]
pick[, collection := factor(collection, COLL)]
pick[, rank := seq_len(.N), by = cluster]
pick[, y := factor(paste(cluster, rank), levels = rev(paste(rep(CL, each = N_TOP), seq_len(N_TOP))))]
fwrite(pick[, .(cluster, rank, collection, term, overlap, targets, term_size, fold, padj, genes)],
       file.path(IN, "target_enrich_top.csv"))
print(pick[, .(cluster, collection, label, overlap, padj = signif(padj, 2))])

# dot plot as the GO dot plots of fig_ec_deseq2.R / fig_ec_beta.R (user, 2026-11-24): columns = clusters, rows = the
# chosen terms (grouped by the cluster that chose them, then rank), a dot wherever the term has padj < 0.05 in a
# cluster; size = target genes in the term, colour = -log10 padj (capped at CAP so C9's cilium terms do not flatten
# the scale); the database in brackets after each term
CAP <- 10
ord <- unique(pick[order(cluster, rank), .(collection, term, label)], by = c("collection", "term"))
ord[, row := sprintf("%s  [%s]", label, collection)]
x <- r[paste(collection, term) %in% paste(ord$collection, ord$term)]
x <- merge(x, ord, by = c("collection", "term"))
x[, row := factor(row, rev(ord$row))][, cluster := factor(cluster, CL)]
nt <- r[, .(targets = targets[1]), by = cluster]
tg <- fread(file.path(IN, "target_TF_counts.csv"))   # targets per cluster (also clusters without any term)
CLAB <- setNames(sprintf("%s (%d)", tg$cluster, tg$targets), tg$cluster)
p <- ggplot(x, aes(cluster, row)) +
  geom_point(aes(size = overlap, colour = pmin(-log10(padj), CAP))) +
  scale_colour_gradient(low = "#fdbb84", high = "#7f0000", limits = c(-log10(0.05), CAP),
                        name = sprintf("-log10 padj\n(capped at %d)", CAP)) +
  scale_size_area(max_size = 3.2, name = "target genes\nin term") +
  labs(x = NULL, y = NULL, title = "Direct ISL1 targets: GO and pathway enrichment",
       subtitle = sprintf(paste0("Top %d terms per cluster, pooled over GO BP / MF / CC, KEGG,\nReactome and Hallmark ",
                                 "(padj < 0.05; near-duplicates, Jaccard >= %.1f, removed).\nDots wherever padj < 0.05. ",
                                 "Columns: cluster (direct targets).\nUniverse: 14,490 genes tested along the D15 lineages."),
                          N_TOP, JAC)) +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 6), plot.title.position = "plot",
        axis.text.y = element_text(size = 6.3), axis.text.x = element_text(face = "bold"),
        panel.grid = element_blank(), legend.key.size = unit(0.3, "cm"))
# compact core (user, 2026-11-24: "make the dots close and the core plot smaller: remove the grids, very small gap
# between rows and columns"): fixed cell size CELL_W x CELL_H in, axes expanded by half a cell only
CELL_W <- 0.15; CELL_H <- 0.115
p <- p + scale_y_discrete(expand = expansion(add = 0.6)) + theme(axis.text.x = element_text(face = "bold", size = 6, angle = 90, hjust = 1, vjust = 0.5))
p <- p + scale_x_discrete(drop = FALSE, labels = CLAB, expand = expansion(add = 0.6))
g <- ggplotGrob(p); pp <- g$layout[grepl("^panel", g$layout$name), ]
g$widths[pp$l] <- unit(CELL_W * (length(CL) + 0.2), "in"); g$heights[pp$t] <- unit(CELL_H * (nlevels(x$row) + 0.2), "in")
w <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); h <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
f <- file.path(PLOT, "isl1_target_enrich_top.pdf")
pdf(f, width = w, height = h, useDingbats = FALSE, bg = "white"); grid.draw(g); invisible(dev.off())
cat(sprintf("wrote %s (%.2f x %.2f in; core %.2f x %.2f in)\n", f, w, h, CELL_W * (length(CL) + 0.2),
            CELL_H * (nlevels(x$row) + 0.2)))
