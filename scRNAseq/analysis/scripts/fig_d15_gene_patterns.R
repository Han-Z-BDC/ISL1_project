# Figure: D15 trajectory gene patterns, the final 3,118 KO vs control trajectory genes in k = 10 clusters (7.3h on the
# 7.3l pseudobulk + 7.3m tradeSeq set; user, 2026-11-24: "drop the sparse start bins and build the figure").
# Nothing is recomputed: curves and clusters are read from the 7.3h outputs.
# Inputs (reports/trajectory/D15 -> ${OUT_DIR}/day_analysis/D15/trajectory):
#   pseudobulk/gene_patterns/gene_pattern_curves.csv   per gene x lineage x genotype x bin: cells, smoothed (curve-
#                                                      weighted mean log1p, 3-bin moving average; NA when < 15 cells)
#   pseudobulk/gene_patterns/gene_patterns_k10.csv     cluster per gene
#   study_config/annotation/TF_names_v_1.01.txt        Lambert 2018 human TFs (TF count per cluster)
# Columns: 10 equal-width pseudotime bins per lineage (beta, alpha, EC) x genotype. Bins where both genotypes have
#   < 15 cells (the sparse curve starts: bin 1 on every lineage, bin 2 on EC) are dropped. Remaining bins with < 15
#   cells of one genotype are grey: "no KO cells" after the KO stop on beta / alpha (bins 8-10), "< 15 cells" for
#   control EC bin 10. Values = z-score per gene across all kept bins of the six curves (as the clustering input).
# Rows: clusters C1-C10 (7.3h numbering), genes within a cluster ordered by the column of their maximum z (stable).
# Labelled genes: the user's lineage-gene list (fig_d15_trajectory.R GENES) that are in the set.
# Outputs in plot/D15/trajectory/ (git-ignored):
#   gene_patterns_heatmap.pdf   heatmap, cluster strip (genes, TFs), labelled genes at the right; heatmap area 2.5 in wide
#   gene_patterns_curves.pdf    mean z per cluster along beta / alpha / EC (control solid, KO dashed), IQR ribbon;
#                               panels 0.55 x 0.72 in
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_d15_gene_patterns.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
PLOT <- Sys.getenv("PLOT_DIR", "plot/D15/trajectory"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
GP <- "reports/trajectory/D15/pseudobulk/gene_patterns"
MIN_CELLS <- 15; ZLIM <- 2.5
GENO <- c(ctr = "Control", KO = "ISL1 KO")
LIN <- c("beta", "alpha", "EC")
LABEL_GENES <- c("ISL1", "PDX1", "FOXO1", "PBX3", "MEIS3", "PAX6", "LMO2", "NKX6-1", "MNX1", "ARX", "ETV1", "PAX4",
                 "MAFB", "NEUROD1", "NEUROD4", "LHX1", "FEV")
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
LCOL <- c(beta = unname(pal["SC_beta"]), alpha = unname(pal["SC_alpha"]), EC = unname(pal["SC_EC1"]))
n_fmt <- function(n) trimws(format(as.integer(n), big.mark = ","))

cl <- fread(file.path(GP, "gene_patterns_k10.csv"), select = c("gene", "cluster"))
cv <- fread(file.path(GP, "gene_pattern_curves.csv"), select = c("gene", "lineage", "genotype", "bin", "pseudotime",
                                                                 "cells", "smoothed"))
cv <- cv[gene %in% cl$gene]
stopifnot(nrow(cl) == 3118, uniqueN(cv$gene) == nrow(cl), setequal(unique(cv$lineage), LIN))

# ---- bins: drop those sparse in both genotypes; grey reason for the rest
nb <- unique(cv[, .(lineage, genotype, bin, pseudotime, cells)])
keep <- nb[, .(keep = any(cells >= MIN_CELLS)), by = .(lineage, bin)]
print(dcast(nb, lineage + genotype ~ bin, value.var = "cells"))
cat("dropped bins (both genotypes < 15 cells):", keep[!(keep), paste(lineage, bin)], "\n")
cv <- merge(cv, keep, by = c("lineage", "bin"))[(keep)][, keep := NULL]

# ---- z per gene across all kept bins of the six curves (population sd, as 7.3h)
cv[, z := { v <- smoothed; m <- mean(v, na.rm = TRUE); s <- sqrt(mean((v - m)^2, na.rm = TRUE)); (v - m) / (s + 1e-9) },
   by = gene]
cv <- merge(cv, cl, by = "gene")
cv[, curve := factor(sprintf("%s %s", lineage, genotype), levels = as.vector(t(outer(LIN, names(GENO), paste))))]
setorder(cv, curve, bin); cv[, col := rleid(curve, bin)]

# gene order: cluster, then the column of the maximum z
pk <- cv[!is.na(z), .(peak = col[which.max(z)]), by = .(gene, cluster)]
setorderv(pk, c("cluster", "peak", "gene"))
pk[, row := .N - seq_len(.N) + 1]          # C1 at the top
cv <- merge(cv, pk[, .(gene, row)], by = "gene")

# ---- heatmap
tf <- fread("study_config/annotation/TF_names_v_1.01.txt", header = FALSE)$V1
cs <- pk[, .(n = .N, tfs = sum(gene %in% tf), top = max(row), bottom = min(row)), by = cluster]
CLAB <- setNames(sprintf("C%d\n%s genes\n%d TFs", cs$cluster, n_fmt(cs$n), cs$tfs), cs$cluster)
cv[, clab := factor(CLAB[as.character(cluster)], levels = CLAB[order(cs$cluster)])]
cv[, zc := pmax(pmin(z, ZLIM), -ZLIM)]
grey <- cv[is.na(z), .(gene, curve, bin, row, clab, reason = fifelse(genotype == "KO" & lineage != "EC" & cells == 0,
                                                                    "no KO cells", "< 15 cells"))]
CURVE_LAB <- setNames(sprintf("%s\n%s", sub(" .*", "", levels(cv$curve)), GENO[sub(".* ", "", levels(cv$curve))]),
                      levels(cv$curve))
GREY_NOKO <- "#9e9e9e"; GREY_FEW <- "#d9d9d9"
hm <- ggplot(cv, aes(bin, row)) +
  geom_raster(data = cv[!is.na(z)], aes(fill = zc)) +
  geom_raster(data = grey[reason == "no KO cells"], fill = GREY_NOKO) +
  geom_raster(data = grey[reason == "< 15 cells"], fill = GREY_FEW) +
  scale_fill_gradientn(colours = c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"), limits = c(-ZLIM, ZLIM),
                       name = "z-score") +
  facet_grid(clab ~ curve, scales = "free", space = "free", switch = "y", labeller = labeller(curve = CURVE_LAB)) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) + scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  labs(x = "pseudotime (Slingshot, 10 equal-width bins;\nsparse start bins dropped)", y = NULL,
       title = sprintf("D15 KO vs control trajectory genes (%s), k-means k = 10", n_fmt(nrow(cl))),
       subtitle = "z-score per gene across the six curves;\ndark grey = no KO cells (past the KO stop), light grey = < 15 cells") +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), panel.spacing.x = unit(0.06, "cm"), panel.spacing.y = unit(0.04, "cm"),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6), strip.text.x = element_text(face = "bold", size = 6),
        strip.placement = "outside", plot.title = element_text(face = "bold"), plot.title.position = "plot",
        legend.position = "bottom", legend.key.height = unit(0.25, "cm"), legend.key.width = unit(0.6, "cm"))

# labelled genes at the right, one panel per cluster with the same row scale
lab <- pk[gene %in% LABEL_GENES]
cat("labelled:", paste(lab$gene, lab$cluster, collapse = "; "), "| not in set:",
    paste(setdiff(LABEL_GENES, lab$gene), collapse = ", "), "\n")
lab[, clab := factor(CLAB[as.character(cluster)], levels = levels(cv$clab))]
setorder(lab, -row)
lab[, ylab := {   # spread labels inside the cluster so they do not overlap
  mn <- min(pk$row[pk$cluster == cluster[1]]); mx <- max(pk$row[pk$cluster == cluster[1]])
  gap <- max((mx - mn) * 0.16, 24); y <- row
  for (i in seq_along(y)[-1]) y[i] <- min(y[i], y[i - 1] - gap)
  lo <- mn + 0.06 * (mx - mn); hi <- mx - 0.04 * (mx - mn)   # keep labels inside the cluster band
  if (min(y) < lo) y <- y + (lo - min(y))
  pmin(y, hi)
}, by = cluster]
lb <- ggplot(lab) +
  geom_blank(data = cs[, .(clab = factor(CLAB[as.character(cluster)], levels = levels(cv$clab)), top, bottom)],
             aes(y = top)) +
  geom_blank(data = cs[, .(clab = factor(CLAB[as.character(cluster)], levels = levels(cv$clab)), top, bottom)],
             aes(y = bottom)) +
  geom_segment(aes(x = 0, xend = 0.35, y = row, yend = ylab), linewidth = 0.25, colour = "grey40") +
  geom_text(aes(x = 0.4, y = ylab, label = gene), hjust = 0, size = 2.1, fontface = "italic") +
  facet_grid(clab ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 7) + theme(strip.text = element_blank(), panel.spacing.y = unit(0.04, "cm"))
# narrow columns (user, 2026-11-24: "make the heatmap columns narrower too"): heatmap area HM_W wide, label column LAB_W
HM_W <- 2.5; LAB_W <- 0.5   # inches
p <- hm + lb + plot_layout(widths = unit(c(HM_W, LAB_W), "in"))
f <- file.path(PLOT, "gene_patterns_heatmap.pdf")
ggsave(f, p, width = HM_W + LAB_W + 0.9, height = 8.2, useDingbats = FALSE); cat("wrote", f, "\n")

# ---- mean curves per cluster
mc <- cv[!is.na(z), .(mean = mean(z), q1 = quantile(z, 0.25), q3 = quantile(z, 0.75)),
         by = .(cluster, lineage, genotype, bin, pseudotime)]
mc[, `:=`(lineage = factor(lineage, levels = LIN), clab = factor(sprintf("C%d (%s)", cluster, n_fmt(cs$n[match(cluster, cs$cluster)])),
                                                                  levels = sprintf("C%d (%s)", cs$cluster[order(cs$cluster)], n_fmt(cs$n[order(cs$cluster)]))))]
pc <- ggplot(mc, aes(pseudotime, mean, colour = lineage, fill = lineage, linetype = genotype)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_ribbon(aes(ymin = q1, ymax = q3, group = genotype), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.55) +
  facet_grid(clab ~ lineage, scales = "free_x", switch = "y") +
  scale_colour_manual(values = LCOL, guide = "none") + scale_fill_manual(values = LCOL, guide = "none") +
  scale_linetype_manual(values = c(ctr = "solid", KO = "22"), breaks = names(GENO), labels = GENO, name = NULL) +
  labs(x = "Slingshot pseudotime", y = "mean z-score (ribbon = interquartile range)",
       title = "D15 trajectory gene clusters: mean curves",
       subtitle = "control solid, ISL1 KO dashed;\nno line where a genotype has < 15 cells") +
  scale_x_continuous(breaks = c(0, 100, 200)) +
  theme_classic(base_size = 7) +
  theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold"),
        strip.text.x = element_text(face = "bold"), strip.placement = "outside", plot.title = element_text(face = "bold"),
        plot.title.position = "plot", legend.position = "top", panel.spacing.y = unit(0.15, "cm"),
        panel.spacing.x = unit(0.1, "cm"))
f <- file.path(PLOT, "gene_patterns_curves.pdf")
# fixed panel size (user, 2026-11-24: the three lineages together as wide as one lineage was before, ~1.6 in)
CURVE_W <- 0.55; CURVE_H <- 0.72   # inches per panel
g <- ggplotGrob(pc); pp <- g$layout[grepl("^panel", g$layout$name), ]
for (i in unique(pp$l)) g$widths[i] <- unit(CURVE_W, "in")
for (i in unique(pp$t)) g$heights[i] <- unit(CURVE_H, "in")
w <- grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE); h <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
pdf(f, width = w, height = h, useDingbats = FALSE); grid::grid.draw(g); invisible(dev.off())
cat(sprintf("wrote %s: %.2f x %.2f in\n", f, w, h))
