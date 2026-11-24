# Figures: D15 endocrine trajectory, Slingshot (map, where KO stops) and CellRank (fates, end states, alpha / beta bias),
# plus gene dynamics along the lineages (user, 2026-11-23: figures 1-8 from the existing outputs; nothing rerun).
# Style as fig_umap_labels.R: fixed square UMAP panels (PANEL_CM, equal scales), label colours from
# study_config/annotation/label_colours.tsv in lineage order (ENDO_ORDER), control | KO with the other genotype grey.
# Inputs (copied from Alpine, git-ignored):
#   reports/labels/D15/endocrine_subclusters.csv   7.1d subset (sub_r0p6, UMAP1/2) + D15_endocrine_r0.6.tsv names
#   reports/trajectory/D15/                        ${OUT_DIR}/day_analysis/D15/trajectory: slingshot_pc30_{cells,
#     lineages,curves_2d,tests}.csv (7.3a), cellrank_pseudotime{,_ctr,_KO}_{fates,macrostates}.csv and
#     cellrank_pseudotime_fate_summary.csv (7.3d, pseudotime kernel = main), alpha_bias/D15/alpha_bias_cells_pseudotime.csv
#     (7.3i)
#   reports/figure_data/D15_genes.csv.gz           (fig_export_genes.py D15; log1p, all D15 cells)
# Outputs in plot/D15/trajectory/ (git-ignored):
#   1 slingshot_umap.pdf          lineage curves on the endocrine UMAP, cells by label, control | KO
#     slingshot_umap_all.pdf      the same, all cells on one panel (not split by genotype)
#   2 lineage_tree.pdf            subclusters at their UMAP centres (size = cells), Slingshot lineage edges, arrows
#   3b ko_stop_celltypes.pdf   cells per pseudotime bin stacked by subcluster, per lineage, control / KO
#   3 ko_stop_density.pdf         pseudotime density on each lineage (curve weight >= 0.5), control vs KO; condiments
#                                 progression statistic per lineage; dashed = KO 95th percentile
#   4 cellrank_fates_umap.pdf     fate probability per end state (rows) on the UMAP, control | KO; nodes + arrows of
#                                 that end state's Slingshot path per genotype (dashed circle = < 20 cells)
#   5 cellrank_fate_heatmap.pdf   mean fate probability per subcluster x genotype
#   6 cellrank_end_states.pdf     data-driven macrostates of the all-cell, control-only and KO-only runs (terminal filled)
#   7 alpha_bias_umap.pdf, alpha_bias_hist.pdf   P(alpha) / (P(alpha) + P(beta)) of ENP1, ENP2, ENP_FEV_high,
#                                 pre_alpha_KO (7.3i); x = P(alpha) + P(beta) < 0.05
#   png/ UMAP cores (no labels): slingshot_umap_{ctr,KO,all}, lineage_tree, cellrank_fate_{end}_{ctr,KO},
#     alpha_bias_umap_{ctr,KO}
#   8 gene_dynamics.pdf           binned curves along the beta, alpha and EC lineages, control (solid) vs KO (dashed);
#                                 as 7.3h: curve weight >= 0.5, 10 equal-width bins, curve-weighted mean log1p, 3-bin
#                                 moving average, no line where a genotype has < 15 cells in the bin
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_d15_trajectory.R [GENES for figure 8, comma list]
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })
PLOT <- Sys.getenv("PLOT_DIR", "plot/D15/trajectory"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
TR <- "reports/trajectory/D15"
PANEL_CM <- 6; POINT_SIZE <- 1.1; GREY <- "#e3e3e3"; SEED <- 1234
# POINT_SIZE 1.1 (was 0.6; user, 2026-11-23: dots larger, 5,994 cells looked too sparse; as fig_umap_labels.R D15_endocrine)
GENO <- c(ctr = "Control", KO = "ISL1 KO")
GCOL <- c(ctr = "#2166ac", KO = "#d6604d")
ENDO_ORDER <- c("ENP1", "ENP2", "ENP_FEV_high", "ENP_PHOX2A+", "pre_EC", "pre_beta1", "pre_beta2", "pre_alpha",
                "pre_alpha_KO", "pre_delta", "KO_only")
END_NAME <- c(pre_beta2 = "beta", pre_alpha = "alpha", pre_EC = "EC", "ENP_PHOX2A+" = "PHOX2A+", KO_only = "KO_only",
              pre_delta = "delta")
GENES <- c("ISL1", "PDX1", "FOXO1", "PBX3", "MEIS3", "PAX6", "LMO2", "NKX6-1", "MNX1", "ARX", "ETV1", "PAX4",
           "MAFB", "NEUROD1", "NEUROD4", "LHX1", "FEV")
a <- commandArgs(TRUE); if (length(a)) GENES <- strsplit(a[1], ",")[[1]]
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
n_fmt <- function(n) trimws(format(as.integer(n), big.mark = ","))
FATE_SCALE <- scale_colour_gradientn(colours = c("#d9d9d9", "#fdbb84", "#e34a33", "#7f0000"), limits = c(0, 1),
                                     name = "fate\nprobability")
BIAS_SCALE <- scale_colour_gradientn(colours = c("#08519c", "#7b3fa0", "#cb181d"), limits = c(0, 1),
                                     name = "alpha bias\n0 = beta\n1 = alpha")

# ---- cells: UMAP, label, genotype
s <- fread("reports/labels/D15/endocrine_subclusters.csv"); setnames(s, 1, "barcode")
nm <- fread("study_config/annotation/D15_endocrine_r0.6.tsv")
x <- s[, .(barcode, UMAP1, UMAP2, label = nm$label[match(sub_r0p6, nm$subcluster)],
           genotype = fifelse(grepl("-KO_", barcode), "KO", "ctr"))]
stopifnot(nrow(x) == 5994, !anyNA(x$label), all(ENDO_ORDER %in% names(pal)))
x[, label := factor(label, levels = ENDO_ORDER)]
set.seed(SEED); x <- x[sample(.N)]
ng <- x[, .N, by = genotype]; FL <- setNames(sprintf("%s (%s cells)", GENO[ng$genotype], n_fmt(ng$N)), ng$genotype)
x[, panel := factor(FL[genotype], levels = FL[names(GENO)])]
bg <- x[, .(UMAP1, UMAP2, panel = factor(FL[setdiff(names(GENO), genotype)], levels = levels(x$panel))), by = genotype]

half <- max(diff(range(x$UMAP1)), diff(range(x$UMAP2))) / 2 * 1.08
cx <- mean(range(x$UMAP1)); cy <- mean(range(x$UMAP2))
umap_theme <- list(
  scale_x_continuous(limits = cx + c(-half, half), expand = c(0, 0)),
  scale_y_continuous(limits = cy + c(-half, half), expand = c(0, 0)), coord_fixed(), theme_classic(base_size = 8),
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(face = "bold"),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8, hjust = 0),
        legend.key.height = unit(0.32, "cm"), legend.key.width = unit(0.3, "cm"), plot.margin = margin(4, 4, 4, 4),
        panel.spacing = unit(0.4, "cm")))
label_scale <- function() {
  n <- table(x$label)
  list(scale_colour_manual(values = pal, name = NULL, labels = function(b) sprintf("%s (%s)", b, n_fmt(n[b]))),
       guides(colour = guide_legend(override.aes = list(size = 2.5))))
}
save_fixed <- function(p, f) {   # every panel PANEL_CM x PANEL_CM (as fig_umap_labels.R)
  g <- ggplotGrob(p)
  pp <- g$layout[grepl("^panel", g$layout$name), ]
  for (i in unique(pp$l)) g$widths[i] <- unit(PANEL_CM, "cm")
  for (i in unique(pp$t)) g$heights[i] <- unit(PANEL_CM, "cm")
  w <- convertWidth(sum(g$widths), "in", valueOnly = TRUE); h <- convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  pdf(f, width = w, height = h, useDingbats = FALSE); grid.draw(g); invisible(dev.off())
  cat(sprintf("%s: page %.2f x %.2f in, %d panel(s)\n", f, w, h, nrow(pp)))
}
# UMAP cores (user, 2026-11-23: "png of the UMAP core without any labels"): one PNG per panel, the plotting square only
# (no axes, title, legend, strip or text labels; curves / tree arrows kept), PANEL_CM x PANEL_CM at PNG_DPI, white
# background, same limits and point size as the PDF panels, in plot/D15/trajectory/png/
PNG_DPI <- 600
core_png <- function(layers, f) {
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  p <- ggplot() + layers + umap_theme + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  ggsave(f, p, width = PANEL_CM, height = PANEL_CM, units = "cm", dpi = PNG_DPI, bg = "white")
}
PNG <- file.path(PLOT, "png")

# ---- Slingshot inputs
sl <- fread(file.path(TR, "slingshot_pc30_cells.csv")); setnames(sl, 1, "barcode")
ln <- fread(file.path(TR, "slingshot_pc30_lineages.csv"))
paths <- setNames(strsplit(ln$clusters, " -> "), ln$lineage)
lend <- sapply(paths, function(p) p[length(p)])
lname <- setNames(sprintf("%s (%s)", END_NAME[lend], sub("Lineage", "L", names(lend))), names(lend))
cv <- fread(file.path(TR, "slingshot_pc30_curves_2d.csv"))
stopifnot(nrow(cv) > 0, setequal(unique(cv$lineage), names(paths)))

# 1. lineage curves on the UMAP
ce <- cv[, .SD[.N], by = lineage]
ce[, name := END_NAME[lend[lineage]]]
p <- ggplot(x, aes(UMAP1, UMAP2)) +
  geom_point(data = bg, colour = GREY, size = POINT_SIZE, stroke = 0) +
  geom_point(aes(colour = label), size = POINT_SIZE, stroke = 0) +
  geom_path(data = cv, aes(group = lineage), colour = "black", linewidth = 0.55,
            arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_label(data = ce, aes(label = name), size = 2.2, label.size = 0, label.padding = unit(0.06, "lines"),
             fill = "white", alpha = 0.85, nudge_y = 0.35) +
  facet_wrap(~panel, nrow = 1) + label_scale() + umap_theme +
  labs(title = "D15 endocrine: Slingshot lineages (PCs 1-30, start ENP1, fixed ends)")
save_fixed(p, file.path(PLOT, "slingshot_umap.pdf"))
for (gt in names(GENO)) core_png(list(
  geom_point(data = x[genotype != gt], aes(UMAP1, UMAP2), colour = GREY, size = POINT_SIZE, stroke = 0),
  geom_point(data = x[genotype == gt], aes(UMAP1, UMAP2, colour = label), size = POINT_SIZE, stroke = 0),
  geom_path(data = cv, aes(UMAP1, UMAP2, group = lineage), colour = "black", linewidth = 0.55,
            arrow = arrow(length = unit(0.12, "cm"), type = "closed")),
  scale_colour_manual(values = pal)), file.path(PNG, sprintf("slingshot_umap_%s.png", gt)))
# 1b. same, all cells together (user, 2026-11-23: "full UMAP, not split by genotype", plus PNG core)
p <- ggplot(x, aes(UMAP1, UMAP2)) +
  geom_point(aes(colour = label), size = POINT_SIZE, stroke = 0) +
  geom_path(data = cv, aes(group = lineage), colour = "black", linewidth = 0.55,
            arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_label(data = ce, aes(label = name), size = 2.2, label.size = 0, label.padding = unit(0.06, "lines"),
             fill = "white", alpha = 0.85, nudge_y = 0.35) +
  label_scale() + umap_theme +
  labs(title = sprintf("D15 endocrine: Slingshot lineages, all cells (%s)", n_fmt(nrow(x))))
save_fixed(p, file.path(PLOT, "slingshot_umap_all.pdf"))
core_png(list(
  geom_point(data = x, aes(UMAP1, UMAP2, colour = label), size = POINT_SIZE, stroke = 0),
  geom_path(data = cv, aes(UMAP1, UMAP2, group = lineage), colour = "black", linewidth = 0.55,
            arrow = arrow(length = unit(0.12, "cm"), type = "closed")),
  scale_colour_manual(values = pal)), file.path(PNG, "slingshot_umap_all.png"))

# 2. lineage tree: subclusters at their UMAP medians, edges along the Slingshot paths
cen <- x[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), n = .N), by = label]
ed <- unique(rbindlist(lapply(paths, function(p) data.table(from = head(p, -1), to = p[-1]))))
ed <- merge(merge(ed, cen[, .(from = label, x0 = UMAP1, y0 = UMAP2)], by = "from"),
            cen[, .(to = label, x1 = UMAP1, y1 = UMAP2)], by = "to")
ed[, `:=`(x1s = x0 + 0.82 * (x1 - x0), y1s = y0 + 0.82 * (y1 - y0))]
# background cells #d0d0d0 and nodes with a white outline (user, 2026-11-23: background "more obvious"; chose the
# darker grey + outlined nodes of four previewed options)
TREE_BG <- "#d0d0d0"
node <- function() geom_point(data = cen, aes(UMAP1, UMAP2, fill = label, size = n), shape = 21, colour = "white",
                              stroke = 0.8)
nt <- table(x$label)
p <- ggplot() +
  geom_point(data = x, aes(UMAP1, UMAP2), colour = TREE_BG, size = POINT_SIZE, stroke = 0) +
  geom_segment(data = ed, aes(x0, y0, xend = x1s, yend = y1s), linewidth = 0.6, colour = "grey30",
               arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
  node() +
  geom_label(data = cen, aes(UMAP1, UMAP2, label = label), size = 2.1, vjust = -0.9, label.size = 0,
             label.padding = unit(0.05, "lines"), fill = "white", alpha = 0.8) +
  scale_size_area(max_size = 7, guide = "none") +
  scale_fill_manual(values = pal, name = NULL, labels = function(b) sprintf("%s (%s)", b, n_fmt(nt[b]))) +
  guides(fill = guide_legend(override.aes = list(size = 2.5, stroke = 0.4))) + umap_theme +
  labs(title = "Slingshot lineage tree", subtitle = "subclusters at their UMAP centres, size = cells")
save_fixed(p, file.path(PLOT, "lineage_tree.pdf"))
core_png(list(
  geom_point(data = x, aes(UMAP1, UMAP2), colour = TREE_BG, size = POINT_SIZE, stroke = 0),
  geom_segment(data = ed, aes(x0, y0, xend = x1s, yend = y1s), linewidth = 0.6, colour = "grey30",
               arrow = arrow(length = unit(0.14, "cm"), type = "closed")),
  node(), scale_size_area(max_size = 7), scale_fill_manual(values = pal)), file.path(PNG, "lineage_tree.png"))

# 3. where KO stops: pseudotime density per lineage, control vs KO
te <- fread(file.path(TR, "slingshot_pc30_tests.csv"))
pr <- te[test == "progression" & lineage != "All"]
pr[, lin_id := paste0("Lineage", lineage)]
prs <- setNames(pr$statistic, pr$lin_id); prp <- setNames(pr$p.value, pr$lin_id)
dd <- rbindlist(lapply(names(paths), function(L) {
  sl[get(paste0("w_", L)) >= 0.5 & !is.na(get(paste0("pt_", L))),
     .(lineage = L, genotype, pt = get(paste0("pt_", L)))]
}))
lab_l <- sapply(names(paths), function(L) sprintf("%s: %s\nprogression stat %.2f, p %s", lname[L],
  paste(tail(paths[[L]], 2), collapse = " > "), prs[[L]], format.pval(prp[[L]], digits = 2)))
dd[, lin := factor(lab_l[lineage], levels = lab_l)]
q95 <- dd[genotype == "KO", .(q = quantile(pt, 0.95)), by = lin]
cnt <- dd[, .N, by = .(lin, genotype)]
cnt[, txt := sprintf("%s %s", GENO[genotype], n_fmt(N))]
cnt <- cnt[, .(txt = paste(txt, collapse = "\n")), by = lin]
dall <- te[test == "differentiation" & pair == "All"]
p <- ggplot(dd, aes(pt, fill = genotype, colour = genotype)) +
  geom_density(alpha = 0.3, adjust = 0.8, linewidth = 0.4) +
  geom_vline(data = q95, aes(xintercept = q), linetype = "dashed", colour = GCOL[["KO"]], linewidth = 0.4) +
  geom_text(data = cnt, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = -0.05, vjust = 1.2, size = 2.1) +
  facet_wrap(~lin, ncol = 3, scales = "free") +
  scale_fill_manual(values = GCOL, labels = GENO, name = NULL) + scale_colour_manual(values = GCOL, labels = GENO, name = NULL) +
  labs(x = "Slingshot pseudotime (cells with curve weight >= 0.5)", y = "density",
       title = "Where the KO stops: cells along each lineage",
       subtitle = sprintf("condiments differentiation test, all lineages: stat %.2f, p %s; dashed = KO 95th percentile",
                          dall$statistic, format.pval(dall$p.value, digits = 2))) +
  theme_classic(base_size = 8) + theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 7),
                                       plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(PLOT, "ko_stop_density.pdf"), p, width = 7.8, height = 5.2, useDingbats = FALSE)
cat("wrote", file.path(PLOT, "ko_stop_density.pdf"), "\n")

# 3b. cell types along each lineage (user, 2026-11-23: "distribution of cell types along the pseudotime for ctrl and
# KO"; option 1 of three previewed): same cells as 3 (curve weight >= 0.5, so the shared ENP trunk is in every
# lineage), Slingshot pseudotime in 20 equal-width bins per lineage, bars = cells per bin stacked by subcluster,
# control (top) and KO (bottom)
NB <- 20
dl <- rbindlist(lapply(names(paths), function(L) {
  sl[get(paste0("w_", L)) >= 0.5 & !is.na(get(paste0("pt_", L))), .(barcode, lineage = L, genotype, pt = get(paste0("pt_", L)))]
}))
dl <- merge(dl, x[, .(barcode, label)], by = "barcode")
dl[, `:=`(lo = min(pt), w = (max(pt) - min(pt)) / NB), by = lineage]
dl[, bin := pmin(NB, floor((pt - lo) / w) + 1)][, mid := lo + (bin - 0.5) * w]
cb <- dl[, .N, by = .(lineage, genotype, bin, mid, w, label)]
cb[, lin := factor(lname[lineage], lname[names(paths)])][, geno := factor(GENO[genotype], GENO)]
# stop lines (user, 2026-11-23): where the KO cells stop on the beta, alpha and delta lineages and where the control
# cells stop on the EC and KO_only lineages = 95th percentile of that genotype's pseudotime (as the dashed line in 3),
# drawn across both rows, coloured by genotype
STOP_BY <- c(beta = "KO", alpha = "KO", delta = "KO", EC = "ctr", KO_only = "ctr")
sp <- rbindlist(lapply(names(paths), function(L) {
  g <- STOP_BY[END_NAME[[lend[[L]]]]]; if (is.na(g)) return(NULL)
  data.table(lineage = L, genotype = g, q = quantile(dl[lineage == L & genotype == g, pt], 0.95))
}))
sp[, lin := factor(lname[lineage], levels(cb$lin))]
p <- ggplot(cb, aes(mid, N, fill = label)) + suppressWarnings(geom_col(aes(width = w * 0.95))) +   # width is used
  geom_vline(data = sp, aes(xintercept = q, colour = genotype), linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = GCOL, breaks = names(GENO), labels = c(ctr = "control stops (95th pct)", KO = "KO stops (95th pct)"),
                      name = NULL) +
  facet_grid(geno ~ lin, scales = "free_x") + scale_fill_manual(values = pal, name = NULL) +
  labs(x = "Slingshot pseudotime (cells with curve weight >= 0.5)", y = "cells per bin",
       title = "Cell types along each lineage, control vs KO",
       subtitle = sprintf("%d equal-width pseudotime bins per lineage; the shared ENP trunk appears in every lineage", NB)) +
  theme_classic(base_size = 8) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
                                       plot.title = element_text(face = "bold"), legend.key.size = unit(0.3, "cm"))
ggsave(file.path(PLOT, "ko_stop_celltypes.pdf"), p, width = 11, height = 4.2, useDingbats = FALSE)
cat("wrote", file.path(PLOT, "ko_stop_celltypes.pdf"), "\n")

# 4. CellRank fates on the UMAP (pseudotime kernel, run on all cells)
fa <- fread(file.path(TR, "cellrank_pseudotime_fates.csv")); setnames(fa, 1, "barcode")
ends <- intersect(names(END_NAME), names(fa))
fl <- melt(merge(x[, .(barcode, UMAP1, UMAP2, genotype)], fa[, c("barcode", ends), with = FALSE], by = "barcode"),
           measure.vars = ends, variable.name = "fate", value.name = "p")
fl[, fate := factor(sprintf("%s (%s)", END_NAME[as.character(fate)], fate), levels = sprintf("%s (%s)", END_NAME[ends], ends))]
fl[, gpanel := factor(GENO[genotype], levels = GENO)]
fbg <- fl[, .(UMAP1, UMAP2, gpanel = factor(GENO[setdiff(names(GENO), genotype)], levels = GENO)), by = .(genotype, fate)]
setorder(fl, p)
# lineage nodes on the fate panels (user, 2026-11-24): only the Slingshot path ending at that panel's end state, per
# genotype: node = that genotype's subcluster at its median UMAP position, area = that genotype's cells (one scale for
# all panels), label colour with a white outline, dark arrows with a white outline (user, 2026-11-24: bigger and
# clearer when small: nodes max size 10, dashed circles 1 mm line); a subcluster with < MIN_NODE cells of
# that genotype (KO: pre_alpha, pre_beta2, pre_delta = 0; control: KO_only = 11) = dashed empty circle at the other
# genotype's position, with a dashed arrow into it
MIN_NODE <- 20
gc <- x[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), n = .N), by = .(label, genotype)]
NODE_MAX <- gc[, max(n)]
lin_of <- setNames(names(lend), lend)   # end subcluster -> lineage
nodes <- rbindlist(lapply(ends, function(e) rbindlist(lapply(names(GENO), function(gt) {
  pth <- paths[[lin_of[[e]]]]
  rbindlist(lapply(seq_along(pth), function(i) {
    own <- gc[label == pth[i] & genotype == gt]; oth <- gc[label == pth[i] & genotype != gt]
    ok <- nrow(own) == 1 && own$n >= MIN_NODE
    src <- if (ok) own else oth
    data.table(fate = sprintf("%s (%s)", END_NAME[[e]], e), genotype = gt, step = i, label = pth[i],
               UMAP1 = src$UMAP1, UMAP2 = src$UMAP2, n = if (ok) own$n else 0L, present = ok)
  }))
}))))
nodes[, `:=`(fate = factor(fate, levels(fl$fate)), gpanel = factor(GENO[genotype], levels = GENO))]
R_MISS <- 2.4 / 60 * 2 * half   # dashed circle radius: 2.4 mm on the 6 cm panel
NODE_SIZE_MAX <- 10
MM <- 2 * half / (PANEL_CM * 10)   # data units per mm on the panel
# arrows run from edge to edge of the nodes (node radius ~ 0.41 mm per size unit + outline, 0.5 mm gap)
nodes[, r := fifelse(present, (0.41 * NODE_SIZE_MAX * sqrt(n / NODE_MAX) + 0.4) * MM, R_MISS)]
edges <- nodes[, {
  i <- seq_len(.N - 1)
  dx <- UMAP1[i + 1] - UMAP1[i]; dy <- UMAP2[i + 1] - UMAP2[i]; d <- sqrt(dx^2 + dy^2)
  s0 <- (r[i] + 0.3 * MM) / d; s1 <- 1 - (r[i + 1] + 0.5 * MM) / d
  .(x0 = UMAP1[i] + s0 * dx, y0 = UMAP2[i] + s0 * dy, x1 = UMAP1[i] + s1 * dx, y1 = UMAP2[i] + s1 * dy,
    dashed = !present[i + 1], keep = s1 > s0)
}, by = .(fate, genotype, gpanel)][keep == TRUE]   # nodes that touch get no arrow
node_layers <- function(nd, ed) list(
  # arrows with a white outline (user, 2026-11-24): a wider white copy underneath, then the dark arrow
  geom_segment(data = ed, aes(x0, y0, xend = x1, yend = y1), linewidth = 2.1, colour = "white", lineend = "round",
               arrow = arrow(length = unit(0.26, "cm"), type = "closed")),
  geom_segment(data = ed[dashed == FALSE], aes(x0, y0, xend = x1, yend = y1), linewidth = 0.9, colour = "grey15",
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")),
  geom_segment(data = ed[dashed == TRUE], aes(x0, y0, xend = x1, yend = y1), linewidth = 0.8, colour = "grey15",
               linetype = "22", arrow = arrow(length = unit(0.2, "cm"), type = "open")),
  ggforce::geom_circle(data = nd[present == FALSE], aes(x0 = UMAP1, y0 = UMAP2, r = R_MISS), inherit.aes = FALSE,
                       colour = "grey10", linetype = "22", linewidth = 1, fill = NA),
  geom_point(data = nd[present == TRUE], aes(UMAP1, UMAP2, fill = label, size = n), shape = 21, colour = "white",
             stroke = 1.1),
  scale_fill_manual(values = pal, guide = "none"),
  scale_size_area(max_size = NODE_SIZE_MAX, limits = c(0, NODE_MAX), guide = "none"))
p <- ggplot(fl, aes(UMAP1, UMAP2)) +
  geom_point(data = fbg, colour = GREY, size = POINT_SIZE, stroke = 0) +
  geom_point(aes(colour = p), size = POINT_SIZE, stroke = 0) + node_layers(nodes, edges) +
  facet_grid(fate ~ gpanel, switch = "y") + FATE_SCALE + umap_theme +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
  labs(title = "CellRank fate probabilities (pseudotime kernel, run on all cells)",
       subtitle = sprintf(paste0("nodes: the Slingshot path to that end state, per genotype (median position, area = cells)\n",
                                 "dashed circle = < %d cells of that genotype"), MIN_NODE))
save_fixed(p, file.path(PLOT, "cellrank_fates_umap.pdf"))
for (f in levels(fl$fate)) for (gt in names(GENO)) core_png(list(
  geom_point(data = fl[fate == f & genotype != gt], aes(UMAP1, UMAP2), colour = GREY, size = POINT_SIZE,
             stroke = 0),
  geom_point(data = fl[fate == f & genotype == gt], aes(UMAP1, UMAP2, colour = p), size = POINT_SIZE, stroke = 0),
  node_layers(nodes[fate == f & genotype == gt], edges[fate == f & genotype == gt]),
  FATE_SCALE), file.path(PNG, sprintf("cellrank_fate_%s_%s.png", gsub("[^A-Za-z0-9_]", "", sub(" .*", "", f)), gt)))

# 5. mean fate per subcluster x genotype
fs <- fread(file.path(TR, "cellrank_pseudotime_fate_summary.csv"))
fh <- melt(fs, id.vars = c("label", "genotype", "cells"), measure.vars = ends, variable.name = "fate", value.name = "p")
fh[, fate := factor(END_NAME[as.character(fate)], levels = END_NAME[ends])]
fh[, label := factor(label, levels = rev(ENDO_ORDER))]
fh[, gpanel := factor(GENO[genotype], levels = GENO)]
fh <- fh[cells >= 20]
p <- ggplot(fh, aes(fate, label, fill = p)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p), colour = p > 0.55), size = 2.1, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_gradientn(colours = c("#f7f7f7", "#fdbb84", "#e34a33", "#7f0000"), limits = c(0, 1), name = "mean fate\nprobability") +
  facet_wrap(~gpanel, nrow = 1) + coord_fixed() +
  labs(x = "end state", y = NULL, title = "Mean CellRank fate probability per subcluster",
       subtitle = "pseudotime kernel; groups with >= 20 cells") +
  theme_classic(base_size = 8) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                                       strip.background = element_blank(), strip.text = element_text(face = "bold"),
                                       axis.text.x = element_text(angle = 40, hjust = 1), plot.title = element_text(face = "bold"))
ggsave(file.path(PLOT, "cellrank_fate_heatmap.pdf"), p, width = 6.4, height = 3.6, useDingbats = FALSE)
cat("wrote", file.path(PLOT, "cellrank_fate_heatmap.pdf"), "\n")

# 6. data-driven end states per run
runs <- c(pseudotime = "all cells", pseudotime_ctr = "control only", pseudotime_KO = "KO only")
ms <- rbindlist(lapply(names(runs), function(r) {
  f <- file.path(TR, sprintf("cellrank_%s_macrostates.csv", r))
  if (!file.exists(f)) return(NULL)
  m <- fread(f); setnames(m, 1, "macrostate"); m[, run := runs[[r]]]
}), fill = TRUE)
ms[, run := factor(run, levels = runs)]
ms[, label := factor(main_label, levels = rev(ENDO_ORDER))]
ms <- ms[, .(terminal = any(terminal), n = .N), by = .(run, label)]
p <- ggplot(ms, aes(run, label)) +
  geom_point(aes(fill = label, shape = terminal), size = 4.2, colour = "black", stroke = 0.4) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), labels = c(`TRUE` = "terminal", `FALSE` = "macrostate, not terminal"),
                     name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "macrostate (main subcluster)", title = "End states CellRank finds without being told",
       subtitle = "GPCCA macrostates, pseudotime kernel, per run\n(circle = predicted terminal, triangle = not terminal)") +
  theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
                                       panel.grid.major = element_line(colour = "#eeeeee"))
ggsave(file.path(PLOT, "cellrank_end_states.pdf"), p, width = 4.4, height = 3.8, useDingbats = FALSE)
cat("wrote", file.path(PLOT, "cellrank_end_states.pdf"), "\n")

# 7. alpha vs beta bias (7.3i)
ab <- fread(file.path(TR, "alpha_bias", "D15", "alpha_bias_cells_pseudotime.csv")); setnames(ab, 1, "barcode")
ab <- merge(x[, .(barcode, panel)], ab[, .(barcode, lab, genotype, bias, S)], by = "barcode")
ab[, lab := factor(lab, levels = intersect(ENDO_ORDER, unique(lab)))]
setorder(ab, -S)
p <- ggplot(x, aes(UMAP1, UMAP2)) +
  geom_point(colour = GREY, size = POINT_SIZE, stroke = 0) +
  geom_point(data = merge(ab[S >= 0.05], x[, .(barcode, UMAP1, UMAP2)], by = "barcode"), aes(colour = bias),
             size = POINT_SIZE, stroke = 0) +
  geom_point(data = merge(ab[S < 0.05], x[, .(barcode, UMAP1, UMAP2)], by = "barcode"), colour = "black",
             shape = 4, size = 0.6, stroke = 0.3) +
  facet_wrap(~panel, nrow = 1) + BIAS_SCALE + umap_theme +
  labs(title = "Alpha vs beta bias of the progenitors (CellRank, pseudotime kernel)",
       subtitle = "ENP1, ENP2, ENP_FEV_high, pre_alpha_KO; x = P(alpha) + P(beta) < 0.05")
save_fixed(p, file.path(PLOT, "alpha_bias_umap.pdf"))
abu <- merge(ab, x[, .(barcode, UMAP1, UMAP2)], by = "barcode")
for (gt in names(GENO)) core_png(list(
  geom_point(data = x[genotype == gt], aes(UMAP1, UMAP2), colour = GREY, size = POINT_SIZE, stroke = 0),
  geom_point(data = abu[genotype == gt & S >= 0.05], aes(UMAP1, UMAP2, colour = bias), size = POINT_SIZE, stroke = 0),
  geom_point(data = abu[genotype == gt & S < 0.05], aes(UMAP1, UMAP2), colour = "black", shape = 4, size = 0.6,
             stroke = 0.3),
  BIAS_SCALE), file.path(PNG, sprintf("alpha_bias_umap_%s.png", gt)))
md <- ab[, .(m = median(bias, na.rm = TRUE), n = .N), by = .(lab, genotype)]
md[, txt := sprintf("%s median %.2f (n %s)", GENO[genotype], m, n_fmt(n))]
md[, genotype := factor(genotype, levels = names(GENO))]
md <- md[order(lab, genotype), .(txt = paste(txt, collapse = "\n")), by = lab]
strip <- setNames(sprintf("%s\n%s", md$lab, md$txt), md$lab)
ab[, labf := factor(strip[as.character(lab)], levels = strip[levels(lab)])]
p <- ggplot(ab, aes(bias, fill = genotype, colour = genotype)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.025, boundary = 0, alpha = 0.35, position = "identity",
                 linewidth = 0.25) +
  facet_wrap(~labf, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = GCOL, labels = GENO, name = NULL) + scale_colour_manual(values = GCOL, labels = GENO, name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "alpha bias = P(alpha) / (P(alpha) + P(beta))   (0 = beta, 1 = alpha)", y = "density",
       title = "Distribution of the alpha vs beta bias") +
  theme_classic(base_size = 8) + theme(strip.background = element_blank(), strip.text = element_text(size = 6.5),
                                       plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(PLOT, "alpha_bias_hist.pdf"), p, width = 9.4, height = 2.9, useDingbats = FALSE)
cat("wrote", file.path(PLOT, "alpha_bias_hist.pdf"), "\n")

# 8. gene dynamics along the beta, alpha and EC lineages
e <- fread(cmd = "gzcat reports/figure_data/D15_genes.csv.gz")
miss <- setdiff(GENES, names(e)); if (length(miss)) cat("not in D15_genes.csv.gz, skipped:", miss, "\n")
GENES <- intersect(GENES, names(e))
e <- e[match(sl$barcode, barcode)]; stopifnot(!anyNA(e$barcode))
L3 <- names(lend)[lend %in% c("pre_beta2", "pre_alpha", "pre_EC")]
cur <- rbindlist(lapply(L3, function(L) {
  pt <- sl[[paste0("pt_", L)]]; w <- sl[[paste0("w_", L)]]; w[is.na(w)] <- 0
  on <- w >= 0.5 & !is.na(pt)
  br <- seq(min(pt[on]), max(pt[on]), length.out = 11)
  b <- pmin(pmax(findInterval(pt, br, rightmost.closed = TRUE), 1), 10)
  rbindlist(lapply(GENES, function(g) rbindlist(lapply(names(GENO), function(gt) {
    r <- rbindlist(lapply(1:10, function(i) {
      m <- on & sl$genotype == gt & b == i
      data.table(bin = i, mid = (br[i] + br[i + 1]) / 2, n = sum(m),
                 mean = if (any(m)) weighted.mean(e[[g]][m], w[m]) else NA_real_)
    }))
    v <- ifelse(r$n >= 15, r$mean, NA)
    r[, smooth := sapply(seq_len(.N), function(i) { k <- max(1, i - 1):min(.N, i + 1); if (is.na(v[i])) NA else mean(v[k], na.rm = TRUE) })]
    r[, `:=`(gene = g, genotype = gt, lineage = END_NAME[lend[L]])]
  }))))
}))
cur[, gene := factor(gene, levels = GENES)]
cur[, lineage := factor(lineage, levels = c("beta", "alpha", "EC"))]
LCOL <- c(beta = unname(pal["SC_beta"]), alpha = unname(pal["SC_alpha"]), EC = unname(pal["SC_EC1"]))
p <- ggplot(cur, aes(mid, smooth, colour = lineage, linetype = genotype)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  facet_grid(gene ~ lineage, scales = "free", switch = "y") +
  scale_colour_manual(values = LCOL, guide = "none") +
  scale_linetype_manual(values = c(ctr = "solid", KO = "22"), breaks = names(GENO), labels = GENO, name = NULL) +
  labs(x = "Slingshot pseudotime", y = "mean log1p (curve-weighted, 3-bin smoothed)",
       title = "Gene dynamics along the beta, alpha and EC lineages",
       subtitle = "control solid, ISL1 KO dashed; no line where a genotype has < 15 cells in the bin") +
  theme_classic(base_size = 7) + theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, face = "bold"),
                                       strip.text.x = element_text(face = "bold"), strip.placement = "outside",
                                       plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(PLOT, "gene_dynamics.pdf"), p, width = 6.2, height = 0.95 * length(GENES) + 1.1, useDingbats = FALSE,
       limitsize = FALSE)
cat("wrote", file.path(PLOT, "gene_dynamics.pdf"), "\n")
