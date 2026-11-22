# Figures for step 7.1n: joint DESeq2 of the D22 + D40 SC_EC1 / SC_EC2 pseudobulks (user, 2026-11-22). Runs locally.
# Inputs (copied from Alpine, git-ignored): reports/ec_deseq2/ (7.1n outputs)
# Outputs (git-ignored) in plot/:
#   ec_deseq2_celltype_volcano.pdf  EC2 vs EC1, average of the 4 day x genotype groups (model A), core genes coloured,
#                                   EC markers / top TFs labelled
#   ec_deseq2_KO_scatter.pdf        KO vs control log2FC in EC1 (x) vs EC2 (y), average of the days (model B),
#                                   coloured by KO class; TFs of the shared class labelled
#   ec_deseq2_KO_upset.pdf          overlap of the KO calls: EC1 / EC2 x D22 / D40
#   ec_deseq2_heatmap_celltype.pdf  VST z-scores, 16 pseudobulks, top core EC1-high / EC2-high genes
#   ec_deseq2_heatmap_KO.pdf        VST z-scores, top shared KO up / down genes
#   ec_deseq2_GO.pdf                GO BP over-representation, top terms per gene set
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_ec_deseq2.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(pheatmap)
                                 library(UpSetR) })
D <- "reports/ec_deseq2"; dir.create("plot", showWarnings = FALSE)
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
th <- theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold"))
A <- fread(file.path(D, "results_A.csv")); B <- fread(file.path(D, "results_B.csv"))
gA <- fread(file.path(D, "genes_celltype.csv")); gB <- fread(file.path(D, "genes_KO.csv"))
smp <- fread(file.path(D, "samples.csv"))
MARK <- c("TAC1", "TRPA1", "CDX2", "TPH1", "LMX1A", "PCSK1", "NEUROD4", "ONECUT2", "CALCR", "NEUROG3", "SOX2", "SOX1",
          "ASCL1", "PAX6", "ARX", "SLC6A4", "PLCL1", "RUNX1T1", "VIM", "HES6", "LIN28B")

# 1. EC2 vs EC1 volcano
v <- A[!is.na(EC2_vs_EC1_average__padj), .(gene, lfc = EC2_vs_EC1_average__log2FC_shrunk,
                                           p = pmax(EC2_vs_EC1_average__padj, 1e-300), TF)]
v[, cls := fifelse(gene %in% gA[class == "EC2-high", gene], "EC2-high (core)",
          fifelse(gene %in% gA[class == "EC1-high", gene], "EC1-high (core)", "other"))]
lab <- v[gene %in% MARK | (TF & cls != "other" & rank(-abs(lfc)) <= 400)][order(-abs(lfc))][1:min(.N, 40)]
p <- ggplot(v, aes(lfc, -log10(p))) +
  geom_point(data = v[cls == "other"], colour = "grey85", size = 0.4) +
  geom_point(data = v[cls != "other"], aes(colour = cls), size = 0.5) +
  geom_text_repel(data = lab, aes(label = gene), size = 2.1, max.overlaps = 40, segment.size = 0.2,
                  min.segment.length = 0, fontface = "italic") +
  scale_colour_manual(values = c(`EC2-high (core)` = unname(pal["SC_EC2"]), `EC1-high (core)` = unname(pal["SC_EC1"])),
                      name = NULL) +
  labs(x = "log2FC EC2 vs EC1 (shrunk; average of D22/D40 x ctr/KO)", y = "-log10 padj",
       title = sprintf("SC_EC2 vs SC_EC1: %s EC2-high, %s EC1-high core genes",
                       format(gA[class == "EC2-high", .N], big.mark = ","), format(gA[class == "EC1-high", .N], big.mark = ","))) +
  th + theme(legend.position = "bottom")
ggsave("plot/ec_deseq2_celltype_volcano.pdf", p, width = 5, height = 4.6, useDingbats = FALSE)

# 2. KO scatter EC1 vs EC2
s <- B[, .(gene, x = KO_vs_ctr_EC1_average__log2FC_shrunk, y = KO_vs_ctr_EC2_average__log2FC_shrunk, TF)]
s[, KO_class := gB$KO_class[match(gene, gB$gene)]][is.na(KO_class), KO_class := "none"]   # classes in genes_KO.csv
s <- s[!is.na(x) & !is.na(y)]
KC <- c(shared = "#b2182b", `EC1 only` = unname(pal["SC_EC1"]), `EC2 only` = unname(pal["SC_EC2"]),
        `both, opposite direction` = "#e08214", none = "grey85")
s[, KO_class := factor(KO_class, names(KC))]
n <- s[, .N, by = KO_class]; nl <- setNames(sprintf("%s (%s)", n$KO_class, trimws(format(n$N, big.mark = ","))), n$KO_class)
LIM <- 7   # axis limit; the few genes beyond it are drawn at the edge (open triangles)
s[, out := abs(x) > LIM | abs(y) > LIM][, `:=`(x = pmax(pmin(x, LIM), -LIM), y = pmax(pmin(y, LIM), -LIM))]
lab <- s[KO_class == "shared" & (TF | gene %in% MARK)]
p <- ggplot(s, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_abline(slope = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(data = s[KO_class == "none"], colour = KC["none"], size = 0.35) +
  geom_point(data = s[KO_class != "none" & !out], aes(colour = KO_class), size = 0.6) +
  geom_point(data = s[KO_class != "none" & out], aes(colour = KO_class), size = 1.2, shape = 2, show.legend = FALSE) +
  geom_text_repel(data = s[out == TRUE], aes(label = gene), size = 2, fontface = "italic", colour = "grey30") +
  geom_text_repel(data = lab, aes(label = gene), size = 2, max.overlaps = 50, segment.size = 0.2, fontface = "italic",
                  min.segment.length = 0) +
  scale_colour_manual(values = KC, labels = nl, name = NULL) +
  coord_equal(xlim = c(-LIM, LIM), ylim = c(-LIM, LIM)) +
  labs(x = "KO vs control log2FC in SC_EC1 (shrunk, mean of D22/D40)", y = "KO vs control log2FC in SC_EC2",
       title = sprintf("ISL1 KO effect in EC1 vs EC2 (r = %.2f, all genes)", cor(s$x, s$y))) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 2)) +
  th + theme(legend.position = "bottom")
ggsave("plot/ec_deseq2_KO_scatter.pdf", p, width = 5.2, height = 5.8, useDingbats = FALSE)

# 3. UpSet of the KO calls per cell type x day
ups <- list(`EC1 D22` = B[KO_vs_ctr_EC1_D22__call == TRUE, gene], `EC1 D40` = B[KO_vs_ctr_EC1_D40__call == TRUE, gene],
            `EC2 D22` = B[KO_vs_ctr_EC2_D22__call == TRUE, gene], `EC2 D40` = B[KO_vs_ctr_EC2_D40__call == TRUE, gene])
pdf("plot/ec_deseq2_KO_upset.pdf", width = 5.5, height = 3.6, useDingbats = FALSE, onefile = FALSE)
print(upset(fromList(ups), sets = rev(names(ups)), keep.order = TRUE, order.by = "freq", nintersects = 15,
            main.bar.color = "#b2182b", sets.bar.color = "grey40", text.scale = 1.1,
            mainbar.y.label = "KO vs control DE genes", sets.x.label = "DE genes per test"))
invisible(dev.off())

# 4. heatmaps (VST z-scores)
vs <- fread(file.path(D, "vst.csv")); vm <- as.matrix(vs[, -1]); rownames(vm) <- vs$gene
smp[, `:=`(celltype = factor(celltype, c("EC1", "EC2")), genotype = factor(genotype, c("ctr", "KO")), day = factor(day))]
setorder(smp, celltype, genotype, day, rep); vm <- vm[, smp$id]
ann <- data.frame(celltype = smp$celltype, genotype = smp$genotype, day = smp$day, row.names = smp$id)
acol <- list(celltype = c(EC1 = unname(pal["SC_EC1"]), EC2 = unname(pal["SC_EC2"])),
             genotype = c(ctr = "grey80", KO = "grey20"), day = c(D22 = "#a6cee3", D40 = "#1f78b4"))
hm <- function(genes, rowlab, f, title) {
  z <- t(scale(t(vm[genes, ]))); z[z > 2.5] <- 2.5; z[z < -2.5] <- -2.5
  ra <- data.frame(set = rowlab, row.names = genes)
  rc <- setNames(c("#b2182b", "#2166ac", unname(pal["SC_EC2"]), unname(pal["SC_EC1"]))[seq_along(unique(rowlab))], unique(rowlab))
  pheatmap(z, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = ann, annotation_row = ra,
           annotation_colors = c(acol, list(set = rc)), show_colnames = FALSE, fontsize = 6.5, fontsize_row = 5.5,
           color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101), breaks = seq(-2.5, 2.5, length.out = 102),
           gaps_col = c(4, 8, 12), gaps_row = sum(rowlab == rowlab[1]), border_color = NA, main = title,
           filename = f, width = 4.6, height = 1.2 + 0.085 * length(genes))
}
top <- function(dt, n) dt[order(-abs(log2FC_shrunk))][1:min(.N, n), gene]
g2 <- top(gA[class == "EC2-high"], 30); g1 <- top(gA[class == "EC1-high"], 30)
hm(c(g2, g1), c(rep("EC2-high", length(g2)), rep("EC1-high", length(g1))), "plot/ec_deseq2_heatmap_celltype.pdf",
   "Core EC2 vs EC1 genes (VST z-score)")
sh <- gB[KO_class == "shared"][, m := (EC1_log2FC + EC2_log2FC) / 2]
gu <- sh[KO_direction == "up in KO"][order(-m)][1:min(.N, 30), gene]
gd <- sh[KO_direction == "down in KO"][order(m)][1:min(.N, 30), gene]
hm(c(gu, gd), c(rep("up in KO", length(gu)), rep("down in KO", length(gd))), "plot/ec_deseq2_heatmap_KO.pdf",
   "Shared KO genes, EC1 and EC2 (VST z-score)")

# 5. GO dot plot
go <- fread(file.path(D, "go_bp.csv"))
SETS <- c("EC1-high (core)", "EC2-high (core)", "KO down, shared EC1 + EC2", "KO up, shared EC1 + EC2",
          "KO down, EC1 only", "KO up, EC1 only", "KO down, EC2 only", "KO up, EC2 only")
tg <- go[padj < 0.05][order(pval), head(.SD, 6), by = set]
x <- go[term %in% tg$term & set %in% SETS]
x[, set := factor(set, SETS)][, term := factor(term, rev(unique(tg[order(match(set, SETS), pval), term])))]
p <- ggplot(x[padj < 0.05], aes(set, term)) +
  geom_point(aes(size = overlap, colour = -log10(padj))) +
  scale_colour_gradient(low = "#fdbb84", high = "#7f0000", name = "-log10 padj") +
  scale_size_area(max_size = 4, name = "genes") + scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = "GO biological process, top 6 terms per gene set") +
  th + theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(size = 6.5))
ggsave("plot/ec_deseq2_GO.pdf", p, width = 8, height = 1.4 + 0.13 * nlevels(x$term), useDingbats = FALSE)
cat("wrote plot/ec_deseq2_*.pdf\n")
