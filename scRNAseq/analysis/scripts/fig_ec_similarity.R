# Figures for step 7.1m: KO vs control SC_EC1 / SC_EC2 similarity (user, 2026-11-22). Runs locally (R, ggplot2).
# Inputs (copied from Alpine, git-ignored): reports/ec_similarity/ (7.1m outputs)
# Outputs (git-ignored) in plot/:
#   ec_corr_{day}.pdf      pseudobulk correlation (sample x label, gene-centred logCPM on the day's variable genes),
#                          endocrine labels, ordered by label then genotype
#   ec_dotplot_{day}.pdf   EC identity / EC1-EC2 / other-lineage genes, rows = label x genotype (>= 20 cells)
#   ec_de_counts.pdf       DE genes (edgeR QL, FDR < 0.05, |log2FC| >= 0.5): KO vs control within EC1 / EC2 next to
#                          the reference contrasts (EC2 vs EC1, EC1 vs beta) on the same 2 vs 2 design
#   ec_mixing.pdf          kNN mixing (1 = fully mixed, 0 = separate): genotype within EC1 / EC2 vs replicates and
#                          vs EC1 / EC2
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_ec_similarity.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
dir.create("plot", showWarnings = FALSE)
D <- "reports/ec_similarity"
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
GENO <- c(ctr = "Control", KO = "ISL1 KO")
ENDO <- c("SC_alpha", "SC_beta", "SC_delta", "SC_delta_GHRL", "GHRL+", "ENP", "SC_EC1", "SC_EC2")   # SC_delta_GHRL = D40 (2026-11-22)
th <- theme_classic(base_size = 8) + theme(plot.title = element_text(face = "bold"),
                                           strip.background = element_blank(), strip.text = element_text(face = "bold"))
days <- c("D22", "D40")

# 1. correlation heatmaps
for (d in days) {
  c <- fread(file.path(D, sprintf("corr_%s.csv", d))); m <- as.matrix(c[, -1]); rownames(m) <- c$pb
  info <- data.table(pb = rownames(m))   # "<label> <genotype> <rep>"; labels may contain spaces
  info[, `:=`(label = sub(" \\S+ \\S+$", "", pb), geno = sub("^.* (\\S+) \\S+$", "\\1", pb), rep = sub("^.* ", "", pb))]
  info <- info[label %in% ENDO][, label := factor(label, intersect(names(pal), ENDO))]
  setorder(info, label, geno, rep)
  info[, name := sprintf("%s %s %s", label, ifelse(geno == "KO", "KO", "ctr"), rep)]
  x <- as.data.table(as.table(m[info$pb, info$pb]))
  setnames(x, c("a", "b", "r"))
  x[, a := factor(info$name[match(a, info$pb)], info$name)][, b := factor(info$name[match(b, info$pb)], rev(info$name))]
  p <- ggplot(x, aes(a, b, fill = r)) + geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.7) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "Pearson r") +
    coord_fixed() + labs(x = NULL, y = NULL, title = sprintf("%s: pseudobulk profile correlation", d)) +
    th + theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.line = element_blank(),
               axis.ticks = element_blank())
  n <- nrow(info)
  ggsave(file.path("plot", sprintf("ec_corr_%s.pdf", d)), p, width = 1.6 + 0.24 * n, height = 1.2 + 0.24 * n,
         useDingbats = FALSE)
}

# 2. dot plots
GENES <- c("TPH1", "FEV", "LMX1A", "SLC18A1", "DDC", "TAC1", "CHGA", "PLCL1", "RUNX1T1", "AFF3", "NKX6-1", "PRKCH",
           "ISL1", "PAX6", "INS", "GCG", "SST", "ARX", "PDX1", "MAFB")
for (d in days) {
  x <- fread(file.path(D, sprintf("dotplot_%s.csv", d)))[label %in% ENDO & cells >= 20 & gene %in% GENES]
  x[, row := sprintf("%s %s", label, ifelse(geno == "KO", "KO", "ctr"))]
  ro <- unique(x[, .(row, label = factor(label, intersect(names(pal), ENDO)), geno)])[order(label, geno)]$row
  x[, row := factor(row, rev(ro))][, gene := factor(gene, intersect(GENES, gene))]
  x[, scaled := mean_log1p / max(mean_log1p), by = gene]
  p <- ggplot(x, aes(gene, row)) + geom_point(aes(size = pct, colour = scaled)) +
    scale_size_area(max_size = 3.5, name = "% expressing", breaks = c(25, 50, 75, 100)) +
    scale_colour_gradient(low = "grey90", high = "#b2182b", name = "mean log1p\n(scaled per gene)") +
    geom_vline(xintercept = c(7.5, 12.5), colour = "grey80", linewidth = 0.3) +
    labs(x = NULL, y = NULL, title = sprintf("%s: EC identity genes by genotype", d)) +
    th + theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic"))
  ggsave(file.path("plot", sprintf("ec_dotplot_%s.pdf", d)), p, width = 6.2, height = 1.6 + 0.2 * length(ro),
         useDingbats = FALSE)
}

# 3. DE counts
s <- fread(file.path(D, "de_summary.csv"))
TN <- c(KO_vs_ctr_SC_EC1 = "KO vs ctr, SC_EC1", KO_vs_ctr_SC_EC2 = "KO vs ctr, SC_EC2",
        ref_EC2_vs_EC1_ctr = "EC2 vs EC1 (ctr)", ref_EC2_vs_EC1_KO = "EC2 vs EC1 (KO)",
        ref_EC1_vs_beta_ctr = "EC1 vs beta (ctr)")
s[, test := factor(TN[test], rev(TN))][, kind := fifelse(grepl("^KO vs", test), "KO vs control", "reference")]
p <- ggplot(s, aes(DE, test, fill = kind)) + geom_col(width = 0.7) +
  geom_text(aes(label = format(DE, big.mark = ",")), hjust = -0.15, size = 2.3) +
  scale_fill_manual(values = c(`KO vs control` = "#b2182b", reference = "grey65"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) + facet_wrap(~day) +
  labs(x = "DE genes (FDR < 0.05, |log2FC| >= 0.5)", y = NULL, title = "Pseudobulk DE, 2 vs 2 samples") +
  th + theme(legend.position = "bottom")
ggsave("plot/ec_de_counts.pdf", p, width = 5.4, height = 2.6, useDingbats = FALSE)

# 4. mixing
mx <- fread(file.path(D, "mixing.csv"))
mx[, what := fifelse(grepl("^genotype", comparison), sprintf("KO vs ctr, %s", label),
               fifelse(grepl("^replicates", comparison),
                       sprintf("rep1 vs rep2 (%s), %s", sub(".*\\((\\w+) .*", "\\1", comparison), label),
                       sprintf("EC1 vs EC2 (%s)", sub(".*\\((\\w+)\\)", "\\1", comparison))))]
mx[, kind := fcase(grepl("^genotype", comparison), "KO vs control", grepl("^replicates", comparison),
                   "replicates (same cells)", default = "EC1 vs EC2 (different cells)")]
ord <- c(sprintf("KO vs ctr, SC_EC%d", 1:2), sprintf("rep1 vs rep2 (%s), SC_EC%d", rep(c("ctr", "KO"), 2), rep(1:2, each = 2)),
         sprintf("EC1 vs EC2 (%s)", c("ctr", "KO")))
mx[, what := factor(what, rev(ord))]
p <- ggplot(mx, aes(mixing_mean, what, fill = kind)) + geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = mixing_mean - mixing_sd, xmax = mixing_mean + mixing_sd), height = 0.3, linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = c(`KO vs control` = "#b2182b", `replicates (same cells)` = "grey65",
                               `EC1 vs EC2 (different cells)` = "grey30"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) + facet_wrap(~day) +
  labs(x = "kNN mixing (observed / expected share of the other group's neighbours)", y = NULL,
       title = "Neighbourhood mixing on SCT PCs (balanced subsamples, 20 seeds)") +
  th + theme(legend.position = "bottom")
ggsave("plot/ec_mixing.pdf", p, width = 6, height = 3, useDingbats = FALSE)
cat("wrote plot/ec_*.pdf\n")
