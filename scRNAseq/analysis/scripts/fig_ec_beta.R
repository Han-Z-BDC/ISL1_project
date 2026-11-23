# Figures for step 7.1o: D22 KO EC vs control EC vs control beta (user, 2026-11-23). Runs locally (R, ggplot2).
# Inputs (copied from Alpine, git-ignored): reports/ec_beta_D22/ (7.1o outputs: quadrants.csv, beta_shift.csv,
#   pb_corr.csv, cells.csv.gz)
# Outputs (git-ignored) in plot/D22D40/ec_beta/:
#   ec_beta_quadrant.pdf     per EC type, each beta-axis gene (control beta vs control EC, FDR < 0.05, |log2FC| >= 1,
#                            INS left out): x = log2FC control beta vs control EC, y = log2FC KO EC vs control EC;
#                            coloured by quadrant of the KO call; slope (shift toward beta) and Spearman r; markers labelled
#   ec_beta_scores.pdf       per-cell beta score, EC score and beta index (beta - EC) for the 5 groups; dashed line =
#                            99th percentile of the control cells of the EC type (beyond = beta-leaning)
#   ec_beta_vs_INS.pdf       KO EC beta score by INS quintile (INS per 10,000 SoupX counts), control EC / beta medians
#   ec_beta_corr.pdf         correlation of the 5 group means (TMM logCPM, beta-axis genes, centred)
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_ec_beta.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
PLOT <- file.path(Sys.getenv("PLOT_DIR", "plot/D22D40"), "ec_beta"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
D <- "reports/ec_beta_D22"
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
GROUPS <- c("ctr_EC1", "KO_EC1", "ctr_EC2", "KO_EC2", "ctr_beta")
GLAB <- c(ctr_EC1 = "Control EC1", KO_EC1 = "ISL1 KO EC1", ctr_EC2 = "Control EC2", KO_EC2 = "ISL1 KO EC2",
          ctr_beta = "Control beta")
GCOL <- c(ctr_EC1 = pal[["SC_EC1"]], KO_EC1 = pal[["SC_EC1"]], ctr_EC2 = pal[["SC_EC2"]], KO_EC2 = pal[["SC_EC2"]],
          ctr_beta = pal[["SC_beta"]])
GALPHA <- c(ctr_EC1 = 0.9, KO_EC1 = 0.35, ctr_EC2 = 0.9, KO_EC2 = 0.35, ctr_beta = 0.9)
MARK <- c("PDX1", "NKX6-1", "MAFA", "MAFB", "IAPP", "G6PC2", "HADH", "PCSK1", "ERO1B", "ADCYAP1", "NKX2-2", "ISL1",
          "TPH1", "SLC18A1", "DDC", "FEV", "LMX1A", "TAC1", "CHGA", "NEUROD1")
base <- theme_classic(base_size = 8) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
                                             plot.title = element_text(face = "bold"))

# 1. quadrant
q <- fread(file.path(D, "quadrants.csv")); sh <- fread(file.path(D, "beta_shift.csv"))
QCOL <- c("beta-high, up in KO EC" = "#d6604d", "beta-high, down in KO EC (check ambient)" = "#f4a582",
          "EC-high, down in KO EC" = "#2166ac", "EC-high, up in KO EC" = "#92c5de", "not changed in KO" = "#d9d9d9")
q[, quadrant := factor(quadrant, levels = names(QCOL))]
q[, panel := factor(sprintf("%s: shift toward beta %.2f, Spearman %.2f (%s axis genes)", ec, sh$slope[match(ec, sh$ec)],
                            sh$spearman[match(ec, sh$ec)], format(sh$axis_genes[match(ec, sh$ec)], big.mark = ",")))]
setorder(q, -quadrant)
p <- ggplot(q, aes(beta_vs_EC, KO_vs_ctr_EC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.35) +
  geom_point(aes(colour = quadrant), size = 0.7, stroke = 0) +
  geom_abline(data = q[, .(slope = sh$slope[match(ec[1], sh$ec)]), by = panel], aes(slope = slope, intercept = 0),
              colour = "black", linewidth = 0.5) +
  geom_text_repel(data = q[gene %in% MARK], aes(label = gene), size = 2.1, max.overlaps = 30, segment.size = 0.2,
                  min.segment.length = 0) +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = QCOL, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.2), nrow = 2)) +
  labs(x = "log2FC control beta vs control EC", y = "log2FC ISL1 KO EC vs control EC",
       title = "D22: are KO EC shifted toward beta?",
       subtitle = "beta-axis genes (control beta vs control EC, FDR < 0.05, |log2FC| >= 1, INS left out); solid = fitted shift (0 = like control EC, 1 = like beta), dotted = y = x") +
  base + theme(legend.position = "bottom")
ggsave(file.path(PLOT, "ec_beta_quadrant.pdf"), p, width = 9, height = 5, useDingbats = FALSE)

# 2. scores
x <- fread(cmd = "gzcat reports/ec_beta_D22/cells.csv.gz")
x[, group := factor(group, levels = GROUPS)]
sl <- melt(x, id.vars = c("group"), measure.vars = c("beta_score", "EC_score", "beta_index"), variable.name = "score")
sl[, score := factor(score, levels = c("beta_score", "EC_score", "beta_index"),
                     labels = c("beta score", "EC score", "beta index (beta - EC)"))]
thr <- x[group %in% c("ctr_EC1", "ctr_EC2"), .(q99 = quantile(beta_index, 0.99)), by = group]
thr <- rbind(thr[, .(group = as.character(group), q99)], thr[, .(group = sub("ctr", "KO", group), q99)])
thr[, group := factor(group, levels = GROUPS)]
thr[, score := factor("beta index (beta - EC)", levels = levels(sl$score))]
n <- x[, .N, by = group]
p <- ggplot(sl, aes(group, value)) +
  geom_violin(aes(fill = group, alpha = group), colour = NA, scale = "width") +
  geom_boxplot(width = 0.14, outlier.size = 0.25, outlier.stroke = 0, linewidth = 0.25, fill = "white") +
  geom_errorbar(data = thr, aes(x = group, ymin = q99, ymax = q99), inherit.aes = FALSE, width = 0.8,
                linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~score, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = GCOL, guide = "none") + scale_alpha_manual(values = GALPHA, guide = "none") +
  scale_x_discrete(labels = sprintf("%s\n(%s)", GLAB[levels(x$group)], format(n$N[match(levels(x$group), n$group)], big.mark = ","))) +
  labs(x = NULL, y = "score (scanpy score_genes, log1p)", title = "D22: beta and EC identity scores per cell",
       subtitle = "100 beta-high and 100 EC-high genes (control beta vs both control EC types, INS left out); dashed = 99th percentile of the control cells of that EC type") +
  base + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(PLOT, "ec_beta_scores.pdf"), p, width = 9, height = 3.8, useDingbats = FALSE)

# 3. beta score by INS quintile in KO EC
k <- x[group %in% c("KO_EC1", "KO_EC2")]
k[, ins_q := cut(frank(INS_per10k, ties.method = "first"), 5, labels = paste0("Q", 1:5)), by = group]
km <- k[, .(INS = median(INS_per10k)), by = .(group, ins_q)]
ref <- x[group %in% c("ctr_EC1", "ctr_EC2", "ctr_beta"), .(m = median(beta_score)), by = group]
ref <- rbind(ref[group != "ctr_beta", .(panel = sub("ctr", "KO", group), m, ref = "control EC median")],
             ref[group == "ctr_beta", .(panel = c("KO_EC1", "KO_EC2"), m, ref = "control beta median")])
k[, panel := group]; ref[, panel := factor(panel, levels = levels(k$group))]
p <- ggplot(k, aes(ins_q, beta_score)) +
  geom_violin(aes(fill = group), alpha = 0.35, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.25, outlier.stroke = 0, linewidth = 0.25, fill = "white") +
  geom_hline(data = ref, aes(yintercept = m, linetype = ref), linewidth = 0.45) +
  facet_wrap(~panel, nrow = 1, labeller = labeller(panel = GLAB)) +
  scale_fill_manual(values = GCOL, guide = "none") +
  scale_linetype_manual(values = c("control EC median" = "dashed", "control beta median" = "solid"), name = NULL) +
  geom_text(data = km[, .(panel = group, ins_q, INS)], aes(ins_q, -Inf, label = sprintf("%.0f", INS)), vjust = -0.4,
            size = 2.1, colour = "grey30") +
  labs(x = "INS quintile within the KO EC type (numbers = median INS per 10,000 SoupX counts)", y = "beta score",
       title = "D22 KO EC: does the beta score rise with INS?") +
  base + theme(legend.position = "top")
ggsave(file.path(PLOT, "ec_beta_vs_INS.pdf"), p, width = 6.4, height = 3.4, useDingbats = FALSE)
cat("INS medians per quintile:\n"); print(dcast(km, ins_q ~ group, value.var = "INS"))

# 4. group-mean correlation
cm <- fread(file.path(D, "pb_corr.csv")); m <- as.matrix(cm[, -1]); rownames(m) <- cm$group
cl <- as.data.table(as.table(m)); setnames(cl, c("a", "b", "r"))
cl[, `:=`(a = factor(a, levels = GROUPS, labels = GLAB[GROUPS]), b = factor(b, levels = rev(GROUPS), labels = GLAB[rev(GROUPS)]))]
p <- ggplot(cl, aes(a, b, fill = r)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "Pearson r") +
  coord_fixed() + labs(x = NULL, y = NULL, title = "D22 group means, beta-axis genes (centred)") +
  base + theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
ggsave(file.path(PLOT, "ec_beta_corr.pdf"), p, width = 4.4, height = 3.6, useDingbats = FALSE)
cat("wrote", PLOT, "\n")
