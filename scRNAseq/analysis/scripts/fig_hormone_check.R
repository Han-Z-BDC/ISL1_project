# Figure: per-cell hormone reads in SC_EC1 / SC_EC2 by genotype, D22 and D40 (user, 2026-11-22: INS, then GCG).
# Input (copied from Alpine, git-ignored): reports/ec_similarity/hormone_check/cells.csv.gz (7.1m_hormone_check.py)
# Output: plot/{gene}_ec_distribution.pdf, per gene given
#   top: UMIs per cell (SoupX-corrected, log10(x + 1)), violin + box per label x genotype, the hormone's own control
#        cell type as the reference (INS: SC_beta, GCG: SC_alpha); dashed line = its 10th percentile ("hormone-cell level")
#   bottom: the same for raw counts (before SoupX), showing the background SoupX removed
#   numbers: % cells with the gene > 0 / % at hormone-cell level
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_hormone_check.R INS GCG
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
dir.create("plot", showWarnings = FALSE)
REF <- c(INS = "SC_beta", GCG = "SC_alpha")
genes <- commandArgs(TRUE); stopifnot(length(genes) > 0, all(genes %in% names(REF)))
x0 <- fread(cmd = "gzcat reports/ec_similarity/hormone_check/cells.csv.gz")
GC <- c(ctr = "grey60", KO = "#b2182b")
for (g in genes) {
  ref <- REF[[g]]
  x <- x0[label %in% c("SC_EC1", "SC_EC2") | (label == ref & genotype == "ctr")]
  x[, `:=`(sx = get(paste0(g, "_soupx")), raw = get(paste0(g, "_raw")))]
  x[, grp := factor(fifelse(label == ref, sprintf("%s ctr", ref), sprintf("%s %s", label, genotype)),
                    c("SC_EC1 ctr", "SC_EC1 KO", "SC_EC2 ctr", "SC_EC2 KO", sprintf("%s ctr", ref)))]
  thr <- x[label == ref, .(thr = quantile(sx, 0.1)), by = day]
  pct <- x[label != ref, .(lab = sprintf("%d%%\n%d%%", round(100 * mean(sx > 0)),
                                         round(100 * mean(sx >= thr[day == .BY$day, thr])))), by = .(day, grp)]
  long <- rbind(x[, .(day, grp, genotype, value = sx, what = "SoupX-corrected counts")],
                x[, .(day, grp, genotype, value = raw, what = "raw counts (before SoupX)")])
  long[, what := factor(what, c("SoupX-corrected counts", "raw counts (before SoupX)"))]
  pct[, what := factor("SoupX-corrected counts", levels(long$what))]
  top <- ceiling(max(log10(long$value + 1)) * 10) / 10
  p <- ggplot(long, aes(grp, log10(value + 1), fill = genotype)) +
    geom_violin(scale = "width", linewidth = 0.2, colour = "grey30") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
    geom_hline(data = thr, aes(yintercept = log10(thr + 1)), linetype = "dashed", linewidth = 0.3, colour = "grey30") +
    geom_text(data = pct, aes(grp, top + 0.35, label = lab), inherit.aes = FALSE, size = 1.9, lineheight = 0.85) +
    scale_fill_manual(values = GC, labels = c(ctr = "Control", KO = "ISL1 KO"), name = NULL) +
    scale_y_continuous(limits = c(0, top + 0.6)) +
    facet_grid(what ~ day) +
    labs(x = NULL, y = sprintf("%s UMIs per cell, log10(x + 1)", g), title = sprintf("%s in EC cells by genotype", g),
         subtitle = sprintf("numbers: %% cells with %s > 0 / %% at %s level (>= 10th percentile of control %s, dashed)",
                            g, sub("SC_", "", ref), ref)) +
    theme_classic(base_size = 8) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 6.5),
          strip.background = element_blank(), strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")
  f <- file.path("plot", sprintf("%s_ec_distribution.pdf", g))
  ggsave(f, p, width = 6.4, height = 5, useDingbats = FALSE)
  cat("wrote", f, "\n")
}
