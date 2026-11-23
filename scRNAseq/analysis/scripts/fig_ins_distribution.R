# Figure: INS distribution in SC_EC1 / SC_EC2 by genotype, with control SC_beta as reference, D22 and D40 (user,
# 2026-11-23: choose the INS cut-off for "INS-high" KO EC cells). Runs locally (R, ggplot2).
# Input (copied from Alpine, git-ignored): reports/ec_similarity/hormone_check/cells.csv.gz (7.1m_hormone_check.py)
# INS is depth-normalised: INS per 10,000 SoupX counts of the cell (total_soupx), shown as log10(x + 1).
#   rows "SoupX" = INS_soupx, "raw" = INS_raw (before SoupX; same denominator), so the removed background is visible
#   lines: control SC_beta 10th (solid) and 25th (dashed) percentile of the SoupX value; numbers = % of the group's
#   cells at or above the 10th percentile
#   page 2: KO and control EC only, counts (not density), finer bins, to see whether INS-high cells form their own peak
# Output: plot/D22D40/ec_ins/INS_distribution.pdf
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_ins_distribution.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
PLOT <- file.path(Sys.getenv("PLOT_DIR", "plot/D22D40"), "ec_ins"); dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)
x <- fread(cmd = "gzcat reports/ec_similarity/hormone_check/cells.csv.gz")
x[, `:=`(ins_soupx = log10(INS_soupx / total_soupx * 1e4 + 1), ins_raw = log10(INS_raw / total_soupx * 1e4 + 1))]
GRP <- c("KO EC", "control EC", "control SC_beta")
GCOL <- c("KO EC" = "#d6604d", "control EC" = "#2166ac", "control SC_beta" = unname(pal["SC_beta"]))
ec <- x[label %in% c("SC_EC1", "SC_EC2")]
ec[, grp := fifelse(genotype == "KO", "KO EC", "control EC")]
be <- x[label == "SC_beta" & genotype == "ctr"]
# control beta shown in both EC columns as the reference
be2 <- rbindlist(lapply(c("SC_EC1", "SC_EC2"), function(l) copy(be)[, `:=`(col = l, grp = "control SC_beta")]))
ec[, col := label]
d <- rbind(ec[, .(day, col, grp, sample, ins_soupx, ins_raw)], be2[, .(day, col, grp, sample, ins_soupx, ins_raw)])
d <- melt(d, measure.vars = c("ins_soupx", "ins_raw"), variable.name = "counts", value.name = "v")
d[, counts := factor(fifelse(counts == "ins_soupx", "SoupX-corrected", "raw (before SoupX)"),
                     levels = c("SoupX-corrected", "raw (before SoupX)"))]
d[, grp := factor(grp, levels = GRP)]
d[, panel := factor(sprintf("%s %s", day, col))]
thr <- be[, .(p10 = quantile(ins_soupx, 0.10), p25 = quantile(ins_soupx, 0.25)), by = day]
thr2 <- merge(unique(d[, .(day, col, panel, counts)]), thr, by = "day")
pct <- merge(d[counts == "SoupX-corrected"], thr, by = "day")[, .(pct = 100 * mean(v >= p10), n = .N), by = .(panel, grp)]
pct[, txt := sprintf("%s: %.0f%% of %s", grp, pct, format(n, big.mark = ","))]
pct <- pct[order(panel, grp), .(txt = paste(txt, collapse = "\n")), by = panel]
pct[, counts := factor("SoupX-corrected", levels = levels(d$counts))]
base <- theme_classic(base_size = 8) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
                                             plot.title = element_text(face = "bold"), legend.position = "top")
p1 <- ggplot(d, aes(v, colour = grp, fill = grp)) +
  geom_density(alpha = 0.18, adjust = 0.6, linewidth = 0.45) +
  geom_vline(data = thr2, aes(xintercept = p10), linewidth = 0.4) +
  geom_vline(data = thr2, aes(xintercept = p25), linewidth = 0.4, linetype = "dashed") +
  geom_text(data = pct, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = -0.03, vjust = 1.15, size = 2) +
  facet_grid(counts ~ panel, scales = "free_y") +
  scale_colour_manual(values = GCOL, name = NULL) + scale_fill_manual(values = GCOL, name = NULL) +
  labs(x = "INS, log10(UMIs per 10,000 SoupX counts + 1)", y = "density",
       title = "INS in EC cells by genotype, control SC_beta as reference",
       subtitle = "lines: control SC_beta 10th (solid) / 25th (dashed) percentile of SoupX-corrected INS; % = cells at or above the 10th percentile") +
  base
p2 <- ggplot(d[counts == "SoupX-corrected" & grp != "control SC_beta"], aes(v, fill = grp)) +
  geom_histogram(binwidth = 0.05, boundary = 0, alpha = 0.75, position = "identity") +
  geom_vline(data = thr2[counts == "SoupX-corrected"], aes(xintercept = p10), linewidth = 0.4) +
  geom_vline(data = thr2[counts == "SoupX-corrected"], aes(xintercept = p25), linewidth = 0.4, linetype = "dashed") +
  facet_grid(grp ~ panel, scales = "free_y") +
  scale_fill_manual(values = GCOL, guide = "none") +
  labs(x = "INS, log10(UMIs per 10,000 SoupX counts + 1), SoupX-corrected", y = "cells",
       title = "INS in EC cells, counts per bin (is there a separate INS-high peak?)",
       subtitle = "lines: control SC_beta 10th (solid) / 25th (dashed) percentile") +
  base
f <- file.path(PLOT, "INS_distribution.pdf")
pdf(f, width = 9.5, height = 5.2, useDingbats = FALSE); print(p1); print(p2); invisible(dev.off())
cat("wrote", f, "\n")
print(thr)
