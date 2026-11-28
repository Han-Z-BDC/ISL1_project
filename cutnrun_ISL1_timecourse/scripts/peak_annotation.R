# Step 14: one genomic class per peak and the stacked bar plot (see 14_peak_annotation.sb).
# Inputs: <homer annotation dir>/{D15,D18,D22,genome_random}.annotatePeaks.txt (step 10c),
#   <out dir>/{set}.repeats.tsv (window centre in RepeatMasker; peak, name, family, class).
# Output in <out dir>: peaks.tsv (set, peak, HOMER annotation, repeat name / family / class, class),
#   summary.tsv (set, class, n, pct), peak_annotation.pdf (peak counts per day | % per day and random genome windows).
#   Rscript peak_annotation.R <homer annotation dir> <out dir>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(TRUE); A <- a[1]; O <- a[2]
SETS <- c(D15 = "D15", D18 = "D18", D22 = "D22", genome_random = "Random\ngenome")
CLS <- c("Promoter-TSS", "TTS", "Exon", "Repeat", "Intron", "Intergenic")
COL <- setNames(c("#08306b", "#4292c6", "#9ecae1", "#d6604d", "#fdae6b", "grey75"), CLS)

p <- rbindlist(lapply(names(SETS), function(s) {
  h <- fread(file.path(A, paste0(s, ".annotatePeaks.txt")), sep = "\t", quote = "", select = c(1, 8))
  setnames(h, c("peak", "homer"))
  r <- fread(file.path(O, paste0(s, ".repeats.tsv")), header = FALSE, col.names = c("peak", "rep_name", "rep_family", "rep_class"))
  stopifnot(!anyDuplicated(h$peak), all(r$peak %in% h$peak))
  h <- merge(h, r, by = "peak", all.x = TRUE)
  b <- sub(" \\(.*", "", h$homer)
  h[, class := fcase(b == "promoter-TSS", "Promoter-TSS", b == "TTS", "TTS",
                     b %in% c("exon", "5' UTR", "3' UTR", "non-coding"), "Exon",
                     !is.na(rep_class), "Repeat", b == "intron", "Intron", b == "Intergenic", "Intergenic",
                     default = NA_character_)]
  if (anyNA(h$class)) stop("unclassified HOMER annotation in ", s, ": ", paste(unique(b[is.na(h$class)]), collapse = ", "))
  cbind(set = s, h)
}))
fwrite(p[, .(set, peak, homer, rep_name, rep_family, rep_class, class)], file.path(O, "peaks.tsv"), sep = "\t")

s <- p[, .N, by = .(set, class)][CJ(set = names(SETS), class = CLS), on = .(set, class)]
s[is.na(N), N := 0L]; s[, pct := 100 * N / sum(N), by = set]
setorder(s[, `:=`(o1 = match(set, names(SETS)), o2 = match(class, CLS))], o1, o2)
fwrite(s[, .(set, class, n = N, pct = round(pct, 2))], file.path(O, "summary.tsv"), sep = "\t")
print(dcast(s, class ~ set, value.var = "N")[match(CLS, class)])

s[, `:=`(setl = factor(SETS[set], levels = SETS), class = factor(class, levels = rev(CLS)))]
th <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.text = element_text(colour = "black"), legend.key.size = unit(0.3, "cm"), plot.title = element_text(size = 7))
g1 <- ggplot(s[set != "genome_random"], aes(setl, N, fill = class)) + geom_col(width = 0.7, colour = NA) +
  scale_fill_manual(values = COL, breaks = CLS, name = NULL) +
  scale_y_continuous(expand = expansion(c(0, 0.03)), labels = scales::comma) +
  labs(x = NULL, y = "Consensus peaks", title = "Peaks") + th
g2 <- ggplot(s, aes(setl, pct, fill = class)) + geom_col(width = 0.7, colour = NA) +
  scale_fill_manual(values = COL, breaks = CLS, name = NULL) +
  scale_y_continuous(expand = expansion(c(0, 0)), limits = c(0, 100.01)) +
  labs(x = NULL, y = "% of peaks / windows", title = "Fraction") + th
f <- file.path(O, "peak_annotation.pdf")
cairo_pdf(f, width = 4.2, height = 2.1, family = "Arial")
grid::grid.newpage(); lay <- grid::grid.layout(1, 2, widths = grid::unit(c(1.75, 2.45), "in"))
grid::pushViewport(grid::viewport(layout = lay))
print(g1 + theme(legend.position = "none"), vp = grid::viewport(layout.pos.col = 1))
print(g2, vp = grid::viewport(layout.pos.col = 2))
invisible(dev.off()); cat("wrote", f, "\n")
