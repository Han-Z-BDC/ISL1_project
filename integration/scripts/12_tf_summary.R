# Step 12: candidate TFs among the direct ISL1 targets of C1, C4 and C7 + C8, with ISL1 binding per stage (user,
# 2026-11-25: "show candidate TFs per cluster, focus on C1, C4 and C8 ... just a summary of TF names"; layout B = names +
# D15 / D18 / D22 binding dots). TFs chosen by the user from the ranked table (bound at D15 or at >= 2 stages first,
# then best-lineage FDR of the 7.3l pseudobulk KO test): C1 top 5; C4 PAX6, ETV1, ARX, HHEX, CREB3L2; C7 and C8 shown
# together: FEV, ASCL2 (C7), LMX1A, MNX1, HEY2 (C8). ISL1 itself left out (user: many of its binding sites come from the
# ISL1 cDNA).
# Dots: filled = ISL1 site within +/- 1 kb of the TSS that day; open = site linked only through an ABC enhancer; grey =
#   no site that day. Header: direct targets (C1-C10 genes bound at any stage), TFs among them (ISL1 not counted),
#   within-cluster result (within_clusters.tsv: cluster vs the other nine, Fisher two-sided, BH FDR < 0.05).
# Inputs: <step 12 out dir>/target_TFs.csv, target_TF_counts.csv, within_clusters.tsv
# Output: <step 12 out dir>/isl1_target_TF_summary.pdf
# Usage: Rscript 12_tf_summary.R <step 12 out dir>
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
od <- commandArgs(TRUE)[1]
PICK <- list("C1" = c("SIM1", "PDX1", "PBX3", "FOXO1", "TBX3"),
             "C4" = c("PAX6", "ETV1", "ARX", "HHEX", "CREB3L2"),
             "C7 + C8" = c("FEV", "ASCL2", "LMX1A", "MNX1", "HEY2"))
DAYS <- c("D15", "D18", "D22"); COL <- "#08306b"; GREY <- "grey35"
tf <- fread(file.path(od, "target_TFs.csv"))
cn <- fread(file.path(od, "target_TF_counts.csv"))
wi <- fread(file.path(od, "within_clusters.tsv"))
cn[, tfs := TFs - grepl("(^|;)ISL1(;|$)", genes)]
res_txt <- function(cl) {
  w <- wi[group == cl & padj < 0.05]
  if (!nrow(w)) return("n.s.")
  paste(unlist(lapply(c("enriched", "depleted"), function(k) {
    s <- w[(odds_ratio > 1) == (k == "enriched"), stage]; if (length(s)) paste(k, paste(s, collapse = ", ")) }))
    , collapse = "; ")
}
head_txt <- sapply(names(PICK), function(col) {
  cls <- strsplit(col, " + ", fixed = TRUE)[[1]]
  r <- if (length(cls) == 1) res_txt(cls) else paste(sprintf("%s %s", cls, sapply(cls, res_txt)), collapse = "; ")
  sprintf("%d targets, %d TFs\n%s", sum(cn[match(cls, cluster), targets]), sum(cn[match(cls, cluster), tfs]), r)
})
d <- rbindlist(lapply(names(PICK), function(col) {
  x <- tf[match(PICK[[col]], gene)]; stopifnot(!anyNA(x$cluster))
  x[, `:=`(column = col, y = -(seq_len(.N)))]
}))
multi <- sapply(strsplit(names(PICK), " + ", fixed = TRUE), length) > 1
d[, lab := fifelse(column %in% names(PICK)[multi], sprintf('italic("%s")~scriptstyle("  (%s)")', gene, cluster),
                   sprintf('italic("%s")', gene))]
d[column %in% names(PICK)[multi] & cluster == names(which.max(table(cluster[column %in% names(PICK)[multi]]))),
  lab := sprintf('italic("%s")', gene)]                          # tag only the minority cluster of a merged column
dots <- rbindlist(lapply(seq_along(DAYS), function(k) {
  dd <- DAYS[k]
  d[, .(column, y, x = 2.55 + 0.42 * (k - 1),
        state = fifelse(get(paste0("tss_", dd)), "TSS site", fifelse(get(paste0("bound_", dd)), "enhancer only (ABC)",
                                                                       "not bound")))]
}))
dots[, state := factor(state, levels = c("TSS site", "enhancer only (ABC)", "not bound"))]
hd <- data.table(column = names(PICK), h = head_txt)
dl <- CJ(column = names(PICK), k = seq_along(DAYS))[, `:=`(x = 2.55 + 0.42 * (k - 1), lab = DAYS[k])]
for (z in list(d, dots, hd, dl)) z[, column := factor(column, levels = names(PICK))]
p <- ggplot() +
  geom_text(data = hd, aes(0, 2.25, label = column), hjust = 0, vjust = 1, fontface = "bold", size = 2.5) +
  geom_text(data = hd, aes(0, 1.45, label = h), hjust = 0, vjust = 1, size = 1.55, colour = GREY, lineheight = 0.95) +
  geom_text(data = dl, aes(x, 0.05, label = lab), size = 1.45, colour = GREY) +
  geom_text(data = d, aes(0, y, label = lab), parse = TRUE, hjust = 0, size = 2.2) +
  geom_point(data = dots, aes(x, y, fill = state, size = state, colour = state), shape = 21, stroke = 0.35) +
  facet_grid(. ~ column) +
  scale_fill_manual(values = c("TSS site" = COL, "enhancer only (ABC)" = "white", "not bound" = "grey80"), name = NULL) +
  scale_colour_manual(values = c("TSS site" = COL, "enhancer only (ABC)" = COL, "not bound" = "grey80"), name = NULL) +
  scale_size_manual(values = c("TSS site" = 1.6, "enhancer only (ABC)" = 1.6, "not bound" = 0.6), name = NULL) +
  scale_x_continuous(limits = c(0, 3.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-5.5, 2.35), expand = c(0, 0)) +
  theme_void(base_size = 6) +
  theme(strip.text = element_blank(), legend.position = "bottom", legend.text = element_text(size = 5),
        legend.key.size = unit(0.25, "cm"), panel.spacing.x = unit(0.15, "cm"), plot.margin = margin(4, 4, 2, 4))
ggsave(file.path(od, "isl1_target_TF_summary.pdf"), p, width = 3.4, height = 1.45, useDingbats = FALSE, bg = "white")
fwrite(d[, .(column, cluster, gene, bound_D15, bound_D18, bound_D22, tss_D15, tss_D18, tss_D22, first_bound, n_sites)],
       file.path(od, "isl1_target_TF_summary.csv"))
print(hd)
cat("wrote", file.path(od, "isl1_target_TF_summary.pdf"), "\n")
