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

# ---- layout D (user, 2026-11-25: "I also want to try D, I need a plot that is narrow but tall"; then "8 pt, 80 wide,
# 150 tall, upright", no "from D18" labels and no legend; then "horizontal, ISL1 at top and others pointing down"):
# ISL1 -> cluster -> TF tree drawn top-down; then "try the 45 degree version": TF names below their dots at 45 degrees;
# C7 and C8 as separate nodes;
# node colour = within-cluster result (red enriched at any stage, blue depleted, grey n.s.); TF dot filled = TSS site at
# any stage, open = enhancer only (explained in the figure legend, not in the plot). Page 194 x 72 pt, TF names Arial
# regular 8 pt at 45 degrees, node labels Arial bold 5 pt; cairo_pdf with Arial embedded (Arial.ttf / Arial Bold.ttf copied
# from the Mac to ~/.local/share/fonts on Alpine, found through fontconfig).
# Output: <step 12 out dir>/isl1_target_TF_network.pdf
NODE_COL <- c(enriched = "#b2182b", depleted = "#2166ac", n.s. = "grey60")
nd <- rbindlist(lapply(names(PICK), function(col) {
  x <- d[column == col]; x[, .(n = .N, genes = list(gene)), by = cluster]
}))
nd[, state := sapply(cluster, function(cl) { w <- wi[group == cl & padj < 0.05]
  if (!nrow(w)) "n.s." else if (any(w$odds_ratio > 1)) "enriched" else "depleted" })]
W_PT <- 194; H_PT <- 72; PITCH <- 11.5; GAP_BLOCK <- 5.4; GAP_IN <- 3.6   # all in pt; 11.5 x sin 45 = 8.1 pt between names
x <- 14; lv <- list(); col_of <- d[, .(column = column[1]), by = cluster]
for (i in seq_len(nrow(nd))) {
  cl <- nd$cluster[i]; g <- nd$genes[[i]]
  if (i > 1) x <- x + PITCH + (if (col_of[cluster == cl, column] == col_of[cluster == nd$cluster[i - 1], column]) GAP_IN else GAP_BLOCK)
  xx <- x + (seq_along(g) - 1) * PITCH; x <- max(xx)
  lv[[i]] <- data.table(cluster = cl, gene = g, xt = xx, xc = mean(xx))
}
lv <- rbindlist(lv)
stopifnot(min(lv$xt) > 3, max(lv$xt) < W_PT - 3)
lv[tf, on = "gene", tss := i.tss_D15 | i.tss_D18 | i.tss_D22]
nodes <- unique(lv[, .(cluster, xc)])[nd, on = "cluster", state := i.state]
xI <- mean(range(lv$xt)); yI <- 63; yC <- 46.5; yT <- 36
PT <- ggplot2::.pt
pn <- ggplot() +
  geom_segment(data = nodes, aes(x = xI, y = yI, xend = xc, yend = yC), colour = "grey60", linewidth = 0.5 / PT * 1.07) +
  geom_segment(data = lv, aes(x = xc, y = yC, xend = xt, yend = yT), colour = "grey55", linewidth = 0.35 / PT * 1.07) +
  geom_point(aes(xI, yI), size = 17 / PT, colour = "grey15") +
  geom_text(aes(xI, yI, label = "ISL1"), colour = "white", fontface = "bold", size = 5 / PT, family = "Arial") +
  geom_point(data = nodes, aes(xc, yC, colour = state), size = 11 / PT) +
  geom_text(data = nodes, aes(xc, yC, label = cluster), colour = "white", fontface = "bold", size = 5 / PT,
            family = "Arial") +
  geom_point(data = lv, aes(xt, yT, fill = tss), shape = 21, colour = COL, size = 3.2 / PT, stroke = 0.4) +
  geom_text(data = lv, aes(xt + 1, yT - 3, label = gene), angle = 45, hjust = 1, vjust = 1, size = 8 / PT, family = "Arial") +
  scale_colour_manual(values = NODE_COL, guide = "none") +
  scale_fill_manual(values = c("TRUE" = COL, "FALSE" = "white"), guide = "none") +
  scale_x_continuous(limits = c(0, W_PT), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, H_PT), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
stopifnot("Arial" %in% systemfonts::system_fonts()$family)
cairo_pdf(file.path(od, "isl1_target_TF_network.pdf"), width = W_PT / 72, height = H_PT / 72, family = "Arial",
          bg = "white")
print(pn); par(ps = 8)
w <- strwidth(lv$gene, units = "inches", family = "Arial") * 72; invisible(dev.off())
r <- sqrt(0.5); left <- min(lv$xt + 1 - w * r); bottom <- min(yT - 3 - w * r - 6 * r)
cat(sprintf("wrote %s (%d x %d pt; names reach x = %.1f, y = %.1f pt)\n", file.path(od, "isl1_target_TF_network.pdf"),
            W_PT, H_PT, left, bottom))
stopifnot(left > 0, bottom > 0, max(lv$xt) < W_PT - 3)
