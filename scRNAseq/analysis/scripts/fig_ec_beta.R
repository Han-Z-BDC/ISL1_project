# Figures for step 7.1o: D22 KO EC vs control EC vs control beta (user, 2026-11-23). Runs locally (R, ggplot2).
# Inputs (copied from Alpine, git-ignored): reports/ec_beta_D22/ (7.1o outputs: quadrants.csv, beta_shift.csv,
#   pb_corr.csv, cells.csv.gz)
# Outputs (git-ignored) in plot/D22D40/ec_beta/:
#   ec_beta_quadrant.pdf     per EC type, every gene on the beta axis, selected identity genes (MARK) large and labelled (control beta vs control EC, FDR < 0.05, |log2FC| >= 1,
#                            INS left out): x = log2FC control beta vs control EC, y = log2FC KO EC vs control EC;
#                            coloured by quadrant of the KO call; slope (shift toward beta) and Spearman r; markers labelled
#   ec_beta_scores.pdf       per-cell beta score, EC score and beta index (beta - EC) for the 5 groups; dashed line =
#                            99th percentile of the control cells of the EC type (beyond = beta-leaning)
#   ec_beta_vs_INS.pdf       KO EC beta score by INS quintile (INS per 10,000 SoupX counts), control EC / beta medians
#   ec_beta_corr.pdf         correlation of the 5 group means (TMM logCPM, beta-axis genes, centred)
#   ec_beta_identity_heatmap.pdf  curated beta / EC identity genes: scaled group means and log2FC KO EC vs control EC
#   ec_beta_enrich_EC1.pdf   GO BP / MF / CC, KEGG, Reactome, WikiPathways, Hallmark of the EC1 quadrant categories
#                            (enrich_quadrants_EC1.csv, 7.1o_ec_beta_enrich.R)
#   ec_beta_umap_{beta,EC}_score.pdf   beta (blue) / EC (green) score of all D22 cells on the SCT B UMAP, control | KO (needs
#                            reports/labels/D22/sct_umaps.csv.gz and cells_all.csv.gz); high-resolution cores
#                            png/umap_{beta,EC}_score_{ctr,KO}.png (6 cm at 600 dpi, panel only)
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
# quadrant plot genes (user, 2026-11-23: show only these, coloured by their category)
MARK <- c("PDX1", "G6PC2", "IAPP", "PCSK1", "MAFA", "ERO1B", "TPH1", "SLC18A1", "FEV", "LMX1A", "DDC")
base <- theme_classic(base_size = 8) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
                                             plot.title = element_text(face = "bold"))

# 1. quadrant
q <- fread(file.path(D, "quadrants.csv")); sh <- fread(file.path(D, "beta_shift.csv"))
QCOL <- c("beta-high, up in KO EC" = "#d6604d", "beta-high, down in KO EC (check ambient)" = "#f4a582",
          "EC-high, down in KO EC" = "#2166ac", "EC-high, up in KO EC" = "#92c5de", "not changed in KO" = "#d9d9d9")
q[, quadrant := factor(quadrant, levels = names(QCOL))]
q[, panel := factor(sprintf("%s: shift toward beta %.2f, Spearman %.2f (%s axis genes)", ec, sh$slope[match(ec, sh$ec)],
                            sh$spearman[match(ec, sh$ec)], format(sh$axis_genes[match(ec, sh$ec)], big.mark = ",")))]
# all axis genes coloured by quadrant; the selected genes (MARK) larger, outlined and labelled (user, 2026-11-23)
qm <- q[gene %in% MARK]
miss <- lapply(c("EC1", "EC2"), function(k) setdiff(MARK, qm[ec == k, gene]))
cat("quadrant genes not on the axis: EC1", miss[[1]], "| EC2", miss[[2]], "\n")
setorder(q, -quadrant)
p <- ggplot(q, aes(beta_vs_EC, KO_vs_ctr_EC)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.35) +
  geom_point(aes(colour = quadrant), size = 0.7, stroke = 0) +
  geom_abline(data = q[, .(slope = sh$slope[match(ec[1], sh$ec)]), by = panel], aes(slope = slope, intercept = 0),
              colour = "black", linewidth = 0.5) +
  geom_point(data = qm, aes(fill = quadrant), shape = 21, size = 2.4, colour = "black", stroke = 0.4, show.legend = FALSE) +
  geom_text_repel(data = qm, aes(label = gene), size = 2.4, segment.size = 0.25, min.segment.length = 0,
                  box.padding = 0.4, max.overlaps = Inf, seed = 1) +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = QCOL, name = NULL, drop = FALSE) +
  scale_fill_manual(values = QCOL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.2), nrow = 2)) +
  labs(x = "log2FC control beta vs control EC", y = "log2FC ISL1 KO EC vs control EC",
       title = "D22: are KO EC shifted toward beta?",
       subtitle = paste0("beta-axis genes (control beta vs control EC, FDR < 0.05, |log2FC| >= 1, INS left out); colour = KO call ",
                         "(FDR < 0.05, |log2FC| >= 0.5); large outlined = selected identity genes\n",
                         "solid = fitted shift (0 = like control EC, 1 = like beta); dotted = y = x (fully beta-like)")) +
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

# 5. beta and EC scores on the D22 UMAP, all D22 cells, control | KO (user, 2026-11-23: "separately for KO and ctr";
#    "score all D22 cells"; gene expression stays red, scores in other colours: beta = blue, EC = green). Style as the
#    feature plots (fig_feature_dot.R): fixed square panels, the other genotype light grey underneath, low -> high from
#    grey; one colour scale per score for both panels, 1st to 99th percentile of all D22 cells (ends drawn at the limit)
suppressPackageStartupMessages(library(grid))
PANEL_CM <- 6; POINT_SIZE <- 0.6; OTHER <- "#f2f2f2"
SCORE_SCALE <- list(beta_score = c("#d9d9d9", "#9ecae1", "#3182bd", "#08306b"),   # grey -> blue
                    EC_score = c("#d9d9d9", "#a1d99b", "#31a354", "#00441b"))     # grey -> green
SCORE_NAME <- c(beta_score = "beta score", EC_score = "EC score")
SCORE_GENES <- c(beta_score = "100 beta-high genes", EC_score = "100 EC-high genes")
GENO <- c(ctr = "Control", KO = "ISL1 KO")
if (file.exists("reports/ec_beta_D22/cells_all.csv.gz")) {
um <- fread("reports/labels/D22/sct_umaps.csv.gz"); setnames(um, 1, "barcode")
ca <- fread(cmd = "gzcat reports/ec_beta_D22/cells_all.csv.gz"); setnames(ca, 1, "barcode")
u <- merge(ca, um[, .(barcode, UMAP1 = SCT_B_UMAP1, UMAP2 = SCT_B_UMAP2)], by = "barcode")
stopifnot(nrow(u) == nrow(ca))
ng <- u[, .N, by = genotype]; fl <- setNames(sprintf("%s (%s cells)", GENO[ng$genotype], format(ng$N, big.mark = ",")), ng$genotype)
u[, panel := factor(fl[genotype], levels = fl[names(GENO)])]
bgd <- u[, .(UMAP1, UMAP2, panel = factor(fl[setdiff(names(GENO), genotype)], levels = levels(u$panel))), by = genotype]
half <- max(diff(range(u$UMAP1)), diff(range(u$UMAP2))) / 2 * 1.08
cx <- mean(range(u$UMAP1)); cy <- mean(range(u$UMAP2))
for (sc_ in names(SCORE_SCALE)) {
  lim <- quantile(u[[sc_]], c(0.01, 0.99))
  d <- u[, .(UMAP1, UMAP2, panel, genotype, val = pmin(pmax(get(sc_), lim[1]), lim[2]))]
  setorder(d, val)
  p <- ggplot(d, aes(UMAP1, UMAP2)) +
    geom_point(data = bgd, colour = OTHER, size = POINT_SIZE, stroke = 0) +
    geom_point(aes(colour = val), size = POINT_SIZE, stroke = 0) +
    scale_colour_gradientn(colours = SCORE_SCALE[[sc_]], limits = lim, name = sub(" ", "\n", SCORE_NAME[[sc_]]),
                           breaks = function(l) pretty(l, 3)) +
    scale_x_continuous(limits = cx + c(-half, half), expand = c(0, 0)) +
    scale_y_continuous(limits = cy + c(-half, half), expand = c(0, 0)) +
    coord_fixed() + facet_wrap(~panel, nrow = 1) +
    labs(title = sprintf("D22: %s (%s; control beta vs control EC, INS left out)", SCORE_NAME[[sc_]], SCORE_GENES[[sc_]]),
         subtitle = "all D22 cells; scanpy score_genes on log1p") +
    theme_classic(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.title = element_text(face = "bold"),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8, hjust = 0),
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.25, "cm"), panel.spacing = unit(0.4, "cm"),
          plot.margin = margin(4, 4, 4, 4))
  g <- ggplotGrob(p); pp <- g$layout[grepl("^panel", g$layout$name), ]
  for (i in unique(pp$l)) g$widths[i] <- unit(PANEL_CM, "cm")
  for (i in unique(pp$t)) g$heights[i] <- unit(PANEL_CM, "cm")
  f <- file.path(PLOT, sprintf("ec_beta_umap_%s.pdf", sc_))
  pdf(f, width = convertWidth(sum(g$widths), "in", valueOnly = TRUE), height = convertHeight(sum(g$heights), "in", valueOnly = TRUE),
      useDingbats = FALSE); grid.draw(g); invisible(dev.off())
  # high-resolution cores (user, 2026-11-23; as the feature plots): plotting square only, one PNG per genotype, same
  # limits, scale and point size as the PDF, other genotype light grey underneath; the colour bar is in the PDF
  for (gt in names(GENO)) {
    pc <- ggplot(d[genotype == gt], aes(UMAP1, UMAP2)) +
      geom_point(data = d[genotype != gt], colour = OTHER, size = POINT_SIZE, stroke = 0) +
      geom_point(aes(colour = val), size = POINT_SIZE, stroke = 0) +
      scale_colour_gradientn(colours = SCORE_SCALE[[sc_]], limits = lim) +
      scale_x_continuous(limits = cx + c(-half, half), expand = c(0, 0)) +
      scale_y_continuous(limits = cy + c(-half, half), expand = c(0, 0)) +
      coord_fixed() + theme_void() + theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
    ggsave(file.path(PLOT, "png", sprintf("umap_%s_%s.png", sc_, gt)), pc, width = PANEL_CM, height = PANEL_CM,
           units = "cm", dpi = 600, bg = "white", create.dir = TRUE)
  }
}
} else cat("cells_all.csv.gz not there yet: score UMAPs skipped\n")

# 6. identity-gene heatmap (user, 2026-11-23): curated beta and EC identity genes; left = group means of the 7.1m
#    pseudobulks (TMM logCPM, >= 20 cells, scaled per gene across the 5 groups); right = log2FC KO EC vs control EC
#    from the 7.1o model, * = FDR < 0.05 and |log2FC| >= 0.5. INS kept here (ambient INS in control samples: its
#    control EC level is partly background)
suppressPackageStartupMessages(library(edgeR))
ID_GENES <- list(beta = c("INS", "ISL1", "PDX1", "NKX6-1", "MAFA", "MAFB", "NKX2-2", "NEUROD1", "PAX6", "IAPP", "G6PC2",
                          "HADH", "PCSK1", "ERO1B", "GCK", "ABCC8", "SLC30A8", "ADCYAP1"),
                 EC = c("TPH1", "SLC18A1", "DDC", "FEV", "LMX1A", "TAC1", "ADRA2A", "CHGA"))
cnt <- fread("reports/ec_similarity/pb_counts_D22.csv.gz"); gn <- cnt[[1]]; cm2 <- as.matrix(cnt[, -1]); rownames(cm2) <- gn
ps <- fread("reports/ec_similarity/pb_samples_D22.csv")[match(colnames(cm2), pb)]
ps[, grp := paste(genotype, sub("^SC_", "", label), sep = "_")]
k2 <- ps$cells >= 20 & ps$grp %in% GROUPS
lc <- cpm(calcNormFactors(DGEList(cm2[, k2])), log = TRUE, prior.count = 2)
gm <- sapply(GROUPS, function(g) rowMeans(lc[, ps$grp[k2] == g, drop = FALSE]))
hg <- data.table(gene = unlist(ID_GENES), set = rep(names(ID_GENES), lengths(ID_GENES)))
hg <- hg[gene %in% rownames(gm)]
z <- gm[hg$gene, ]; z <- (z - rowMeans(z)) / apply(z, 1, sd)
hz <- as.data.table(z, keep.rownames = "gene")
hz <- melt(hz, id.vars = "gene", variable.name = "group", value.name = "z")
hz[, group := factor(group, levels = GROUPS, labels = GLAB[GROUPS])]
hz[, gene := factor(gene, levels = rev(hg$gene))]
w <- fread(file.path(D, "de_wide.csv"))
fc <- rbindlist(lapply(c("EC1", "EC2"), function(k) w[gene %in% hg$gene, .(gene, col = sprintf("KO vs control %s", k),
  lfc = get(sprintf("lfc_KO_%s_vs_ctr_%s", k, k)), FDR = get(sprintf("FDR_KO_%s_vs_ctr_%s", k, k)))]))
fc[, star := fifelse(FDR < 0.05 & abs(lfc) >= 0.5, "*", "")]
fc[, gene := factor(gene, levels = rev(hg$gene))]
nb <- sum(hg$set == "beta")
th <- theme_classic(base_size = 8) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
  axis.text.x = element_text(angle = 40, hjust = 1), plot.title = element_text(face = "bold"), legend.key.width = unit(0.25, "cm"))
p1 <- ggplot(hz, aes(group, gene, fill = z)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = length(hg$gene) - nb + 0.5, linewidth = 0.5) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-2, 2), oob = scales::squish,
                       name = "scaled\nlogCPM") +
  labs(x = NULL, y = NULL, title = "Group means (pseudobulk)") + th
p2 <- ggplot(fc, aes(col, gene, fill = lfc)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%s", lfc, star)), size = 2) +
  geom_hline(yintercept = length(hg$gene) - nb + 0.5, linewidth = 0.5) +
  scale_fill_gradient2(low = "#4d9221", mid = "white", high = "#c51b7d", limits = c(-3, 3), oob = scales::squish,
                       name = "log2FC\nKO - ctr") +
  labs(x = NULL, y = NULL, title = "KO EC vs control EC") + th + theme(axis.text.y = element_blank())
suppressPackageStartupMessages(library(patchwork))
ggsave(file.path(PLOT, "ec_beta_identity_heatmap.pdf"), p1 + p2 + plot_layout(widths = c(5, 2)),
       width = 6.6, height = 0.15 * nrow(hg) + 1.6, useDingbats = FALSE)
cat("identity genes not in the tested set:", setdiff(unlist(ID_GENES), hg$gene), "\n")
# 7. enrichment of the EC1 quadrant categories (user, 2026-11-23; 7.1o_ec_beta_enrich.R): GO BP / MF / CC, KEGG Medicus,
#    KEGG Legacy, Reactome, WikiPathways and Hallmark (MSigDB) in one plot; rows = top 5 terms (padj < 0.05) per
#    category and collection, grouped by collection; columns = the four categories (x labels in the quadrant colours)
ef <- file.path(D, "enrich_quadrants_EC1.csv")
if (file.exists(ef)) {
  en <- fread(ef)
  CATS <- c("beta-high, up in KO EC", "beta-high, down in KO EC (check ambient)", "EC-high, down in KO EC", "EC-high, up in KO EC")
  COLLS <- c("GO BP", "GO MF", "GO CC", "KEGG Medicus", "KEGG Legacy", "Reactome", "WikiPathways", "Hallmark")
  n_set <- en[, .(n = set_size[1]), by = category]
  CLAB <- setNames(sprintf("%s\n(%s genes)", CATS, n_set$n[match(CATS, n_set$category)]), CATS)
  pretty_term <- function(t) {
    t <- sub("^(GOBP|GOMF|GOCC|REACTOME|HALLMARK|WP|KEGG|KEGG_MEDICUS)_", "", t)
    t <- tolower(gsub("_", " ", t)); t <- sub("^(.)", "\\U\\1", t, perl = TRUE)
    ifelse(nchar(t) > 62, paste0(substr(t, 1, 60), "..."), t)
  }
  top <- en[padj < 0.05][order(pval), head(.SD, 5), by = .(category, collection)]
  xe <- en[padj < 0.05 & paste(collection, term) %in% paste(top$collection, top$term)]
  xe[, label := pretty_term(term)]
  xe[, `:=`(category = factor(category, CATS), collection = factor(collection, intersect(COLLS, unique(collection))))]
  ord <- unique(top[order(match(collection, COLLS), match(category, CATS), pval), .(collection, term)])
  xe[, label := factor(paste(collection, term), levels = rev(paste(ord$collection, ord$term)),
                       labels = rev(pretty_term(ord$term)))]
  none <- setdiff(COLLS, unique(as.character(xe$collection)))
  p <- ggplot(xe, aes(category, label)) +
    geom_point(aes(size = overlap, colour = -log10(padj))) +
    facet_grid(collection ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_gradient(low = "#fdbb84", high = "#7f0000", name = "-log10 padj") +
    scale_size_area(max_size = 4, name = "genes") +
    scale_x_discrete(drop = FALSE, labels = CLAB) +
    labs(x = NULL, y = NULL, title = "D22 EC1: enrichment per beta-axis quadrant category",
         subtitle = paste0("MSigDB release 2026.1 (msigdbr); top 5 terms with padj < 0.05 per category and collection; ",
                           "universe = genes tested in the 7.1o model", if (length(none)) paste0("\nno term with padj < 0.05: ", paste(none, collapse = ", ")) else "")) +
    base + theme(plot.title.position = "plot", strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"),
                 axis.text.x = element_text(angle = 30, hjust = 1, colour = QCOL[CATS], face = "bold"),
                 axis.text.y = element_text(size = 6.3), panel.grid.major = element_line(colour = "#eeeeee"),
                 panel.spacing.y = unit(0.15, "cm"))
  ggsave(file.path(PLOT, "ec_beta_enrich_EC1.pdf"), p, width = 7.6, height = 2.2 + 0.13 * nlevels(xe$label), useDingbats = FALSE)
}
cat("wrote", PLOT, "\n")
