# Step 12 (TFs and enrichment of the direct ISL1 targets per trajectory cluster; user, 2026-11-24: "focus on all the
# direct target genes in each cluster: 1. list all TFs in each gene cluster, write a csv table; 2. gene ontology
# analysis in each gene cluster, all three GO categories, KEGG, Reactome and Hallmark").
# Direct ISL1 target = a C1-C10 gene bound at any stage (genes.tsv bound_any: a CUT&RUN site at the TSS or ABC-linked,
#   D15 / D18 / D22); no KO expression criterion (bulk RNA-seq not used in the paper).
# 1. TFs: Lambert 2018 human TFs (scRNAseq study_config/annotation/TF_names_v_1.01.txt, as 7.3h) among the targets.
#    -> target_TFs.csv: cluster, gene, bound / TSS-bound per stage, first-bound stage, sites; target_TF_counts.csv.
# 2. Over-representation (fgsea::fora, hypergeometric) of each cluster's targets in MSigDB 2026.1 sets (the saved
#    scRNAseq table ${SC_OUT}/reference/msigdb_2026.1.Hs.csv.gz, as scRNAseq 7.1o): GO BP, GO MF, GO CC, KEGG Medicus,
#    KEGG Legacy, Reactome, Hallmark. Universe = all genes of genes.tsv (the 14,491 genes tested along the D15 lineages,
#    scRNAseq 7.3l), sets with 15-500 universe genes, BH padj per cluster x collection.
#    -> target_enrich.csv (all tested terms with overlap >= 1), isl1_target_enrich.pdf (one page per collection: top 5
#    terms per cluster with padj < 0.05, dot = overlap, colour = -log10 padj).
# Usage: Rscript 12_target_tf_enrich.R <step 12 out dir> <TF list> <msigdb csv.gz>
suppressPackageStartupMessages({ library(data.table); library(fgsea); library(ggplot2) })
a <- commandArgs(TRUE); od <- a[1]; tff <- a[2]; gsf <- a[3]
COLL <- c("GO BP", "GO MF", "GO CC", "KEGG Medicus", "KEGG Legacy", "Reactome", "Hallmark")
g <- fread(file.path(od, "genes.tsv"))
uni <- unique(g$gene)
tg <- g[grepl("^C[0-9]+$", group) & bound_any == TRUE]
CL <- paste0("C", 1:10)
tg[, group := factor(group, CL)]
cat("universe", length(uni), "genes; direct targets per cluster:\n"); print(tg[, .N, by = group][order(group)])

# ---- 1. TFs
tfs <- readLines(tff)
first <- function(b15, b18, b22) fifelse(b15, "D15", fifelse(b18, "D18", "D22"))
tf <- tg[gene %in% tfs, .(cluster = group, gene, bound_D15, bound_D18, bound_D22, tss_D15, tss_D18, tss_D22,
                          first_bound = first(bound_D15, bound_D18, bound_D22), n_sites)][order(cluster, gene)]
fwrite(tf, file.path(od, "target_TFs.csv"))
cnt <- merge(tg[, .(targets = .N), by = .(cluster = group)], tf[, .(TFs = .N, genes = paste(gene, collapse = ";")), by = cluster],
             by = "cluster", all.x = TRUE)[order(cluster)]
cnt[is.na(TFs), TFs := 0L]
fwrite(cnt, file.path(od, "target_TF_counts.csv"))
print(cnt[, .(cluster, targets, TFs)])

# ---- 2. enrichment
gs <- fread(cmd = sprintf("zcat '%s'", gsf))[collection %in% COLL & gene %in% uni]
print(gs[, .(sets = uniqueN(term)), by = collection])
res <- rbindlist(lapply(COLL, function(cn) {
  pw <- split(gs[collection == cn, gene], gs[collection == cn, term])
  pw <- pw[lengths(pw) >= 15 & lengths(pw) <= 500]
  rbindlist(lapply(CL, function(cl) {
    s <- tg[group == cl, gene]
    r <- as.data.table(fora(pathways = pw, genes = s, universe = uni, minSize = 15, maxSize = 500))
    if (!nrow(r)) return(NULL)
    r[, `:=`(cluster = cl, collection = cn, targets = length(s), terms_tested = length(pw),
             genes = sapply(overlapGenes, paste, collapse = ","))][, overlapGenes := NULL]
  }))
}))
setnames(res, c("pathway", "size"), c("term", "term_size"))
res[, fold := (overlap / targets) / (term_size / length(uni))]
setcolorder(res, c("cluster", "collection", "term", "overlap", "targets", "term_size", "fold", "pval", "padj",
                   "terms_tested", "genes"))
res[, cluster := factor(cluster, CL)]
res <- res[order(cluster, collection, pval)]
fwrite(res, file.path(od, "target_enrich.csv"))
cat("\nterms with padj < 0.05 (cluster x collection):\n")
print(dcast(res[padj < 0.05, .N, by = .(cluster, collection)], cluster ~ collection, value.var = "N", fill = 0))

# ---- figure: one page per collection, top 5 terms per cluster (padj < 0.05)
short <- function(x) { x <- sub("^(GOBP|GOMF|GOCC|KEGG_MEDICUS|KEGG|REACTOME|HALLMARK)_", "", x)
  x <- gsub("_", " ", tolower(x)); ifelse(nchar(x) > 60, paste0(substr(x, 1, 58), "..."), x) }
pdf(file.path(od, "isl1_target_enrich.pdf"), width = 8.5, height = 9, useDingbats = FALSE)
for (cn in COLL) {
  sig <- res[collection == cn & padj < 0.05]
  top <- sig[order(padj), head(.SD, 5), by = cluster]
  if (!nrow(top)) {
    plot.new(); title(sprintf("%s: no term with padj < 0.05 in any cluster", cn)); next
  }
  terms <- unique(top$term)
  d <- res[collection == cn & term %in% terms]
  ord <- top[order(cluster, padj), unique(term)]
  d[, lab := factor(short(term), rev(unique(short(ord))))]
  d[, sig := padj < 0.05]
  p <- ggplot(d[overlap > 0], aes(cluster, lab)) +
    geom_point(aes(size = overlap, colour = pmin(-log10(padj), 10), alpha = sig)) +
    scale_colour_gradient(low = "#fcbba1", high = "#99000d", name = "-log10 padj\n(capped 10)") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), name = "padj < 0.05") +
    scale_size_area(max_size = 5, name = "target genes\nin term") + scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = sprintf("%s: direct ISL1 targets per D15 trajectory cluster", cn),
         subtitle = sprintf("top 5 terms per cluster (padj < 0.05); universe %s genes tested along the D15 lineages",
                            format(length(uni), big.mark = ","))) +
    theme_bw(base_size = 8) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  print(p)
}
invisible(dev.off())
cat("wrote", file.path(od, c("target_TFs.csv", "target_TF_counts.csv", "target_enrich.csv", "isl1_target_enrich.pdf")),
    sep = "\n")
