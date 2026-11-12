#!/usr/bin/env Rscript
# Step 02: ISL1 direct targets per bulk KO-effect group.
#   Groups (mutually exclusive; padj NA counts as not significant)
#     cluster1-6     LRT_interaction padj < 0.05, cluster from the bulk report
#     stays          LRT_genotype padj < 0.05, not in a cluster, shrunk KO log2FC with one sign at D15, D18, D22
#     other_KO       LRT_genotype padj < 0.05, not in a cluster, KO log2FC changes sign
#     no_KO_effect   LRT_genotype padj >= 0.05 (reference group)
#   Direct target  >= 1 ISL1 site (any day) in the gene's GREAT domain (step 01)
#   Check          >= 1 site whose nearest TSS (any gene in the GTF) is this gene's, within 50 kb
#   Tests per group vs no_KO_effect: Fisher's exact test, and logistic regression bound ~ group +
#   log10(domain length) (longer domains collect more sites by chance); BH over groups.
# usage: Rscript 02_direct_targets.R <bulk DE_DESeq2 dir> <results/01_great> <sites.bed> <out_dir>
# Output in <out_dir>: gene_table.tsv, enrichment.tsv

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
bulk <- args[1]; great <- args[2]; sites_bed <- args[3]; od <- args[4]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) fread(f, sep = "\t", na.strings = "NA", quote = "")

# ---- KO-effect groups
g <- rd(file.path(bulk, "LRT_genotype.tsv"))[, .(gene_id, gene_name, gene_type, padj_genotype = padj)]
g[rd(file.path(bulk, "LRT_interaction.tsv")), padj_interaction := i.padj, on = "gene_id"]
for (d in c("D15", "D18", "D22")) {
  ko <- rd(file.path(bulk, sprintf("KO_vs_Ctr_%s.tsv", d)))
  g[ko, (paste0("ko_lfc_", d)) := i.log2FC_shrunk, on = "gene_id"]
  g[ko, (paste0("ko_padj_", d)) := i.padj, on = "gene_id"]
}
g[rd(file.path(bulk, "timecourse_clusters_LRT_interaction.tsv")), cluster := i.cluster, on = "gene_id"]
sig_int <- !is.na(g$padj_interaction) & g$padj_interaction < 0.05
stopifnot(identical(sig_int, !is.na(g$cluster)))
s <- sign(as.matrix(g[, .(ko_lfc_D15, ko_lfc_D18, ko_lfc_D22)]))
one_sign <- s[, 1] != 0 & s[, 1] == s[, 2] & s[, 1] == s[, 3]
sig_gen <- g$padj_genotype < 0.05
g[, group := fifelse(sig_int, paste0("cluster", cluster),
             fifelse(sig_gen & one_sign, "stays", fifelse(sig_gen, "other_KO", "no_KO_effect")))]
lv <- c("no_KO_effect", paste0("cluster", 1:6), "stays", "other_KO")
g[, group := factor(group, levels = lv)]

# ---- binding: GREAT domains (step 01)
reg <- fread(file.path(great, "gene_regions.bed"), header = FALSE,
             col.names = c("chr", "start", "end", "gene_id", "gene_name", "strand", "tss", "basal_start",
                           "basal_end", "gene_type", "tested"))
g[reg, `:=`(domain_kb = (i.end - i.start) / 1000), on = "gene_id"]
n_no_domain <- g[is.na(domain_kb), .N]
g <- g[!is.na(domain_kb)]
sg <- rd(file.path(great, "site_gene.tsv"))
g[, n_sites := 0L][sg[, .N, by = gene_id], n_sites := i.N, on = "gene_id"]
g[, bound := n_sites > 0]

# ---- check: nearest TSS within 50 kb
sites <- fread(sites_bed, header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "site"))
sites[, pos := start + (end - start) %/% 2]
tss <- reg[, .(chr, pos = tss, tss_pos = tss, gene_id)]
setkey(tss, chr, pos)
near <- tss[sites[, .(chr, pos, site)], on = .(chr, pos), roll = "nearest"]
near <- near[!is.na(gene_id) & abs(pos - tss_pos) <= 50000]
g[, n_sites_near50 := 0L][near[, .N, by = gene_id], n_sites_near50 := i.N, on = "gene_id"]
g[, bound_near50 := n_sites_near50 > 0]

setorder(g, group, gene_id)
fwrite(g, file.path(od, "gene_table.tsv"), sep = "\t", na = "NA", quote = FALSE)

# ---- enrichment per group
ref <- g[group == "no_KO_effect"]
res <- rbindlist(lapply(c(great = "bound", nearest_tss_50kb = "bound_near50"), function(y) {
  fit <- glm(reformulate(c("group", "log10(domain_kb)"), y), data = g, family = binomial)
  co <- summary(fit)$coefficients
  rbindlist(lapply(lv, function(l) {
    sub <- g[group == l]
    r <- data.table(group = l, n_genes = nrow(sub), n_bound = sum(sub[[y]]),
                    pct_bound = round(100 * mean(sub[[y]]), 1),
                    median_domain_kb = round(median(sub$domain_kb), 1))
    if (l != "no_KO_effect") {
      ft <- fisher.test(matrix(c(sum(sub[[y]]), sum(!sub[[y]]), sum(ref[[y]]), sum(!ref[[y]])), 2))
      cn <- paste0("group", l)
      r[, `:=`(odds_ratio = signif(unname(ft$estimate), 3), fisher_p = signif(ft$p.value, 3),
               adj_odds_ratio = signif(exp(co[cn, 1]), 3), adj_p = signif(co[cn, 4], 3))]
    }
    r
  }), fill = TRUE)
}), idcol = "assignment", fill = TRUE)
res[, adj_fdr := signif(p.adjust(adj_p, "BH"), 3), by = assignment]
fwrite(res, file.path(od, "enrichment.tsv"), sep = "\t", na = "NA", quote = FALSE)

message(sprintf("%d tested genes without a GREAT domain left out", n_no_domain))
print(res)
