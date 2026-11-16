#!/usr/bin/env Rscript
# Step 02: ISL1 direct targets among the genes the ISL1 KO changes in the bulk RNA-seq.
#   Differential   padj < 0.05 in any KO-effect table: KO vs Ctr at D15, D18 or D22 (DESeq2 Wald), or LRT_interaction
#                  (the KO effect changes over time). padj NA counts as not significant. The within-genotype time
#                  tables are not used (they do not involve the KO).
#   Direct target  differential and >= 1 ISL1 site (consensus peak on any day; 14,431 sites) in the gene's GREAT
#                  domain (step 01). When the peak appears and whether it changes are not considered.
#                  Check: >= 1 site whose nearest TSS is the gene's, within 50 kb.
#   Label          per day with padj < 0.05: activated = lower in KO (ISL1 activates it), repressed = higher in KO.
#                  Over the days: activated or repressed if every significant day agrees, mixed if they disagree.
#                  Significant only in LRT_interaction: the direction of the day with the largest |shrunk KO log2FC|,
#                  flagged interaction_only. Sign from the unshrunk DESeq2 log2FC (the one the padj belongs to).
#   KO-effect group (used by the later steps; mutually exclusive)
#     cluster1-6     LRT_interaction padj < 0.05, cluster from the bulk report
#     stays          LRT_genotype padj < 0.05, not in a cluster, shrunk KO log2FC with one sign at D15, D18, D22
#     other_KO       LRT_genotype padj < 0.05, not in a cluster, KO log2FC changes sign
#     no_KO_effect   LRT_genotype padj >= 0.05
# usage: Rscript 02_direct_targets.R <bulk DE_DESeq2 dir> <results/01_great> <sites.bed> <out_dir>
# Output in <out_dir>: gene_table.tsv (all bulk-tested genes with a GREAT domain), direct_targets.tsv (differential
#   and bound genes, with labels), summary.tsv (counts per label and source)

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
  g[ko, (paste0("ko_lfc_unshrunk_", d)) := i.log2FoldChange, on = "gene_id"]
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
n_no_domain_diff <- g[is.na(domain_kb) & (!is.na(padj_interaction) & padj_interaction < 0.05 |
                      (!is.na(ko_padj_D15) & ko_padj_D15 < 0.05) | (!is.na(ko_padj_D18) & ko_padj_D18 < 0.05) |
                      (!is.na(ko_padj_D22) & ko_padj_D22 < 0.05)), .N]
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

# ---- differential genes and their direction
sig_int <- !is.na(g$padj_interaction) & g$padj_interaction < 0.05
days <- c("D15", "D18", "D22")
sig <- function(p) !is.na(p) & p < 0.05
for (d in days) {
  p <- g[[paste0("ko_padj_", d)]]; l <- g[[paste0("ko_lfc_unshrunk_", d)]]
  g[, (paste0("label_", d)) := fifelse(sig(p) & l < 0, "activated", fifelse(sig(p) & l > 0, "repressed", ""))]
}
g[, sig_days := apply(.SD, 1, function(v) paste(days[v != ""], collapse = ",")), .SDcols = paste0("label_", days)]
g[, n_act := rowSums(.SD == "activated"), .SDcols = paste0("label_", days)]
g[, n_rep := rowSums(.SD == "repressed"), .SDcols = paste0("label_", days)]
g[, differential := sig_days != "" | sig_int]
g[, interaction_only := sig_int & sig_days == ""]
big <- apply(abs(as.matrix(g[, paste0("ko_lfc_", days), with = FALSE])), 1, which.max)
big_lfc <- as.matrix(g[, paste0("ko_lfc_", days), with = FALSE])[cbind(seq_len(nrow(g)), big)]
g[, label := fifelse(n_act > 0 & n_rep > 0, "mixed", fifelse(n_act > 0, "activated", fifelse(n_rep > 0, "repressed",
             fifelse(interaction_only, fifelse(big_lfc < 0, "activated", "repressed"), ""))))]
g[, direction_day := fifelse(interaction_only, days[big], "")]
g[, direct := differential & bound]
g[, direct_label := fifelse(direct, paste0("direct-", label), "")]
g[, c("n_act", "n_rep") := NULL]
stopifnot(g[differential == TRUE, all(label != "")], g[differential == FALSE, all(label == "")])

setorder(g, group, gene_id)
fwrite(g, file.path(od, "gene_table.tsv"), sep = "\t", na = "NA", quote = FALSE)
dt <- g[direct == TRUE, c("gene_id", "gene_name", "gene_type", "direct_label", "interaction_only", "direction_day",
                          "sig_days", paste0("label_", days), paste0("ko_lfc_", days), paste0("ko_padj_", days),
                          "padj_interaction", "group", "n_sites", "bound_near50", "n_sites_near50", "domain_kb"),
        with = FALSE]
dt[, lab_order := match(direct_label, paste0("direct-", c("activated", "repressed", "mixed")))]
setorder(dt, lab_order, gene_name, gene_id)
dt[, lab_order := NULL]
fwrite(dt, file.path(od, "direct_targets.tsv"), sep = "\t", na = "NA", quote = FALSE)

# ---- counts: per source (each day, interaction only, any) and label
src <- c(setNames(paste0("label_", days), paste0("KO vs Ctr ", days)))
sm <- rbindlist(c(
  lapply(names(src), function(sn) {
    x <- g[get(src[[sn]]) != ""]
    data.table(source = sn, label = c("activated", "repressed"),
               differential = c(x[get(src[[sn]]) == "activated", .N], x[get(src[[sn]]) == "repressed", .N]),
               direct = c(x[get(src[[sn]]) == "activated" & bound, .N], x[get(src[[sn]]) == "repressed" & bound, .N]))
  }),
  list(g[interaction_only == TRUE, .(source = "LRT_interaction only", differential = .N, direct = sum(bound)), by = label],
       g[differential == TRUE, .(source = "any (overall label)", differential = .N, direct = sum(bound)), by = label])),
  use.names = TRUE)
sm[, label := factor(label, levels = c("activated", "repressed", "mixed"))]
sm[, src_order := match(source, c(names(src), "LRT_interaction only", "any (overall label)"))]
setorder(sm, src_order, label)
sm[, src_order := NULL]
sm[, pct_direct := round(100 * direct / differential, 1)]
fwrite(sm, file.path(od, "summary.tsv"), sep = "\t", quote = FALSE)

message(sprintf("%d tested genes without a GREAT domain left out (%d of them differential)", n_no_domain, n_no_domain_diff))
print(sm)
print(dt[gene_name %in% c("INS", "GCG", "PAX6", "ARX", "PAX4", "MAFA", "MAFB", "PDX1", "NKX6-1", "ISL1", "NEUROG3"),
         .(gene_name, direct_label, sig_days, ko_lfc_D15, ko_lfc_D18, ko_lfc_D22)])
