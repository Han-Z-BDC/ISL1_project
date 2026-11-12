#!/usr/bin/env Rscript
# Step 04: candidate time-dependent direct ISL1 targets: a site with differential binding between days whose
# timing matches the gene's KO-effect pattern.
#   Matching rules (KO effect grows when binding grows; sign of the KO effect is not used):
#     cluster1       KO effect strongest at D18           -> a site gained in D18 vs D15
#     cluster3, 4    KO effect emerging from D18, growing -> a site gained in D18 vs D15 or D22 vs D15
#     cluster5       KO effect late (D22)                 -> a site gained in D22 vs D15
#     cluster2, 6    KO effect early, fading by D22       -> a site lost in D18 vs D15 or D22 vs D15
#   Gained / lost: step 27 DiffBind, background normalization, FDR < 0.05 (step 03 site_stats.tsv).
#   Sites link to genes through GREAT domains (step 01).
#   Also reported per candidate: whether the gene's own KO effect is significant (padj < 0.05) on the day of the
#   matching binding change (D18 for a D18 vs D15 change, D22 for a D22 vs D15 change).
#   Gene-level check: share of genes with >= 1 gained site per KO-effect group against genes without a KO effect,
#   logistic regression with log10(domain length) and log2(1 + sites in the domain) (BH over groups).
# usage: Rscript 04_candidates.R <integration results dir> <out_dir>
# Output in <out_dir>: candidates.tsv (one row per candidate gene and matching site), candidate_genes.tsv
#   (one row per candidate gene), summary.tsv (per group), gained_enrichment.tsv

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
res <- args[1]; od <- args[2]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
TH <- 0.05
rd <- function(f) fread(f, sep = "\t", quote = "", na.strings = "NA")

sg <- rd(file.path(res, "01_great", "site_gene.tsv"))
gt <- rd(file.path(res, "02_direct_targets", "gene_table.tsv"))
st <- rd(file.path(res, "03_binding_change", "site_stats.tsv"))
n_genes_per_site <- sg[, .(n_genes_linked = .N), by = site]

# per site-gene link: which matching changes the site has
x <- merge(sg[tested_in_bulk == 1, .(site, chr, summit, gene_id, distance_to_tss, in_basal)],
           st[, .(site, db_log2FC_D18_vs_D15, db_FDR_D18_vs_D15, db_log2FC_D22_vs_D15, db_FDR_D22_vs_D15)], by = "site")
x <- merge(x, n_genes_per_site, by = "site")
x <- merge(x, gt[, .(gene_id, gene_name, group, ko_lfc_D15, ko_padj_D15, ko_lfc_D18, ko_padj_D18, ko_lfc_D22,
                     ko_padj_D22, domain_kb, n_sites)], by = "gene_id")
x[, `:=`(up18 = db_FDR_D18_vs_D15 < TH & db_log2FC_D18_vs_D15 > 0, up22 = db_FDR_D22_vs_D15 < TH & db_log2FC_D22_vs_D15 > 0,
         dn18 = db_FDR_D18_vs_D15 < TH & db_log2FC_D18_vs_D15 < 0, dn22 = db_FDR_D22_vs_D15 < TH & db_log2FC_D22_vs_D15 < 0)]
sig <- function(p) !is.na(p) & p < TH
x[, change := fcase(
  group == "cluster1" & up18, "gained D18 vs D15",
  group %in% c("cluster3", "cluster4") & up18 & up22, "gained D18 and D22 vs D15",
  group %in% c("cluster3", "cluster4") & up18, "gained D18 vs D15",
  group %in% c("cluster3", "cluster4") & up22, "gained D22 vs D15",
  group == "cluster5" & up22, "gained D22 vs D15",
  group %in% c("cluster2", "cluster6") & (dn18 | dn22), "lost",
  default = NA_character_)]
x[, ko_sig_on_change_day := fcase(
  change == "gained D18 vs D15", sig(ko_padj_D18),
  change == "gained D22 vs D15", sig(ko_padj_D22),
  change == "gained D18 and D22 vs D15", sig(ko_padj_D18) | sig(ko_padj_D22),
  change == "lost", fifelse(dn18, sig(ko_padj_D18), FALSE) | fifelse(dn22, sig(ko_padj_D22), FALSE),
  default = NA)]

cand <- x[!is.na(change)]
cand[, best_FDR := pmin(fifelse(up18 | dn18, db_FDR_D18_vs_D15, 1), fifelse(up22 | dn22, db_FDR_D22_vs_D15, 1))]
setorder(cand, group, best_FDR, gene_id, site)
cols <- c("group", "gene_id", "gene_name", "site", "chr", "summit", "distance_to_tss", "in_basal", "n_genes_linked",
          "change", "db_log2FC_D18_vs_D15", "db_FDR_D18_vs_D15", "db_log2FC_D22_vs_D15", "db_FDR_D22_vs_D15",
          "ko_lfc_D15", "ko_padj_D15", "ko_lfc_D18", "ko_padj_D18", "ko_lfc_D22", "ko_padj_D22", "ko_sig_on_change_day")
fwrite(cand[, ..cols], file.path(od, "candidates.tsv"), sep = "\t", na = "NA", quote = FALSE)

genes <- cand[, .(n_matching_sites = .N, best_site = site[1], best_change = change[1], best_FDR = min(best_FDR),
                  max_abs_log2FC = round(max(pmax(abs(db_log2FC_D18_vs_D15 * (up18 | dn18)),
                                                  abs(db_log2FC_D22_vs_D15 * (up22 | dn22)))), 3),
                  any_in_basal = any(in_basal == 1), ko_sig_on_change_day = any(ko_sig_on_change_day),
                  ko_lfc_D15 = ko_lfc_D15[1], ko_lfc_D18 = ko_lfc_D18[1], ko_lfc_D22 = ko_lfc_D22[1],
                  n_sites_in_domain = n_sites[1]),
              by = .(group, gene_id, gene_name)]
setorder(genes, group, -ko_sig_on_change_day, best_FDR, gene_id)
fwrite(genes, file.path(od, "candidate_genes.tsv"), sep = "\t", na = "NA", quote = FALSE)

# ---- per-group summary
gx <- x[, .(any_gained = any(up18 | up22), any_lost = any(dn18 | dn22)), by = gene_id]
g <- merge(gt[, .(gene_id, group, bound, domain_kb, n_sites)], gx, by = "gene_id", all.x = TRUE)
g[is.na(any_gained), `:=`(any_gained = FALSE, any_lost = FALSE)]
lv <- c("no_KO_effect", paste0("cluster", 1:6), "stays", "other_KO")
g[, group := factor(group, levels = lv)]
summ <- g[, .(n_genes = .N, n_bound = sum(bound), n_with_gained_site = sum(any_gained), n_with_lost_site = sum(any_lost)),
          by = group]
summ <- merge(summ, genes[, .(n_candidates = .N, n_candidates_ko_sig_on_day = sum(ko_sig_on_change_day)), by = group],
              by = "group", all.x = TRUE)
summ[is.na(n_candidates), `:=`(n_candidates = 0L, n_candidates_ko_sig_on_day = 0L)]
summ[, pct_with_gained_site := round(100 * n_with_gained_site / n_genes, 1)]
setorder(summ, group)
fwrite(summ, file.path(od, "summary.tsv"), sep = "\t", quote = FALSE)

# ---- gene-level check: genes with >= 1 gained site, per group vs no KO effect
fit <- glm(any_gained ~ group + log10(domain_kb) + log2(1 + n_sites), data = g, family = binomial)
co <- summary(fit)$coefficients
ref <- g[group == "no_KO_effect"]
enr <- rbindlist(lapply(lv[-1], function(l) {
  sub <- g[group == l]
  ft <- fisher.test(matrix(c(sum(sub$any_gained), sum(!sub$any_gained), sum(ref$any_gained), sum(!ref$any_gained)), 2))
  cn <- paste0("group", l)
  data.table(group = l, n_genes = nrow(sub), pct_with_gained_site = round(100 * mean(sub$any_gained), 1),
             pct_reference = round(100 * mean(ref$any_gained), 1), odds_ratio = signif(unname(ft$estimate), 3),
             fisher_p = signif(ft$p.value, 3), adj_odds_ratio = signif(exp(co[cn, 1]), 3), adj_p = signif(co[cn, 4], 3))
}))
enr[, adj_fdr := signif(p.adjust(adj_p, "BH"), 3)]
fwrite(enr, file.path(od, "gained_enrichment.tsv"), sep = "\t", quote = FALSE)

print(summ)
print(enr)
print(genes[, .N, by = .(group, best_change)])
