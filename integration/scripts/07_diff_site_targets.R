#!/usr/bin/env Rscript
# Step 07: target genes of the differentially bound ISL1 sites, and their KO effect.
#   Differential sites: step 27 DiffBind (background normalization), FDR < 0.05; gained = up in D18 vs D15 or
#   D22 vs D15, lost = down in either (step 03 site_stats.tsv; D22 vs D18 has no significant site).
#   Target genes: genes whose GREAT domain contains the site (step 01), limited to the genes tested in the bulk
#   analysis. Gene binding classes: gained (>= 1 gained site), lost (>= 1 lost site, no gained site), stable
#   (bound, no differential site), no site.
#   KO effect (step 02 gene_table.tsv): any day = LRT_genotype padj < 0.05; per day = KO vs Ctr padj < 0.05 (down =
#   log2FC < 0, lower in the KO: ISL1 activates; up = ISL1 represses); changes over time = LRT_interaction padj < 0.05.
#   Tests, gained against stable genes: logistic regression outcome ~ gained + log10(domain length) +
#   log2(1 + sites in the domain) (larger domains and more sites raise the chance of both a gained site and a KO
#   effect). Direction: among genes significant on a day, down ~ gained + the same covariates.
# usage: Rscript 07_diff_site_targets.R <integration results dir> <out_dir>
# Output in <out_dir>: diff_site_genes.tsv (one row per differential site and target gene), gene_classes.tsv,
#   summary.tsv (per binding class), tests.tsv

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
st[, change := fcase(
  gained & db_FDR_D18_vs_D15 < TH & db_log2FC_D18_vs_D15 > 0 & db_FDR_D22_vs_D15 < TH & db_log2FC_D22_vs_D15 > 0,
  "gained D18 and D22 vs D15",
  gained & db_FDR_D18_vs_D15 < TH & db_log2FC_D18_vs_D15 > 0, "gained D18 vs D15",
  gained, "gained D22 vs D15",
  lost, "lost", default = NA_character_)]
n_genes_per_site <- sg[, .(n_genes_linked = .N), by = site]

# ---- one row per differential site and target gene
ds <- merge(st[!is.na(change), .(site, change, db_log2FC_D18_vs_D15, db_FDR_D18_vs_D15, db_log2FC_D22_vs_D15,
                                 db_FDR_D22_vs_D15, db_log2FC_D22_vs_D18, db_FDR_D22_vs_D18)],
            sg[tested_in_bulk == 1, .(site, chr, summit, gene_id, distance_to_tss, in_basal)], by = "site")
ds <- merge(ds, n_genes_per_site, by = "site")
ds <- merge(ds, gt[, .(gene_id, gene_name, gene_type, group, padj_genotype, padj_interaction, ko_lfc_D15, ko_padj_D15,
                       ko_lfc_D18, ko_padj_D18, ko_lfc_D22, ko_padj_D22)], by = "gene_id")
setorder(ds, change, db_FDR_D18_vs_D15, site, gene_id)
setcolorder(ds, c("site", "chr", "summit", "change", "gene_id", "gene_name"))
fwrite(ds, file.path(od, "diff_site_genes.tsv"), sep = "\t", na = "NA", quote = FALSE)

# ---- gene binding classes
x <- merge(sg[tested_in_bulk == 1, .(site, gene_id)], st[, .(site, gained, lost)], by = "site")
gx <- x[, .(any_gained = any(gained), any_lost = any(lost)), by = gene_id]
g <- merge(gt, gx, by = "gene_id", all.x = TRUE)
g[is.na(any_gained), `:=`(any_gained = FALSE, any_lost = FALSE)]
g[, binding := fcase(any_gained, "gained", any_lost, "lost", bound, "stable", default = "no site")]
sig <- function(p) !is.na(p) & p < TH
g[, `:=`(ko_any = sig(padj_genotype), ko_changes = sig(padj_interaction))]
for (d in c("D15", "D18", "D22")) {
  g[, (paste0("sig_", d)) := sig(get(paste0("ko_padj_", d)))]
  g[, (paste0("down_", d)) := sig(get(paste0("ko_padj_", d))) & get(paste0("ko_lfc_", d)) < 0]
  g[, (paste0("up_", d)) := sig(get(paste0("ko_padj_", d))) & get(paste0("ko_lfc_", d)) > 0]
}
fwrite(g[, .(gene_id, gene_name, binding, n_sites, domain_kb, group, ko_any, ko_changes, ko_lfc_D15, ko_padj_D15,
             ko_lfc_D18, ko_padj_D18, ko_lfc_D22, ko_padj_D22)],
       file.path(od, "gene_classes.tsv"), sep = "\t", na = "NA", quote = FALSE)

lv <- c("gained", "stable", "no site", "lost")
summ <- g[, c(list(n_genes = .N, pct_ko_any = round(100 * mean(ko_any), 1),
                   pct_ko_changes = round(100 * mean(ko_changes), 1)),
              unlist(lapply(c("D15", "D18", "D22"), function(d) setNames(
                list(round(100 * mean(get(paste0("down_", d))), 1), round(100 * mean(get(paste0("up_", d))), 1)),
                c(paste0("pct_down_", d), paste0("pct_up_", d)))), recursive = FALSE)),
          by = binding][match(lv, binding)]
fwrite(summ, file.path(od, "summary.tsv"), sep = "\t", quote = FALSE)

# ---- tests: gained vs stable genes, adjusted for domain length and number of sites
b <- g[binding %in% c("gained", "stable")]
b[, gained := binding == "gained"]
cov <- "log10(domain_kb) + log2(1 + n_sites)"
one <- function(data, outcome, label) {
  fit <- glm(as.formula(paste(outcome, "~ gained +", cov)), data = data, family = binomial)
  co <- summary(fit)$coefficients["gainedTRUE", ]
  ft <- fisher.test(table(factor(data$gained, c(TRUE, FALSE)), factor(data[[outcome]], c(TRUE, FALSE))))
  data.table(test = label, n_gained = sum(data$gained), n_stable = sum(!data$gained),
             pct_gained = round(100 * mean(data[[outcome]][data$gained]), 1),
             pct_stable = round(100 * mean(data[[outcome]][!data$gained]), 1),
             odds_ratio = signif(unname(ft$estimate), 3), fisher_p = signif(ft$p.value, 3),
             adj_odds_ratio = signif(exp(co[1]), 3), adj_p = signif(co[4], 3))
}
tests <- rbind(
  one(b, "ko_any", "KO effect on any day"),
  one(b, "ko_changes", "KO effect changes over time"),
  rbindlist(lapply(c("D15", "D18", "D22"), function(d) one(b, paste0("sig_", d), paste("KO effect on", d)))),
  rbindlist(lapply(c("D15", "D18", "D22"), function(d)
    one(b[get(paste0("sig_", d))], paste0("down_", d), paste("Down in KO, among genes significant on", d)))))
tests[, adj_fdr := signif(p.adjust(adj_p, "BH"), 3)]
fwrite(tests, file.path(od, "tests.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("differential sites: %d gained, %d lost; site-gene rows: %d; genes: %d\n", st[gained == TRUE, .N],
            st[lost == TRUE, .N], nrow(ds), uniqueN(ds$gene_id)))
print(summ)
print(tests)
