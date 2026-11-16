#!/usr/bin/env Rscript
# Step 09: target genes of the overlap classes of the per-day ISL1 peak calls, against the bulk KO-effect groups.
#   Classes    from CUT&RUN step 12 (peak_overlap/sites_days.tsv): each of the 14,431 sites carries the days whose
#              consensus peaks (macs3_ge2) contain it; the 7 combinations of days are the overlap classes. The peak-call
#              analysis itself (library call rates, presence against the DiffBind counts) is CUT&RUN step 12.
#   Step 08    per class: batch-corrected VST per day (median) and the share of sites in each step 08 LRT /
#              binding-pattern group.
#   Genes      linked genes tested in the bulk RNA-seq (step 02), two links as in step 08: GREAT domains (main) and every
#              bulk-tested gene with its TSS within 50 kb of the site centre (check). A gene is in a class if at least one
#              of its linked sites is (a gene can be in several classes); "only" genes have all their sites in the class.
#              Per class and KO-effect group: genes in the class against all other linked genes; Fisher test, and
#              logistic regression in_class ~ in_KO_group + log2(1 + n_sites) (+ log10(domain_kb) for GREAT), since
#              genes with more sites are more often in any class and more often KO-affected; BH per link.
# usage: Rscript 09_peak_overlap.R <CUT&RUN peak_overlap dir> <integration results dir> <markers.tsv> <out_dir>
# Output in <out_dir>: class_groups.tsv (per class: step 08 VST and groups), gene_classes.tsv (per link and gene: KO
#   group, marker, sites per class), gene_class_links.tsv (per link, class and KO group)

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
ovl <- args[1]; res <- args[2]; mfile <- args[3]; od <- args[4]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
CLASSES <- c("D15,D18,D22", "D18,D22", "D15,D18", "D15,D22", "D18", "D22", "D15")

si <- fread(file.path(ovl, "sites_days.tsv"))[, .(site, chr, start, end, days_called)]
stopifnot(all(si$days_called %in% CLASSES))
si[, class := factor(days_called, levels = CLASSES)]

# ---- per class: step 08 profile and group
sc <- fread(file.path(res, "08_site_clusters", "site_clusters.tsv"))[
  , .(site, varying, vst_D15, vst_D18, vst_D22, group = k4)]
n0 <- nrow(si)
sv <- merge(si, sc, by = "site")
stopifnot(nrow(sv) == n0, !anyNA(sv$group))
cg <- sv[, .(sites = .N, vst_D15 = round(median(vst_D15), 3), vst_D18 = round(median(vst_D18), 3),
             vst_D22 = round(median(vst_D22), 3), pct_varying = round(100 * mean(varying), 1),
             pct_group1 = round(100 * mean(group == 1), 1), pct_group2 = round(100 * mean(group == 2), 1),
             pct_group3 = round(100 * mean(group == 3), 1), pct_group4 = round(100 * mean(group == 4), 1)),
         by = class][order(class)]
fwrite(cg, file.path(od, "class_groups.tsv"), sep = "\t", quote = FALSE)
print(cg)

# ---- gene level: overlap classes of each gene's sites against its KO-effect group
sgl <- fread(file.path(res, "01_great", "site_gene.tsv"), sep = "\t", quote = "")
gt <- fread(file.path(res, "02_direct_targets", "gene_table.tsv"), sep = "\t", quote = "", na.strings = "NA")
mk <- fread(mfile, skip = "cell_type")
reg <- fread(file.path(res, "01_great", "gene_regions.bed"), header = FALSE,
             col.names = c("chr", "start", "end", "gene_id", "gene_name", "strand", "tss", "basal_start",
                           "basal_end", "gene_type", "tested"))
sb <- si[, .(chr, s0 = start + (end - start) %/% 2L, s1 = start + (end - start) %/% 2L, site)]
tw <- reg[gene_id %in% gt$gene_id, .(chr, s0 = pmax(tss - 50000L, 0L), s1 = tss + 50000L, gene_id)]
setkey(tw, chr, s0, s1)
w50 <- unique(foverlaps(sb, tw, by.x = c("chr", "s0", "s1"), nomatch = NULL)[, .(site, gene_id)])
links_by <- list(GREAT = sgl[tested_in_bulk == 1, .(site, gene_id)], within_50kb = w50)
kgroups <- c(paste0("cluster", 1:6), "stays", "other_KO", "no_KO_effect")
ccol <- setNames(paste0("sites_", gsub(",", "_", CLASSES)), CLASSES)
gg_all <- list(); gl_all <- list()
for (ln in names(links_by)) {
  gs <- merge(links_by[[ln]], si[, .(site, class)], by = "site")
  gw <- dcast(gs, gene_id ~ class, fun.aggregate = length, value.var = "site", drop = FALSE)
  setnames(gw, CLASSES, ccol[CLASSES])
  gg <- merge(gt[, .(gene_id, gene_name, group, domain_kb)], gw, by = "gene_id")
  gg[, n_sites := rowSums(.SD), .SDcols = unname(ccol)]
  stopifnot(!anyDuplicated(gg$gene_id), all(gg$n_sites > 0))
  if (ln == "GREAT") stopifnot(all(gg$n_sites == gt[match(gg$gene_id, gene_id), n_sites]))
  gg[mk, marker := i.cell_type, on = c(gene_name = "gene")]
  gg[, only_class := apply(as.matrix(.SD), 1, function(v) if (sum(v > 0) == 1) CLASSES[v > 0] else ""), .SDcols = unname(ccol)]
  gg[, link := ln]
  cov <- if (ln == "GREAT") "+ log2(1 + n_sites) + log10(domain_kb)" else "+ log2(1 + n_sites)"
  gl_all[[ln]] <- rbindlist(lapply(CLASSES, function(cl) {
    gg[, y := as.integer(get(ccol[[cl]]) > 0)]
    rbindlist(lapply(kgroups, function(kg) {
      gg[, x := as.integer(group == kg)]
      ft <- fisher.test(table(factor(gg$x, 1:0), factor(gg$y, 1:0)))
      fit <- glm(as.formula(paste("y ~ x", cov)), data = gg, family = binomial)
      cf <- summary(fit)$coefficients["x", ]
      data.table(link = ln, class = cl, ko_group = kg, genes = gg[x == 1 & y == 1, .N], genes_in_class = sum(gg$y),
                 pct_in_class = round(100 * gg[x == 1 & y == 1, .N] / sum(gg$y), 2),
                 pct_other_genes = round(100 * gg[x == 1 & y == 0, .N] / sum(gg$y == 0), 2),
                 genes_only = gg[x == 1 & only_class == cl, .N], genes_only_class = gg[only_class == cl, .N],
                 odds_ratio = signif(unname(ft$estimate), 3), p = signif(ft$p.value, 3),
                 adj_odds_ratio = signif(exp(cf[["Estimate"]]), 3), adj_p = signif(cf[["Pr(>|z|)"]], 3))
    }))
  }))
  gg[, c("x", "y") := NULL]
  gg_all[[ln]] <- gg
  cat(sprintf("%s: linked genes %d; with all sites in one class %d\n", ln, nrow(gg), sum(gg$only_class != "")))
}
gg <- rbindlist(gg_all, use.names = TRUE)
setcolorder(gg, c("link", "gene_id", "gene_name", "group", "marker", "n_sites", "only_class"))
setorder(gg, link, group, gene_name)
fwrite(gg, file.path(od, "gene_classes.tsv"), sep = "\t", quote = FALSE)
gl <- rbindlist(gl_all)
gl[, `:=`(fdr = signif(p.adjust(p, "BH"), 3), adj_fdr = signif(p.adjust(adj_p, "BH"), 3)), by = link]
fwrite(gl, file.path(od, "gene_class_links.tsv"), sep = "\t", quote = FALSE)
print(gl[adj_fdr < 0.1 | fdr < 0.05])
print(gg[gene_name %in% c("INS", "PAX6", "ARX", "PAX4", "MAFA", "MAFB", "PDX1", "NKX6-1", "ISL1", "NEUROG3")])
