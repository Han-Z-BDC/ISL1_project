#!/usr/bin/env Rscript
# Step 09: overlap of the per-day ISL1 peak calls (D15 / D18 / D22), checked against the count-based tests.
#   Presence   each of the 14,431 sites (CUT&RUN step 06 sites after the step 07 IgG greylist) carries the days on which it lies in
#              that day's consensus peaks (macs3_ge2: MACS3 q 0.01 peaks in >= 2 of the day's 3 libraries, CUT&RUN step 05).
#              The 7 combinations of days are the overlap classes.
#   Libraries  each site is also overlapped with each library's own MACS3 peaks: number of a day's libraries that call
#              it (0-3), and whether the batch 3 library does. Per library: peaks, sites called, library size and
#              share of reads in the 14,431 sites (DiffBind counts, CUT&RUN step 07), to show how call rates follow the library.
#   Counts     per class: batch-corrected VST per day and the step 08 LRT / binding-pattern group; DiffBind
#              (background normalization, CUT&RUN step 07 via step 03) log2FC and FDR per contrast.
#   Check      per contrast A -> B: sites called on B but not A, on A but not B, on both, on neither; for each, the
#              median DiffBind log2FC and the share with FDR < 0.05 up / down. If presence differences reflect binding
#              changes, the sites that appear should mostly go up in the counts.
#   Genes      linked genes tested in the bulk RNA-seq (step 02), two links as in step 08: GREAT domains (main) and every
#              bulk-tested gene with its TSS within 50 kb of the site centre (check). A gene is in a class if at least one
#              of its linked sites is (a gene can be in several classes); "only" genes have all their sites in the class.
#              Per class and KO-effect group: genes in the class against all other linked genes; Fisher test, and
#              logistic regression in_class ~ in_KO_group + log2(1 + n_sites) (+ log10(domain_kb) for GREAT), since
#              genes with more sites are more often in any class and more often KO-affected; BH per link.
# usage: Rscript 09_peak_overlap.R <CUT&RUN diff dir> <integration results dir> <markers.tsv> <out_dir>
# Output in <out_dir>: sites_days.tsv (per site), classes.tsv (per class), class_days.tsv (per class and day),
#   presence_vs_counts.tsv (per contrast and presence category), libraries.tsv (per library), gene_classes.tsv (per
#   link and gene: KO group, marker, sites per class), gene_class_links.tsv (per link, class and KO group)

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
cnr <- args[1]; res <- args[2]; mfile <- args[3]; od <- args[4]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
days <- c("D15", "D18", "D22"); TH <- 0.05
CLASSES <- c("D15,D18,D22", "D18,D22", "D15,D18", "D15,D22", "D18", "D22", "D15")
root <- dirname(cnr)

si <- fread(file.path(cnr, "sites_used.bed"), header = FALSE,
            col.names = c("chr", "start", "end", "site", "score", "strand", "summit_day", "days_called", "n_days"))
stopifnot(all(si$days_called %in% CLASSES))
for (d in days) si[, (paste0("called_", d)) := grepl(d, days_called, fixed = TRUE)]

# ---- per-library calls
smp <- fread(file.path(cnr, "samples.tsv"))
smp[, peaks := file.path(root, "macs3", paste0(sample, "_peaks.filtered.narrowPeak"))]   # CUT&RUN step 04
stopifnot(all(file.exists(smp$peaks)))
sk <- si[, .(chr, s = start + 1L, e = end, site)]
setkey(sk, chr, s, e)
lib <- vector("list", nrow(smp))
for (i in seq_len(nrow(smp))) {
  p <- fread(smp$peaks[i], header = FALSE, select = 1:3, col.names = c("chr", "s", "e"))
  p[, s := s + 1L]
  setkey(p, chr, s, e)
  hit <- unique(foverlaps(sk, p, nomatch = NULL)$site)
  si[, (smp$sample[i]) := site %in% hit]
  lib[[i]] <- data.table(sample = smp$sample[i], day = smp$day[i], batch = smp$batch[i], peaks = nrow(p),
                         sites_called = length(hit))
}
for (d in days) {
  si[, (paste0("libs_", d)) := rowSums(.SD), .SDcols = smp[day == d, sample]]
  si[, (paste0("b3_", d)) := get(smp[day == d & batch == "batch3", sample])]
}
cat("consensus call vs per-library calls (share of sites called on the day with >= 2 libraries calling; share of",
    "sites not called with <= 1):\n")
for (d in days) {
  cd <- si[[paste0("called_", d)]]; nl <- si[[paste0("libs_", d)]]
  cat(sprintf("  %s: %.1f%% / %.1f%%\n", d, 100 * mean(nl[cd] >= 2), 100 * mean(nl[!cd] <= 1)))
}

nf <- fread(file.path(cnr, "norm_factors.tsv"))[method == "diffbind" & normalization == "background"][, .(sample, lib_size, norm_factor)]
cr <- fread(file.path(cnr, "counts_diffbind.tsv"))
rs <- colSums(as.matrix(cr[, smp$sample, with = FALSE]))
libs <- rbindlist(lib)[nf, on = "sample", nomatch = NULL]
stopifnot(nrow(libs) == nrow(smp))
libs[, `:=`(reads_in_sites = rs[sample], pct_reads_in_sites = round(100 * rs[sample] / lib_size, 2),
            pct_sites_called = round(100 * sites_called / nrow(si), 1))]
libs[, norm_factor := signif(norm_factor, 4)]
setcolorder(libs, c("sample", "day", "batch", "lib_size", "norm_factor", "peaks", "sites_called", "pct_sites_called",
                    "reads_in_sites", "pct_reads_in_sites"))
print(libs)

# ---- counts: step 08 profile and group, step 03 DiffBind results
sc <- fread(file.path(res, "08_site_clusters", "site_clusters.tsv"))[
  , .(site, lrt_p, varying, vst_D15, vst_D18, vst_D22, group = k4)]
ss <- fread(file.path(res, "03_binding_change", "site_stats.tsv"))
ss <- ss[, c("site", grep("^db_", names(ss), value = TRUE), "gained", "lost"), with = FALSE]
n0 <- nrow(si)
si <- merge(merge(si, sc, by = "site"), ss, by = "site")
stopifnot(nrow(si) == n0, !anyNA(si$group))
si[, class := factor(days_called, levels = CLASSES)]
setorder(si, chr, start)

fwrite(si[, c("site", "chr", "start", "end", "summit_day", "days_called", paste0("libs_", days), paste0("b3_", days),
              "lrt_p", "varying", "group", paste0("vst_", days), grep("^db_", names(si), value = TRUE), "gained",
              "lost"), with = FALSE],
       file.path(od, "sites_days.tsv"), sep = "\t", quote = FALSE)

# ---- per class
cl <- si[, .(sites = .N, pct_sites = round(100 * .N / nrow(si), 1),
             vst_D15 = round(median(vst_D15), 3), vst_D18 = round(median(vst_D18), 3), vst_D22 = round(median(vst_D22), 3),
             pct_varying = round(100 * mean(varying), 1),
             pct_group1 = round(100 * mean(group == 1), 1), pct_group2 = round(100 * mean(group == 2), 1),
             pct_group3 = round(100 * mean(group == 3), 1), pct_group4 = round(100 * mean(group == 4), 1),
             log2FC_D18_vs_D15 = round(median(db_log2FC_D18_vs_D15), 3),
             log2FC_D22_vs_D15 = round(median(db_log2FC_D22_vs_D15), 3),
             log2FC_D22_vs_D18 = round(median(db_log2FC_D22_vs_D18), 3),
             pct_gained = round(100 * mean(gained), 1), pct_lost = round(100 * mean(lost), 1)), by = class][order(class)]
fwrite(cl, file.path(od, "classes.tsv"), sep = "\t", quote = FALSE)
print(cl)

# ---- per class and day: library support on the days a site is and is not called
cd <- rbindlist(lapply(days, function(d) {
  si[, .(day = d, called = get(paste0("called_", d))[1], sites = .N,
         pct_libs0 = round(100 * mean(get(paste0("libs_", d)) == 0), 1),
         pct_libs1 = round(100 * mean(get(paste0("libs_", d)) == 1), 1),
         pct_libs2 = round(100 * mean(get(paste0("libs_", d)) == 2), 1),
         pct_libs3 = round(100 * mean(get(paste0("libs_", d)) == 3), 1),
         pct_batch3_lib = round(100 * mean(get(paste0("b3_", d))), 1),
         median_vst = round(median(get(paste0("vst_", d))), 3)), by = class]
}))
setorder(cd, class, day)
fwrite(cd, file.path(od, "class_days.tsv"), sep = "\t", quote = FALSE)
print(cd)

# ---- presence against counts, per contrast
ct <- list(c("D15", "D18"), c("D15", "D22"), c("D18", "D22"))
pc <- rbindlist(lapply(ct, function(ab) {
  a <- ab[1]; b <- ab[2]; cn <- paste0(b, "_vs_", a)
  x <- data.table(inA = si[[paste0("called_", a)]], inB = si[[paste0("called_", b)]],
                  fc = si[[paste0("db_log2FC_", cn)]], fdr = si[[paste0("db_FDR_", cn)]])
  x[, presence := fifelse(inB & !inA, paste0("called on ", b, ", not ", a),
                  fifelse(inA & !inB, paste0("called on ", a, ", not ", b),
                  fifelse(inA & inB, "called on both", "called on neither (third day only)")))]
  x[, .(contrast = cn, sites = .N, median_log2FC = round(median(fc), 3), pct_log2FC_positive = round(100 * mean(fc > 0), 1),
        pct_up_fdr = round(100 * mean(fdr < TH & fc > 0), 1), pct_down_fdr = round(100 * mean(fdr < TH & fc < 0), 1)),
    by = presence]
}))
pc[, rk := fifelse(grepl("^called on D", presence), fifelse(substr(presence, 11, 13) == sub(".*_vs_", "", contrast), 2L, 1L),
                   fifelse(presence == "called on both", 3L, 4L))]
setorder(pc, contrast, rk)
pc[, rk := NULL]
setcolorder(pc, c("contrast", "presence"))
fwrite(pc, file.path(od, "presence_vs_counts.tsv"), sep = "\t", quote = FALSE)
print(pc)
fwrite(libs, file.path(od, "libraries.tsv"), sep = "\t", quote = FALSE)

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
