#!/usr/bin/env Rscript
# Step 10: pre-specified test of ISL1 binding at beta against alpha marker genes over time.
#
# PRE-SPECIFICATION (fixed 2026-11-15, before the unified CUT&RUN time course is analysed; decisions by the user)
#   Hypothesis   ISL1 binding at beta-program genes gains on binding at alpha-program genes late in differentiation.
#   Data         primary test: the unified CUT&RUN time course (cutnrun_ISL1_timecourse; Ctr ISL1, D15 / D18 / D22).
#                The current 9 libraries (cutnrun_ISL1_KO_histone step 27) give a DRY RUN only: step 06 already
#                examined beta vs alpha marker sites in these libraries, so a test on them is not pre-specified.
#   Gene sets    integration/markers.tsv as it stands (28 alpha, 34 beta human islet marker genes).
#   Sites        every ISL1 site whose centre lies within +-50 kb of the TSS of a marker gene of the set (TSS from
#                step 01 gene_regions.bed). A site near both an alpha and a beta gene counts in both sets.
#   Signal       per library and site: DiffBind reads / DiffBind background size factor (step 00).
#   Outcome      per library: y = log2(median signal over the beta sites / median signal over the alpha sites).
#   Model        lm(y ~ batch + day) when the libraries come from more than one batch, else lm(y ~ day).
#   Primary      contrast D22 - D18 of the day effect; H1: > 0 (beta gains on alpha); one-sided t test, alpha 0.05.
#   Secondary    (reported, not the claim) D22 - D15 and D18 - D15 with the same sites; all three contrasts with
#                GREAT-linked sites (step 01); site-level check: DiffBind log2FC (step 03) of the beta against the
#                alpha sites per contrast, Wilcoxon rank-sum (descriptive; sites are not independent replicates).
#
# usage: Rscript 10_beta_alpha_contrast.R <isl1_diff dir> <integration results dir> <markers.tsv> <mode> <out_dir>
#   mode: dry_run (current libraries) or primary (unified time course); written into every output table.
# Output in <out_dir>: sets.tsv (per link and set: genes with sites, sites), set_sites.tsv (site membership),
#   library_values.tsv (per link and library: medians and y), contrasts.tsv (per link and contrast),
#   site_level.tsv (per contrast: DiffBind log2FC of beta vs alpha sites)

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
cnr <- args[1]; res <- args[2]; mfile <- args[3]; mode <- args[4]; od <- args[5]
stopifnot(mode %in% c("dry_run", "primary"))
dir.create(od, showWarnings = FALSE, recursive = TRUE)
days <- c("D15", "D18", "D22"); WIN <- 50000L

# ---- signal
cs <- fread(file.path(res, "00_diffbind_counts", "counts_raw.tsv"))
m <- as.matrix(cs[, -1]); rownames(m) <- cs$name
smp <- fread(file.path(cnr, "samples.tsv"))
stopifnot(identical(colnames(m), smp$sample), all(smp$day %in% days))
nf <- fread(file.path(res, "00_diffbind_counts", "norm_factors.tsv"))
sf <- nf$norm_factor[match(colnames(m), nf$sample)]
stopifnot(!anyNA(sf))
nm <- sweep(m, 2, sf, "/")

# ---- site sets
mk <- fread(mfile, skip = "cell_type")
stopifnot(all(mk$cell_type %in% c("alpha", "beta")))
reg <- fread(file.path(res, "01_great", "gene_regions.bed"), header = FALSE,
             col.names = c("chr", "start", "end", "gene_id", "gene_name", "strand", "tss", "basal_start",
                           "basal_end", "gene_type", "tested"))
miss <- setdiff(mk$gene, reg$gene_name)
if (length(miss)) cat("marker genes not in the annotation:", paste(miss, collapse = ", "), "\n")
si <- fread(file.path(cnr, "sites_used.bed"), header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "site"))
stopifnot(setequal(si$site, rownames(nm)))
si[, centre := start + (end - start) %/% 2L]
mreg <- merge(reg[gene_name %in% mk$gene, .(chr, tss, gene_name)], mk, by.x = "gene_name", by.y = "gene")
w50 <- mreg[si, on = .(chr), allow.cartesian = TRUE, nomatch = NULL][abs(centre - tss) <= WIN,
                                                                      .(site, gene_name, cell_type)]
sg <- fread(file.path(res, "01_great", "site_gene.tsv"), sep = "\t", quote = "")
grt <- merge(unique(sg[, .(site, gene_name)]), mk, by.x = "gene_name", by.y = "gene")[, .(site, gene_name, cell_type)]
links <- list(within_50kb = unique(w50), GREAT = unique(grt))
memb <- rbindlist(lapply(names(links), function(ln) links[[ln]][, .(link = ln, site, gene_name, cell_type)]))
setorder(memb, link, cell_type, gene_name, site)
fwrite(memb, file.path(od, "set_sites.tsv"), sep = "\t", quote = FALSE)
sets <- memb[, .(genes_with_sites = uniqueN(gene_name), genes_in_set = mk[cell_type == .BY$cell_type, .N],
                 sites = uniqueN(site), genes = paste(sort(unique(gene_name)), collapse = ", ")), by = .(link, cell_type)]
shared <- memb[, .(n = uniqueN(cell_type)), by = .(link, site)][n == 2, .N, by = link]
sets[, mode := mode]
fwrite(sets, file.path(od, "sets.tsv"), sep = "\t", quote = FALSE)
print(sets[, -"genes"]); cat("sites in both sets:\n"); print(shared)

# ---- per library outcome and model
batches <- unique(smp$batch)
lv <- rbindlist(lapply(names(links), function(ln) {
  a <- unique(links[[ln]][cell_type == "alpha", site]); b <- unique(links[[ln]][cell_type == "beta", site])
  ma <- apply(nm[a, , drop = FALSE], 2, median); mb <- apply(nm[b, , drop = FALSE], 2, median)
  data.table(link = ln, sample = smp$sample, day = smp$day, batch = smp$batch, alpha_sites = length(a),
             beta_sites = length(b), median_alpha = round(ma, 3), median_beta = round(mb, 3),
             log2_beta_alpha = round(log2(mb / ma), 4), mode = mode)
}))
fwrite(lv, file.path(od, "library_values.tsv"), sep = "\t", quote = FALSE)

CT <- list(D22_vs_D18 = c(D18 = -1, D22 = 1), D22_vs_D15 = c(D18 = 0, D22 = 1), D18_vs_D15 = c(D18 = 1, D22 = 0))
cts <- rbindlist(lapply(names(links), function(ln) {
  d <- lv[link == ln]
  d <- data.frame(y = d$log2_beta_alpha, day = factor(d$day, levels = days), batch = factor(d$batch))
  fit <- if (length(batches) > 1) lm(y ~ batch + day, data = d) else lm(y ~ day, data = d)
  cf <- coef(fit); V <- vcov(fit)
  rbindlist(lapply(names(CT), function(cn) {
    w <- setNames(numeric(length(cf)), names(cf)); w["dayD18"] <- CT[[cn]][["D18"]]; w["dayD22"] <- CT[[cn]][["D22"]]
    est <- sum(w * cf); se <- sqrt(drop(t(w) %*% V %*% w)); tt <- est / se; df <- fit$df.residual
    data.table(link = ln, contrast = cn,
               role = if (ln == "within_50kb" && cn == "D22_vs_D18") "primary" else "secondary",
               model = if (length(batches) > 1) "y ~ batch + day" else "y ~ day",
               estimate = round(est, 4), se = round(se, 4), df = df, t = round(tt, 3),
               p_one_sided = signif(pt(tt, df, lower.tail = FALSE), 3), p_two_sided = signif(2 * pt(-abs(tt), df), 3),
               mode = mode)
  }))
}))
fwrite(cts, file.path(od, "contrasts.tsv"), sep = "\t", quote = FALSE)
print(cts)

# ---- site-level check (descriptive)
ss <- fread(file.path(res, "03_binding_change", "site_stats.tsv"))
slv <- rbindlist(lapply(names(links), function(ln) {
  a <- unique(links[[ln]][cell_type == "alpha", site]); b <- unique(links[[ln]][cell_type == "beta", site])
  rbindlist(lapply(names(CT), function(cn) {
    fc <- setNames(ss[[paste0("db_log2FC_", cn)]], ss$site)
    data.table(link = ln, contrast = cn, alpha_sites = length(a), beta_sites = length(b),
               median_log2FC_alpha = round(median(fc[a]), 3), median_log2FC_beta = round(median(fc[b]), 3),
               difference = round(median(fc[b]) - median(fc[a]), 3),
               p_wilcoxon = signif(wilcox.test(fc[b], fc[a])$p.value, 3), mode = mode)
  }))
}))
fwrite(slv, file.path(od, "site_level.tsv"), sep = "\t", quote = FALSE)
print(slv)
