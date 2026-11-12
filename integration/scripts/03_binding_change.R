#!/usr/bin/env Rscript
# Step 03: does ISL1 binding near each bulk KO-effect group's genes change over time, relative to all sites?
#   Sites      step 27 sites (14,431); a site belongs to every KO-effect group that has a gene linked to it in
#              step 01 (GREAT), so a site can be in more than one group. Groups from step 02.
#   Counts     step 27 featureCounts reads per site (the counts of its DESeq2 check), 9 control ISL1 libraries;
#              scaling from step 27's background size factors (DESeq2 on 15 kb bins), times the geometric mean of
#              the reads in sites so voom sees counts on their real scale.
#   Model      limma voom, ~ batch + day (batch = batch 3 vs earlier libraries; day D15 / D18 / D22);
#              contrasts D18 vs D15, D22 vs D15 and D22 vs D18.
#   Test       camera per group, competitive: the group's sites against all other sites. Main: limma's default
#              inter-site correlation (0.01); check: correlation estimated from the data (inter.gene.cor = NA).
#   Descriptive per group: median step 27 DiffBind (background) log2FC, share of gained sites (FDR < 0.05 up in
#              D18 vs D15 or D22 vs D15), Wilcoxon test of the DiffBind log2FC against all other sites.
# usage: Rscript 03_binding_change.R <step 27 isl1_diff dir> <integration results dir> <out_dir>
# Output in <out_dir>: camera.tsv, site_stats.tsv (per site: DiffBind and voom log2FC, FDR, gained),
#   site_groups.tsv (site, group)

suppressPackageStartupMessages({ library(data.table); library(limma) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
cnr <- args[1]; res <- args[2]; od <- args[3]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
TH <- 0.05

# ---- counts, libraries, scaling
cs <- fread(file.path(cnr, "counts_sites.tsv"))
m <- as.matrix(cs[, -1]); rownames(m) <- cs$name
smp <- fread(file.path(cnr, "samples.tsv"))
stopifnot(identical(colnames(m), smp$sample))
nf <- fread(file.path(cnr, "norm_factors.tsv"))[method == "deseq2" & normalization == "background"]
sf <- nf$norm_factor[match(colnames(m), nf$sample)]
stopifnot(!anyNA(sf))
lib <- sf * exp(mean(log(colSums(m))))
day <- factor(smp$day, levels = c("D15", "D18", "D22"))
batch <- factor(smp$batch, levels = c("earlier", "batch3"))
design <- model.matrix(~ batch + day)
cm <- makeContrasts(D18_vs_D15 = dayD18, D22_vs_D15 = dayD22, D22_vs_D18 = dayD22 - dayD18,
                    levels = make.names(colnames(design)))
rownames(cm) <- colnames(design)
contrasts <- setNames(lapply(colnames(cm), function(ct) cm[, ct]), colnames(cm))
v <- voom(m, design, lib.size = lib)
fit <- eBayes(contrasts.fit(lmFit(v, design), cm))

# ---- per-site statistics
st <- data.table(site = rownames(m))
for (ct in names(contrasts)) {
  tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
  st[, (paste0("voom_log2FC_", ct)) := tt$logFC]
  db <- fread(file.path(cnr, "tables", sprintf("diffbind_background_%s.tsv", ct)))
  st[db, (paste0("db_log2FC_", ct)) := i.log2FC, on = c(site = "name")]
  st[db, (paste0("db_FDR_", ct)) := i.FDR, on = c(site = "name")]
}
stopifnot(!anyNA(st$db_log2FC_D18_vs_D15), !anyNA(st$db_log2FC_D22_vs_D15))
st[, gained := (db_FDR_D18_vs_D15 < TH & db_log2FC_D18_vs_D15 > 0) | (db_FDR_D22_vs_D15 < TH & db_log2FC_D22_vs_D15 > 0)]
st[, lost := (db_FDR_D18_vs_D15 < TH & db_log2FC_D18_vs_D15 < 0) | (db_FDR_D22_vs_D15 < TH & db_log2FC_D22_vs_D15 < 0)]
fwrite(st, file.path(od, "site_stats.tsv"), sep = "\t", quote = FALSE)

# ---- site sets per KO-effect group
sg <- fread(file.path(res, "01_great", "site_gene.tsv"), sep = "\t", quote = "")
gt <- fread(file.path(res, "02_direct_targets", "gene_table.tsv"), sep = "\t", quote = "", na.strings = "NA")
sg[gt, group := i.group, on = "gene_id"]
sgr <- unique(sg[!is.na(group), .(site, group)])
lv <- c("no_KO_effect", paste0("cluster", 1:6), "stays", "other_KO")
sgr[, group := factor(group, levels = lv)]
setorder(sgr, group, site)
fwrite(sgr, file.path(od, "site_groups.tsv"), sep = "\t", quote = FALSE)
idx <- lapply(setNames(lv, lv), function(l) match(sgr[group == l, site], rownames(v)))
stopifnot(!anyNA(unlist(idx)))

# ---- camera and descriptive statistics
out <- rbindlist(lapply(names(contrasts), function(ct) {
  cam <- camera(v, idx, design, contrast = contrasts[[ct]], sort = FALSE)
  cam_est <- camera(v, idx, design, contrast = contrasts[[ct]], inter.gene.cor = NA, sort = FALSE)
  lfc <- st[[paste0("db_log2FC_", ct)]]
  rbindlist(lapply(lv, function(l) {
    i <- idx[[l]]
    data.table(contrast = ct, group = l, n_sites = length(i),
               median_log2FC = round(median(lfc[i]), 3), median_log2FC_other = round(median(lfc[-i]), 3),
               pct_gained = round(100 * mean(st$gained[i]), 1), pct_gained_other = round(100 * mean(st$gained[-i]), 1),
               wilcoxon_p = signif(wilcox.test(lfc[i], lfc[-i])$p.value, 3),
               camera_direction = cam[l, "Direction"], camera_p = signif(cam[l, "PValue"], 3),
               camera_fdr = signif(cam[l, "FDR"], 3),
               est_correlation = signif(cam_est[l, "Correlation"], 3), est_direction = cam_est[l, "Direction"],
               est_p = signif(cam_est[l, "PValue"], 3), est_fdr = signif(cam_est[l, "FDR"], 3))
  }))
}))
fwrite(out, file.path(od, "camera.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("all sites: %d; gained %d, lost %d\n", nrow(st), sum(st$gained), sum(st$lost)))
for (ct in names(contrasts)) cat(sprintf("%s: DiffBind vs voom log2FC, Pearson r = %.3f\n", ct,
  cor(st[[paste0("db_log2FC_", ct)]], st[[paste0("voom_log2FC_", ct)]])))
print(out[, .(contrast, group, n_sites, median_log2FC, median_log2FC_other, pct_gained, pct_gained_other,
              camera_direction, camera_fdr, est_correlation, est_fdr)])
