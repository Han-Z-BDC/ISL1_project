#!/usr/bin/env Rscript
# DESeq2 differential analysis for one histone mark (counts from scripts/16; agreed 2026-11-09).
# Samples: day x genotype x replicate. Regions: those passing the step 16 pre-filter.
# Main analysis, DESeq2 median-of-ratios size factors:
#   LRTs on ~ rep + genotype + day + genotype:day
#     time         reduced ~ rep + genotype          (changes over time in either genotype)
#     genotype     reduced ~ rep + day               (KO differs from Ctr at any day)
#     interaction  reduced ~ rep + genotype + day    (KO effect depends on the day)
#   Wald contrasts on the equivalent ~ rep + group (6 groups), log2FC shrunk with ashr:
#     time within each genotype (D18 vs D15, D22 vs D18, D22 vs D15), KO vs Ctr at each day,
#     and the interaction contrasts (KO - Ctr at day b) - (KO - Ctr at day a)
#   significant: padj < PADJ and |shrunk log2FC| >= LFC (Wald); padj < PADJ (LRT)
#   KO patterns per region from the three KO vs Ctr contrasts (consistent / two days / one day /
#   mixed), with the interaction LRT padj alongside
#   time clusters: regions significant in the time LRT, k-means (K, seed 1) on z-scored VST means
#   of the 6 groups
# Sensitivity checks (summaries only):
#   libsize   size factors = all fragments featureCounts saw per library (geometric mean 1)
#   drop_<s>  without library DROP_SAMPLE (default the weakest H3K27ac library, D18 H3CTR rep2)
# Output in <diff_dir>/deseq2/: tables/*.tsv (all regions, annotated), summary.tsv,
#   ko_patterns.tsv, time_clusters.tsv, size_factors.tsv, *.png, report.html
# usage: Rscript histone_deseq2.R <results/diff/<mark>> [drop_sample]

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2)
dd <- args[1]
mark <- basename(dd)
drop_sample <- if (length(args) == 2) args[2] else NA
out <- file.path(dd, "deseq2"); tab <- file.path(out, "tables")
dir.create(tab, showWarnings = FALSE, recursive = TRUE)
PADJ <- 0.05; LFC <- 1; K <- 6
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
register(MulticoreParam(ncores))
set.seed(1)

# ---- inputs
fc <- fread(file.path(dd, "counts.txt"), skip = "Geneid")
m <- as.matrix(fc[, 7:ncol(fc)]); rownames(m) <- fc$Geneid
colnames(m) <- sub("\\.filtered\\.bam$", "", basename(colnames(m)))
samples <- fread(file.path(dd, "samples.tsv"))
stopifnot(identical(colnames(m), samples$sample))
summ <- fread(file.path(dd, "counts.txt.summary"))
total <- colSums(as.matrix(summ[, -1])); names(total) <- sub("\\.filtered\\.bam$", "", basename(names(total)))
flt <- fread(file.path(dd, "filter.tsv"))
m <- m[flt[pass == 1, region], ]
regions <- fread(file.path(dd, "regions.bed"), header = FALSE,
                 col.names = c("chr", "start", "end", "region", "n_groups", "groups"))
regions <- regions[region %in% rownames(m)]
cat(sprintf("%s: %d libraries, %d regions after the pre-filter\n", mark, ncol(m), nrow(m)))

days <- intersect(c("D9", "D15", "D18", "D22"), unique(samples$day))
genos <- c("H3CTR", "ISL1KO")
glev <- c(paste0(days, "_H3CTR"), paste0(days, "_ISL1KO"))
cd <- data.frame(row.names = samples$sample,
                 day = factor(samples$day, levels = days),
                 genotype = factor(samples$genotype, levels = genos),
                 rep = factor(samples$replicate),
                 group = factor(samples$group, levels = glev))

# ---- annotation: nearest gene to the region centre
gr <- GRanges(regions$chr, IRanges(as.integer((regions$start + regions$end) / 2) + 1, width = 1), region = regions$region)
ann <- as.data.table(as.data.frame(annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                                annoDb = "org.Hs.eg.db", verbose = FALSE)))
regions <- merge(regions, ann[, .(region, annotation, SYMBOL, distanceToTSS)], by = "region", all.x = TRUE)
regions[, annotation := sub(" \\(.*", "", annotation)]

# ---- contrasts on ~ rep + group (vectors over resultsNames)
gv <- function(rn, g) {
  v <- setNames(numeric(length(rn)), rn)
  if (g != glev[1]) v[paste0("group_", g, "_vs_", glev[1])] <- 1
  v
}
pairs_ <- list()
for (g in genos) for (ab in list(c("D18", "D15"), c("D22", "D18"), c("D22", "D15")))
  if (all(ab %in% days)) pairs_[[sprintf("time_%s_%s_vs_%s", g, ab[1], ab[2])]] <- list(plus = paste0(ab[1], "_", g), minus = paste0(ab[2], "_", g))
for (d in days) pairs_[[sprintf("KO_vs_Ctr_%s", d)]] <- list(plus = paste0(d, "_ISL1KO"), minus = paste0(d, "_H3CTR"))
for (ab in list(c("D18", "D15"), c("D22", "D15"), c("D22", "D18")))
  if (all(ab %in% days)) pairs_[[sprintf("interaction_%s_vs_%s", ab[1], ab[2])]] <-
    list(plus = c(paste0(ab[1], "_ISL1KO"), paste0(ab[2], "_H3CTR")), minus = c(paste0(ab[1], "_H3CTR"), paste0(ab[2], "_ISL1KO")))
cvec <- function(rn, p) Reduce(`+`, lapply(p$plus, gv, rn = rn)) - Reduce(`+`, lapply(p$minus, gv, rn = rn))

wald <- function(dds, name, shrink = TRUE) {
  rn <- resultsNames(dds)
  r <- results(dds, contrast = cvec(rn, pairs_[[name]]), alpha = PADJ, parallel = TRUE)
  s <- if (shrink) lfcShrink(dds, res = r, type = "ashr", quiet = TRUE)$log2FoldChange else r$log2FoldChange
  data.table(region = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange, log2FC_shrunk = s,
             lfcSE = r$lfcSE, pvalue = r$pvalue, padj = r$padj)
}
sig_dir <- function(d) ifelse(!is.na(d$padj) & d$padj < PADJ & abs(d$log2FC_shrunk) >= LFC, sign(d$log2FC_shrunk), 0)
summ_row <- function(d, name, analysis) {
  s <- sig_dir(d)
  data.table(analysis = analysis, test = name, tested = sum(!is.na(d$padj)), up = sum(s > 0), down = sum(s < 0),
             median_log2FC = round(median(d$log2FC, na.rm = TRUE), 3))
}
fit_group <- function(mm, cdd, sf = NULL) {
  dds <- DESeqDataSetFromMatrix(mm, cdd, design = ~ rep + group)
  if (!is.null(sf)) sizeFactors(dds) <- sf
  DESeq(dds, quiet = TRUE, parallel = TRUE)
}

# ---- main: Wald contrasts
dds <- fit_group(m, cd)
saveRDS(dds, file.path(out, "dds_group.rds"))
summary_rows <- list(); main <- list()
for (nm in names(pairs_)) {
  d <- wald(dds, nm)
  main[[nm]] <- d
  fwrite(merge(regions, d, by = "region")[order(padj, na.last = TRUE)], file.path(tab, paste0(nm, ".tsv")), sep = "\t")
  summary_rows[[length(summary_rows) + 1]] <- summ_row(d, nm, "main")
  cat(sprintf("Wald %s: %d up, %d down\n", nm, sum(sig_dir(d) > 0), sum(sig_dir(d) < 0)))
}

# ---- main: LRTs
full <- ~ rep + genotype + day + genotype:day
lrts <- list(time = ~ rep + genotype, genotype = ~ rep + day, interaction = ~ rep + genotype + day)
lrt_res <- list()
for (nm in names(lrts)) {
  dl <- DESeq(DESeqDataSetFromMatrix(m, cd, design = full), test = "LRT", reduced = lrts[[nm]], quiet = TRUE, parallel = TRUE)
  r <- results(dl, alpha = PADJ)
  d <- data.table(region = rownames(r), baseMean = r$baseMean, LRT_stat = r$stat, pvalue = r$pvalue, padj = r$padj)
  lrt_res[[nm]] <- d
  fwrite(merge(regions, d, by = "region")[order(padj, na.last = TRUE)], file.path(tab, paste0("LRT_", nm, ".tsv")), sep = "\t")
  n <- sum(!is.na(d$padj) & d$padj < PADJ)
  summary_rows[[length(summary_rows) + 1]] <- data.table(analysis = "main", test = paste0("LRT_", nm),
    tested = sum(!is.na(d$padj)), up = NA_integer_, down = NA_integer_, median_log2FC = NA_real_, significant = n)
  cat(sprintf("LRT %s: %d significant\n", nm, n))
}

# ---- KO patterns
ko <- grep("^KO_vs_Ctr_", names(pairs_), value = TRUE)
S <- sapply(ko, function(nm) sig_dir(main[[nm]])); rownames(S) <- main[[ko[1]]]$region
L <- sapply(ko, function(nm) main[[nm]]$log2FC_shrunk)
lab <- apply(S, 1, function(s) {
  n <- sum(s != 0); dd_ <- sub("KO_vs_Ctr_", "", ko[s != 0])
  if (n == 0) return("none")
  if (any(s > 0) && any(s < 0)) return("mixed")
  dir <- if (any(s > 0)) "up" else "down"
  if (n == length(ko)) return(paste0("all_days_", dir))
  paste0(paste(dd_, collapse = "+"), "_only_", dir)
})
pat <- data.table(region = rownames(S), pattern = lab, S, L)
setnames(pat, c("region", "pattern", paste0("sig_", sub("KO_vs_Ctr_", "", ko)), paste0("lfc_", sub("KO_vs_Ctr_", "", ko))))
pat <- merge(pat, lrt_res$interaction[, .(region, interaction_padj = padj)], by = "region")
pat <- merge(regions, pat[pattern != "none"], by = "region")
fwrite(pat[order(pattern, interaction_padj)], file.path(out, "ko_patterns.tsv"), sep = "\t")
pat_summ <- data.table(table(lab))[order(-N)]; setnames(pat_summ, c("pattern", "regions"))
pat_summ[, interaction_padj_lt_0.05 := sapply(pattern, function(p) pat[pattern == p, sum(interaction_padj < PADJ, na.rm = TRUE)])]
fwrite(pat_summ, file.path(out, "ko_pattern_summary.tsv"), sep = "\t")
print(pat_summ)

# ---- time clusters
vsd <- vst(dds, blind = FALSE)
v <- assay(vsd)
tsig <- lrt_res$time[!is.na(padj) & padj < PADJ, region]
gm <- sapply(glev, function(g) rowMeans(v[tsig, cd$group == g, drop = FALSE]))
z <- t(scale(t(gm))); z[is.na(z)] <- 0
km <- kmeans(z, centers = K, nstart = 25, iter.max = 100)
cl_ord <- order(apply(km$centers, 1, which.max))
cl <- match(km$cluster, cl_ord)
tc <- merge(regions, data.table(region = tsig, cluster = cl, round(gm, 3)), by = "region")
fwrite(tc[order(cluster)], file.path(out, "time_clusters.tsv"), sep = "\t")
cent <- t(sapply(seq_len(K), function(k) colMeans(z[cl == k, , drop = FALSE])))
png(file.path(out, "time_clusters.png"), width = 1500, height = 800, res = 130)
par(mfrow = c(2, ceiling(K / 2)), mar = c(4, 4, 3, 1))
for (k in seq_len(K)) {
  plot(NA, xlim = c(1, length(days)), ylim = range(cent), xaxt = "n", xlab = "", ylab = "mean z (VST)",
       main = sprintf("cluster %d (n = %s)", k, format(sum(cl == k), big.mark = ",")))
  axis(1, at = seq_along(days), labels = days)
  lines(seq_along(days), cent[k, paste0(days, "_H3CTR")], type = "b", pch = 16, col = "#4C72B0", lwd = 2)
  lines(seq_along(days), cent[k, paste0(days, "_ISL1KO")], type = "b", pch = 17, col = "#C44E52", lwd = 2)
  if (k == 1) legend("topleft", c("H3CTR", "ISL1KO"), col = c("#4C72B0", "#C44E52"), pch = c(16, 17), lwd = 2, bty = "n", cex = 0.8)
}
dev.off()

# ---- MA plots
ma_png <- function(res, file, title) {
  n <- length(res); nc <- 3; nr <- ceiling(n / nc)
  png(file, width = 1500, height = 430 * nr, res = 120)
  par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
  for (nm in names(res)) {
    d <- res[[nm]]; s <- sig_dir(d)
    plot(log10(d$baseMean + 1), d$log2FC_shrunk, pch = ".", col = ifelse(s > 0, "#C44E52", ifelse(s < 0, "#4C72B0", "grey60")),
         xlab = "log10 mean count", ylab = "shrunk log2FC", main = sprintf("%s\nup %d, down %d", nm, sum(s > 0), sum(s < 0)), cex.main = 0.9)
    abline(h = 0, col = "black")
  }
  dev.off()
}
ma_png(main, file.path(out, "ma_main.png"), "")

# ---- sensitivity: library-size normalization
sf_lib <- total[colnames(m)] / exp(mean(log(total[colnames(m)])))
dds_l <- fit_group(m, cd, sf_lib)
sens <- grep("^(KO_vs_Ctr|time_H3CTR)", names(pairs_), value = TRUE)
lib <- list()
for (nm in sens) { lib[[nm]] <- wald(dds_l, nm); summary_rows[[length(summary_rows) + 1]] <- summ_row(lib[[nm]], nm, "libsize") }
ma_png(lib, file.path(out, "ma_libsize.png"), "")
fwrite(data.table(sample = colnames(m), group = as.character(cd$group), fragments = total[colnames(m)],
                  counted_in_regions = colSums(m), size_factor_median_of_ratios = sizeFactors(dds), size_factor_libsize = sf_lib),
       file.path(out, "size_factors.tsv"), sep = "\t")

# ---- sensitivity: drop one library
if (!is.na(drop_sample) && drop_sample %in% colnames(m)) {
  keep <- colnames(m) != drop_sample
  cd2 <- droplevels(cd[keep, ]); cd2$group <- factor(as.character(cd2$group), levels = glev)
  dds_d <- fit_group(m[, keep], cd2)
  dgrp <- as.character(cd[drop_sample, "group"])
  for (nm in names(pairs_)) {
    p <- pairs_[[nm]]
    if (!(dgrp %in% c(p$plus, p$minus))) next
    d <- wald(dds_d, nm)
    a <- sig_dir(main[[nm]]); b <- sig_dir(d)
    r <- summ_row(d, nm, paste0("drop_", drop_sample))
    r[, kept_from_main := sprintf("%d of %d", sum(a != 0 & a == b), sum(a != 0))]
    summary_rows[[length(summary_rows) + 1]] <- r
  }
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(out, "summary.tsv"), sep = "\t")
print(summary_dt)

# ---- report
img <- function(f) sprintf('<img src="data:image/png;base64,%s">', base64_enc(readBin(f, "raw", file.info(f)$size)))
htab <- function(d) {
  d <- as.data.frame(d); d[is.na(d)] <- ""
  paste0("<table><tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>",
         paste0(apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")), collapse = ""), "</table>")
}
cl_tab <- data.table(cluster = seq_len(K), regions = as.integer(table(factor(cl, levels = seq_len(K)))), round(cent, 2))
# (paste0, not sprintf: sprintf limits each inserted string to 8192 bytes)
html <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>', mark, ' differential</title>
<style>body{font-family:sans-serif;max-width:1200px;margin:24px auto;padding:0 16px}table{border-collapse:collapse;font-size:12px;margin-bottom:16px}
td,th{border:1px solid #ccc;padding:3px 6px;text-align:right}img{max-width:100%}</style></head><body>
<h1>', mark, ': differential analysis (DESeq2)</h1>
<p>', ncol(m), ' libraries, ', nrow(m), ' regions (union of replicated peaks, pre-filtered). Main analysis: median-of-ratios size factors,
~ rep + genotype + day + genotype:day. Wald: padj &lt; ', PADJ, ' and |ashr-shrunk log2FC| &ge; ', LFC, '. LRT: padj &lt; ', PADJ, '.</p>
<h2>Summary</h2>', htab(summary_dt), '
<h2>KO vs Ctr patterns (main)</h2><p>From the three per-day KO vs Ctr contrasts; interaction LRT significant count alongside.</p>', htab(pat_summ), '
<h2>MA plots, main</h2>', img(file.path(out, "ma_main.png")), '
<h2>Time clusters</h2><p>', length(tsig), ' regions significant in the time LRT, k-means on z-scored group means (VST).</p>',
  htab(cl_tab), img(file.path(out, "time_clusters.png")), '
<h2>MA plots, library-size normalization (sensitivity)</h2>', img(file.path(out, "ma_libsize.png")), '
<h2>Size factors</h2>', htab(fread(file.path(out, "size_factors.tsv"))), '
</body></html>')
writeLines(html, file.path(out, "report.html"))
cat("done\n")
