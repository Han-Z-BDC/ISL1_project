#!/usr/bin/env Rscript
# More sensitive KO vs Ctr tests for one histone mark (scripts/19; agreed 2026-11-09).
# Same samples and normalization (median-of-ratios) as histone_deseq2.R. For regions (step 16
# counts) and for 500 bp windows (step 19 counts), pre-filter >= 10 fragments in >= 2 libraries:
#   per day   ~ rep + group, KO vs Ctr at each day (Wald)
#   shared    ~ rep + day + genotype, genotype (Wald): one KO effect estimated from all days
# Each test twice: on padj alone (no fold-change cut), and as a 1.5-fold lfcThreshold test
# (results(lfcThreshold = 0.585, altHypothesis = "greaterAbs"); suffix _fc1.5).
# Windows: p values combined per region with Simes, BH across regions (csaw-style region FDR);
# direction and log2FC from the region's best window.
# Output in <diff_dir>/sensitive/: summary.tsv, regions_all_tests.tsv (every region; empty = NA),
#   significant_any.tsv (regions significant in any test), windows_all_tests.tsv, report.html
# usage: Rscript histone_deseq2_sensitive.R <results/diff/<mark>> <results/diff/<mark>/win500>

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dd <- args[1]; wd <- args[2]
mark <- basename(dd)
out <- file.path(dd, "sensitive")
dir.create(out, showWarnings = FALSE)
PADJ <- 0.05; LFC <- 0.585; MIN_COUNT <- 10; MIN_LIBS <- 2
register(MulticoreParam(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))))

samples <- fread(file.path(dd, "samples.tsv"))
days <- intersect(c("D9", "D15", "D18", "D22"), unique(samples$day))
glev <- c(paste0(days, "_H3CTR"), paste0(days, "_ISL1KO"))
cd <- data.frame(row.names = samples$sample,
                 day = factor(samples$day, levels = days),
                 genotype = factor(samples$genotype, levels = c("H3CTR", "ISL1KO")),
                 rep = factor(samples$replicate),
                 group = factor(samples$group, levels = glev))
ann <- fread(file.path(dd, "deseq2", "tables", paste0("KO_vs_Ctr_", days[length(days)], ".tsv")))[
  , .(region, chr, start, end, n_groups, groups, annotation, SYMBOL, distanceToTSS)]

load_counts <- function(f) {
  fc <- fread(f, skip = "Geneid")
  m <- as.matrix(fc[, 7:ncol(fc)]); rownames(m) <- fc$Geneid
  colnames(m) <- sub("\\.filtered\\.bam$", "", basename(colnames(m)))
  stopifnot(identical(colnames(m), samples$sample))
  m[rowSums(m >= MIN_COUNT) >= MIN_LIBS, , drop = FALSE]
}
gv <- function(rn, g) { v <- setNames(numeric(length(rn)), rn); if (g != glev[1]) v[paste0("group_", g, "_vs_", glev[1])] <- 1; v }
tests <- c(days, "shared")
tnames <- c(tests, paste0(tests, "_fc1.5"))

analyze <- function(m) {
  dg <- DESeq(DESeqDataSetFromMatrix(m, cd, design = ~ rep + group), quiet = TRUE, parallel = TRUE)
  da <- DESeq(DESeqDataSetFromMatrix(m, cd, design = ~ rep + day + genotype), quiet = TRUE, parallel = TRUE)
  rn <- resultsNames(dg)
  res <- data.table(id = rownames(m), baseMean = rowMeans(counts(dg, normalized = TRUE)))
  for (t in tests) {
    get_res <- function(thr) {
      if (t == "shared") results(da, name = "genotype_ISL1KO_vs_H3CTR", alpha = PADJ, lfcThreshold = thr, altHypothesis = "greaterAbs", parallel = TRUE)
      else results(dg, contrast = gv(rn, paste0(t, "_ISL1KO")) - gv(rn, paste0(t, "_H3CTR")), alpha = PADJ,
                   lfcThreshold = thr, altHypothesis = "greaterAbs", parallel = TRUE)
    }
    r0 <- get_res(0); r1 <- get_res(LFC)
    res[, (paste0(c("lfc_", "pvalue_", "padj_"), t)) := list(r0$log2FoldChange, r0$pvalue, r0$padj)]
    res[, (paste0(c("lfc_", "pvalue_", "padj_"), t, "_fc1.5")) := list(r1$log2FoldChange, r1$pvalue, r1$padj)]
  }
  res
}
count_sig <- function(p, l) c(tested = sum(!is.na(p)), up = sum(!is.na(p) & p < PADJ & l > 0), down = sum(!is.na(p) & p < PADJ & l < 0))

# ---- regions
mr <- load_counts(file.path(dd, "counts.txt"))
cat(sprintf("%s regions: %d after the pre-filter\n", mark, nrow(mr)))
rr <- analyze(mr); setnames(rr, "id", "region")
summ <- list()
for (t in tnames) { x <- count_sig(rr[[paste0("padj_", t)]], rr[[paste0("lfc_", t)]])
  summ[[length(summ) + 1]] <- data.table(level = "region", test = t, tested = x[1], up = x[2], down = x[3]) }

# ---- windows
win <- fread(file.path(wd, "windows.bed"), header = FALSE, col.names = c("chr", "start", "end", "window", "region"))
mw <- load_counts(file.path(wd, "counts.txt"))
cat(sprintf("%s %s windows: %d of %d after the pre-filter\n", mark, basename(wd), nrow(mw), nrow(win)))
wr <- analyze(mw); setnames(wr, "id", "window")
wr <- merge(win, wr, by = "window")
fwrite(wr, file.path(out, "windows_all_tests.tsv"), sep = "\t")
simes <- function(p) { p <- sort(p[!is.na(p)]); if (!length(p)) NA_real_ else min(1, min(length(p) * p / seq_along(p))) }
agg <- wr[, .(n_windows_tested = .N), by = region]
for (t in tnames) {
  x <- count_sig(wr[[paste0("padj_", t)]], wr[[paste0("lfc_", t)]])
  summ[[length(summ) + 1]] <- data.table(level = "window", test = t, tested = x[1], up = x[2], down = x[3])
  tmp <- wr[, .(region, p = get(paste0("pvalue_", t)), l = get(paste0("lfc_", t)))]
  a <- tmp[, .(simes_p = simes(p), best_lfc = if (all(is.na(p))) NA_real_ else l[which.min(p)]), by = region]
  a[, padj := p.adjust(simes_p, "BH")]
  setnames(a, c("region", paste0(c("win_simes_p_", "win_lfc_", "win_padj_"), t)))
  agg <- merge(agg, a, by = "region", all.x = TRUE)
  x <- count_sig(agg[[paste0("win_padj_", t)]], agg[[paste0("win_lfc_", t)]])
  summ[[length(summ) + 1]] <- data.table(level = "region_from_windows", test = t, tested = x[1], up = x[2], down = x[3])
}

# ---- per region, all tests
all_ <- merge(ann, rr, by = "region", all.x = TRUE)
all_ <- merge(all_, agg, by = "region", all.x = TRUE)
sig_cols <- c(paste0("padj_", tnames), paste0("win_padj_", tnames))
all_[, n_tests_significant := rowSums(as.matrix(.SD) < PADJ, na.rm = TRUE), .SDcols = sig_cols]
fwrite(all_, file.path(out, "regions_all_tests.tsv"), sep = "\t")
fwrite(all_[n_tests_significant > 0][order(-n_tests_significant)], file.path(out, "significant_any.tsv"), sep = "\t")

# Overlap of region-level and window-level calls, per test (padj only)
ov <- rbindlist(lapply(tests, function(t) {
  a <- !is.na(all_[[paste0("padj_", t)]]) & all_[[paste0("padj_", t)]] < PADJ
  b <- !is.na(all_[[paste0("win_padj_", t)]]) & all_[[paste0("win_padj_", t)]] < PADJ
  data.table(test = t, region_level = sum(a), windows = sum(b), both = sum(a & b), region_only = sum(a & !b), windows_only = sum(b & !a))
}))
summ_dt <- rbindlist(summ)
fwrite(summ_dt, file.path(out, "summary.tsv"), sep = "\t")
fwrite(ov, file.path(out, "region_vs_window_overlap.tsv"), sep = "\t")
print(summ_dt); print(ov)

htab <- function(d) {
  d <- as.data.frame(d); d[is.na(d)] <- ""
  paste0("<table><tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>",
         paste0(apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")), collapse = ""), "</table>")
}
writeLines(paste0('<!doctype html><html><head><meta charset="utf-8"><title>', mark, ' KO vs Ctr, sensitive</title>
<style>body{font-family:sans-serif;max-width:1100px;margin:24px auto;padding:0 16px}table{border-collapse:collapse;font-size:12px;margin-bottom:16px}
td,th{border:1px solid #ccc;padding:3px 6px;text-align:right}</style></head><body>
<h1>', mark, ': KO vs Ctr, more sensitive tests</h1>
<p>Significant: padj &lt; ', PADJ, '. Tests: KO vs Ctr per day (~ rep + group) and shared across days (~ rep + day + genotype);
_fc1.5 = DESeq2 lfcThreshold test of |log2FC| &gt; ', LFC, '. Levels: region (', nrow(mr), ' regions), window (', nrow(mw),
' windows of ', basename(wd), '), region_from_windows (window p values combined per region with Simes, BH across regions).
up = higher in KO.</p><h2>Counts</h2>', htab(summ_dt), '
<h2>Region-level vs window-level calls (padj only)</h2>', htab(ov), '
</body></html>'), file.path(out, "report.html"))
cat("done\n")
