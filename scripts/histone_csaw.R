#!/usr/bin/env Rscript
# csaw cross-check for one histone mark (scripts/22; agreed 2026-11-09), following the csawBook
# workflows (window counts genome-wide, no peak calls), edgeR quasi-likelihood.
# Same 12 libraries as steps 16-17; paired-end fragments up to 1000 bp; chr1-22/X/Y; blacklist
# discarded; duplicates kept.
#   windows     WIDTHS (default 150 bp every 50 bp and 1000 bp every 250 bp; set CSAW_WIDTHS e.g. "2000,5000"
#               with spacing a quarter of the width for broad marks)
#   filter      filterWindowsGlobal against 10 kb bins, keep windows >= FC_BG x background
#   normalization  composition (TMM on 10 kb bins, csaw's choice for global shifts) and efficiency
#               (TMM on the kept windows)
#   tests       ~ 0 + group + rep: KO vs Ctr at each day, time within each genotype;
#               ~ rep + day + genotype: shared KO vs Ctr (genotype coefficient)
#   regions     the window sizes merged with mergeResultsList (tol 100, max.width MAX_WIDTH,
#               equiweight), region FDR from the combined test (csaw)
# Compared with our DESeq2 region calls (padj only) as in histone_diffbind.R.
# Output in <diff_dir>/csaw/: summary.tsv, comparison.tsv, regions/<norm>_<test>.tsv, report.html
# usage: Rscript histone_csaw.R <results/diff/<mark>> <blacklist.bed>

suppressPackageStartupMessages({
  library(csaw)
  library(edgeR)
  library(data.table)
  library(GenomicRanges)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dd <- args[1]; bl_file <- args[2]
mark <- basename(dd)
out <- file.path(dd, "csaw_edgeR"); dir.create(file.path(out, "regions"), recursive = TRUE, showWarnings = FALSE)
TH <- 0.05
widths <- as.integer(strsplit(Sys.getenv("CSAW_WIDTHS", "150,1000"), ",")[[1]])
spacing <- ifelse(widths <= 150, 50L, as.integer(widths / 4))
FC_BG <- as.numeric(Sys.getenv("CSAW_FC_BG", "3"))
MAX_WIDTH <- as.integer(Sys.getenv("CSAW_MAX_WIDTH", "5000"))
bp <- MulticoreParam(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

samples <- fread(file.path(dd, "samples.tsv"))
days <- intersect(c("D9", "D15", "D18", "D22"), unique(samples$day))
glev <- c(paste0(days, "_H3CTR"), paste0(days, "_ISL1KO"))
cd <- data.frame(group = factor(samples$group, levels = glev), rep = factor(samples$replicate),
                 day = factor(samples$day, levels = days), genotype = factor(samples$genotype, levels = c("H3CTR", "ISL1KO")))
bams <- samples$bam
bl <- fread(bl_file, header = FALSE)
blgr <- GRanges(bl$V1, IRanges(bl$V2 + 1, bl$V3))
chroms <- paste0("chr", c(1:22, "X", "Y"))
param <- readParam(pe = "both", max.frag = 1000, restrict = chroms, discard = blgr, dedup = FALSE, minq = 0)

bins <- windowCounts(bams, bin = TRUE, width = 10000, param = param, BPPARAM = bp)
cat(sprintf("10 kb bins: %d\n", nrow(bins)))
filt <- list()
for (i in seq_along(widths)) {
  w <- widths[i]
  win <- windowCounts(bams, width = w, spacing = spacing[i], param = param, BPPARAM = bp)
  fs <- filterWindowsGlobal(win, bins)
  keep <- fs$filter > log2(FC_BG)
  filt[[as.character(w)]] <- win[keep, ]
  cat(sprintf("%d bp windows (spacing %d): %d counted, %d kept (>= %gx background)\n", w, spacing[i], nrow(win), sum(keep), FC_BG))
  rm(win); gc()
}

# contrasts on ~ 0 + group + rep
X1 <- model.matrix(~ 0 + group + rep, cd); colnames(X1) <- sub("^group", "", colnames(X1))
con <- list()
for (d in days) con[[paste0("KO_vs_Ctr_", d)]] <- c(paste0(d, "_ISL1KO"), paste0(d, "_H3CTR"))
for (g in c("H3CTR", "ISL1KO")) for (ab in list(c("D18", "D15"), c("D22", "D18"), c("D22", "D15")))
  if (all(ab %in% days)) con[[sprintf("time_%s_%s_vs_%s", g, ab[1], ab[2])]] <- paste0(ab, "_", g)
cvec <- function(p) { v <- setNames(numeric(ncol(X1)), colnames(X1)); v[p[1]] <- 1; v[p[2]] <- -1; v }
X2 <- model.matrix(~ rep + day + genotype, cd)
tests <- c(names(con), "KO_vs_Ctr_shared")

fit_size <- function(se) {
  y <- asDGEList(se)
  f1 <- glmQLFit(estimateDisp(y, X1), X1, robust = TRUE)
  f2 <- glmQLFit(estimateDisp(y, X2), X2, robust = TRUE)
  res <- lapply(names(con), function(t) glmQLFTest(f1, contrast = cvec(con[[t]]))$table)
  names(res) <- names(con)
  res$KO_vs_Ctr_shared <- glmQLFTest(f2, coef = "genotypeISL1KO")$table
  res
}

summ <- list(); sig_regions <- list()
for (nm in c("composition", "efficiency")) {
  se <- lapply(filt, function(s) if (nm == "composition") normFactors(bins, se.out = s) else normFactors(s, se.out = s))
  cat(sprintf("%s norm factors: %s\n", nm, paste(round(se[[1]]$norm.factors, 3), collapse = " ")))
  R <- lapply(se, fit_size)
  for (t in tests) {
    tabs <- lapply(R, `[[`, t)
    # unnamed lists: with names, csaw's internal do.call(c, ...) falls back to base c() and returns a list
    m <- mergeResultsList(unname(lapply(se, rowRanges)), tab.list = unname(tabs), equiweight = TRUE, tol = 100,
                          merge.args = list(max.width = MAX_WIDTH))
    reg <- m$regions; comb <- m$combined
    dt <- data.table(chr = as.character(seqnames(reg)), start = start(reg) - 1L, end = end(reg),
                     n_windows = comb$num.tests, PValue = comb$PValue, FDR = comb$FDR, direction = comb$direction,
                     best_logFC = if ("rep.logFC" %in% names(comb)) comb$rep.logFC else NA_real_)
    fwrite(dt[order(FDR)], file.path(out, "regions", paste0(nm, "_", t, ".tsv")), sep = "\t")
    s <- dt[FDR < TH]
    nwin <- sapply(names(tabs), function(k) sum(p.adjust(tabs[[k]]$PValue, "BH") < TH))
    summ[[length(summ) + 1]] <- data.table(normalization = nm, test = t, regions_tested = nrow(dt),
      up = sum(s$direction == "up"), down = sum(s$direction == "down"), mixed = sum(s$direction == "mixed"),
      windows_sig = paste(paste0(names(nwin), "bp:", nwin), collapse = " "))
    sig_regions[[paste(nm, t)]] <- GRanges(s$chr, IRanges(s$start + 1, s$end))
    cat(sprintf("%s %s: %d up, %d down, %d mixed regions (FDR < %.2f)\n", nm, t, sum(s$direction == "up"), sum(s$direction == "down"), sum(s$direction == "mixed"), TH))
  }
}
summ_dt <- rbindlist(summ)
fwrite(summ_dt, file.path(out, "summary.tsv"), sep = "\t")

# comparison with our DESeq2 region calls (padj only)
regions <- fread(file.path(dd, "regions.bed"), header = FALSE, col.names = c("chr", "start", "end", "region", "n_groups", "groups"))
reg_gr <- GRanges(regions$chr, IRanges(regions$start + 1, regions$end), region = regions$region)
st <- fread(file.path(dd, "sensitive", "regions_all_tests.tsv"))
ours <- list()
for (t in tests) {
  if (grepl("^KO_vs_Ctr_", t)) { k <- sub("KO_vs_Ctr_", "", t); ours[[t]] <- st[!is.na(get(paste0("padj_", k))) & get(paste0("padj_", k)) < TH, region] }
  else { f <- file.path(dd, "deseq2", "tables", paste0(t, ".tsv")); ours[[t]] <- if (file.exists(f)) fread(f)[!is.na(padj) & padj < TH, region] else character(0) }
}
comp <- rbindlist(lapply(names(sig_regions), function(k) {
  p <- strsplit(k, " ")[[1]]; o <- reg_gr[reg_gr$region %in% ours[[p[2]]]]; s <- sig_regions[[k]]
  data.table(normalization = p[1], test = p[2], csaw_sig = length(s), ours_sig = length(o),
             csaw_in_ours_regions = sum(overlapsAny(s, o)), csaw_in_any_union_region = sum(overlapsAny(s, reg_gr)),
             ours_in_csaw = sum(overlapsAny(o, s)))
}))
fwrite(comp, file.path(out, "comparison.tsv"), sep = "\t")
print(summ_dt); print(comp)

htab <- function(d) {
  d <- as.data.frame(d); d[is.na(d)] <- ""
  paste0("<table><tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>",
         paste0(apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")), collapse = ""), "</table>")
}
writeLines(paste0('<!doctype html><html><head><meta charset="utf-8"><title>', mark, ' csaw</title>
<style>body{font-family:sans-serif;max-width:1200px;margin:24px auto;padding:0 16px}table{border-collapse:collapse;font-size:12px;margin-bottom:16px}
td,th{border:1px solid #ccc;padding:3px 6px;text-align:right}</style></head><body>
<h1>', mark, ': csaw (', as.character(packageVersion("csaw")), ') + edgeR QL</h1>
<p>Genome-wide windows ', paste(widths, collapse = " and "), ' bp, kept at &ge; ', FC_BG, 'x the 10 kb-bin background; sizes merged per region
(tol 100, max width ', MAX_WIDTH, ' bp, equal weight); region FDR &lt; ', TH, '. up = higher in KO, or in the later day.
Normalization: composition (TMM on 10 kb bins) or efficiency (TMM on kept windows).</p>
<h2>Counts</h2>', htab(summ_dt), '
<h2>Overlap with our DESeq2 region calls (padj &lt; ', TH, ')</h2>', htab(comp), '
</body></html>'), file.path(out, "report.html"))
cat("done\n")
