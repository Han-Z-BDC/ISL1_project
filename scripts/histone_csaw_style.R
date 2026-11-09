#!/usr/bin/env Rscript
# csaw-style differential analysis for one histone mark (scripts/20; agreed 2026-11-09), after the
# csawBook H3K27me3 WT vs Ezh2 KO chapter, with DESeq2 as the engine.
#   windows    2 kb and 5 kb, tiled over the step 16 domains
#   filter     >= 10 fragments in >= 2 libraries, and mean CPM per bp >= 2x the background level
#              (median over the 10 kb genome-wide bins of the mean CPM per bp), as filterWindowsGlobal
#   normalization, three versions compared:
#     background        DESeq2 median-of-ratios on the 10 kb genome-wide bins (csaw: TMM on 10 kb bins)
#     libsize           all fragments featureCounts saw per library (geometric mean 1)
#     median_of_ratios  DESeq2 median-of-ratios on the step 16 domain counts (the original method)
#   tests      ~ rep + group: KO vs Ctr at each day, time within each genotype;
#              ~ rep + day + genotype: shared KO vs Ctr. Called on padj alone (< PADJ).
#   domains    per domain and window size, Simes over its windows; the two sizes combined with
#              Simes (equal weight per size, as csaw's equiweight); BH across domains; direction
#              and log2FC from the domain's best window
#   sub-domains  2 kb windows with window padj < PADJ, joined when adjacent, same domain and same
#              direction, capped at 30 kb (csaw max.width)
# Compared with the original method: whole-domain DESeq2 (step 17 tables, padj alone) and 5 kb
# windows with median-of-ratios (step 19, Simes per domain).
# Output in <csaw_dir>/: summary.tsv, comparison.tsv, domains_all_tests.tsv, subdomains/*.bed,
#   size_factors.tsv, report.html
# usage: Rscript histone_csaw_style.R <results/diff/<mark>> <results/diff/<mark>/csaw>

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
dd <- args[1]; cdir <- args[2]
mark <- basename(dd)
PADJ <- 0.05; MAXW <- 30000; MIN_COUNT <- 10; MIN_LIBS <- 2; FC_BG <- 2
register(MulticoreParam(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))))
dir.create(file.path(cdir, "subdomains"), showWarnings = FALSE)

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

load_fc <- function(f) {
  x <- fread(f, skip = "Geneid")
  m <- as.matrix(x[, 7:ncol(x)]); rownames(m) <- x$Geneid
  colnames(m) <- sub("\\.filtered\\.bam$", "", basename(colnames(m)))
  stopifnot(identical(colnames(m), samples$sample))
  m
}
summ0 <- fread(file.path(dd, "counts.txt.summary"))
total <- colSums(as.matrix(summ0[, -1])); names(total) <- sub("\\.filtered\\.bam$", "", basename(names(total)))
total <- total[samples$sample]

# ---- size factors
bins <- load_fc(file.path(cdir, "bins10k.counts.txt"))
dom <- load_fc(file.path(dd, "counts.txt"))
dom <- dom[rowSums(dom >= MIN_COUNT) >= MIN_LIBS, ]
SF <- list(background = estimateSizeFactorsForMatrix(bins[apply(bins, 1, min) > 0, ]),
           libsize = total / exp(mean(log(total))),
           median_of_ratios = estimateSizeFactorsForMatrix(dom))
sf_tab <- data.table(sample = samples$sample, group = samples$group, fragments = total,
                     sf_background = SF$background, sf_libsize = SF$libsize, sf_median_of_ratios = SF$median_of_ratios)
fwrite(sf_tab, file.path(cdir, "size_factors.tsv"), sep = "\t")
print(sf_tab)
bin_cpm <- rowMeans(sweep(bins, 2, total / 1e6, "/"))
bg_bp <- median(bin_cpm[bin_cpm > 0]) / 1e4
cat(sprintf("background: median mean CPM per 10 kb bin %.4f (%.3g per bp)\n", median(bin_cpm[bin_cpm > 0]), bg_bp))

# ---- windows
W <- list()
for (k in c("2k", "5k")) {
  bed <- fread(file.path(cdir, paste0("windows_", k, ".bed")), header = FALSE, col.names = c("chr", "start", "end", "window", "region"))
  m <- load_fc(file.path(cdir, paste0("windows_", k, ".counts.txt")))
  width <- (bed$end - bed$start)[match(rownames(m), bed$window)]
  cpm_bp <- rowMeans(sweep(m, 2, total / 1e6, "/")) / width
  keep <- rowSums(m >= MIN_COUNT) >= MIN_LIBS & cpm_bp >= FC_BG * bg_bp
  cat(sprintf("%s windows: %d, %d pass the filter (>= %d in >= %d libraries and >= %gx background)\n",
              k, nrow(m), sum(keep), MIN_COUNT, MIN_LIBS, FC_BG))
  W[[k]] <- list(bed = bed, m = m[keep, ])
}

# ---- tests
gv <- function(rn, g) { v <- setNames(numeric(length(rn)), rn); if (g != glev[1]) v[paste0("group_", g, "_vs_", glev[1])] <- 1; v }
pairs_ <- list()
for (d in days) pairs_[[paste0("KO_vs_Ctr_", d)]] <- c(paste0(d, "_ISL1KO"), paste0(d, "_H3CTR"))
for (g in c("H3CTR", "ISL1KO")) for (ab in list(c("D18", "D15"), c("D22", "D18"), c("D22", "D15")))
  if (all(ab %in% days)) pairs_[[sprintf("time_%s_%s_vs_%s", g, ab[1], ab[2])]] <- paste0(ab, "_", g)
tests <- c(paste0("KO_vs_Ctr_", days), "KO_vs_Ctr_shared", grep("^time_", names(pairs_), value = TRUE))

fit_all <- function(m, sf) {
  dg <- DESeqDataSetFromMatrix(m, cd, design = ~ rep + group); sizeFactors(dg) <- sf[colnames(m)]
  dg <- DESeq(dg, quiet = TRUE, parallel = TRUE)
  da <- DESeqDataSetFromMatrix(m, cd, design = ~ rep + day + genotype); sizeFactors(da) <- sf[colnames(m)]
  da <- DESeq(da, quiet = TRUE, parallel = TRUE)
  rn <- resultsNames(dg)
  out <- data.table(window = rownames(m))
  for (t in tests) {
    r <- if (t == "KO_vs_Ctr_shared") results(da, name = "genotype_ISL1KO_vs_H3CTR", alpha = PADJ, parallel = TRUE)
         else results(dg, contrast = gv(rn, pairs_[[t]][1]) - gv(rn, pairs_[[t]][2]), alpha = PADJ, parallel = TRUE)
    out[, (paste0(c("lfc_", "pvalue_", "padj_"), t)) := list(r$log2FoldChange, r$pvalue, r$padj)]
  }
  out
}
simes <- function(p) { p <- sort(p[!is.na(p)]); if (!length(p)) NA_real_ else min(1, min(length(p) * p / seq_along(p))) }
R <- list()
for (nm in names(SF)) for (k in names(W)) {
  cat(sprintf("fitting %s, %s windows\n", nm, k))
  R[[nm]][[k]] <- merge(W[[k]]$bed, fit_all(W[[k]]$m, SF[[nm]]), by = "window")
}

# ---- domains and sub-domains
summ <- list(); dom_all <- copy(ann)
for (nm in names(SF)) for (t in tests) {
  wl <- rbindlist(lapply(names(W), function(k) R[[nm]][[k]][, .(region, size = k, p = get(paste0("pvalue_", t)), l = get(paste0("lfc_", t)))]))
  ps <- wl[, .(ps = simes(p), lbest = if (all(is.na(p))) NA_real_ else l[which.min(p)],
               pmin = if (all(is.na(p))) Inf else min(p, na.rm = TRUE)), by = .(region, size)]
  rg <- ps[, .(p = simes(ps), lfc = lbest[which.min(pmin)]), by = region]
  rg[, padj := p.adjust(p, "BH")]
  setnames(rg, c("region", paste0(c("p_", "lfc_", "padj_"), nm, "_", t)))
  dom_all <- merge(dom_all, rg, by = "region", all.x = TRUE)
  dp <- rg[[paste0("padj_", nm, "_", t)]]; dl <- rg[[paste0("lfc_", nm, "_", t)]]
  s <- R[[nm]][["2k"]][!is.na(get(paste0("padj_", t))) & get(paste0("padj_", t)) < PADJ,
                       .(chr, start, end, region, l = get(paste0("lfc_", t)))][order(chr, start)]
  sub <- data.table()
  if (nrow(s)) {
    s[, dir := sign(l)]
    g <- integer(nrow(s)); g[1] <- 1L; st <- s$start[1]
    if (nrow(s) > 1) for (i in 2:nrow(s)) {
      same <- s$chr[i] == s$chr[i - 1] && s$start[i] <= s$end[i - 1] && s$region[i] == s$region[i - 1] &&
              s$dir[i] == s$dir[i - 1] && s$end[i] - st <= MAXW
      if (same) g[i] <- g[i - 1] else { g[i] <- g[i - 1] + 1L; st <- s$start[i] }
    }
    s[, sd := g]
    sub <- s[, .(chr = chr[1], start = min(start), end = max(end), region = region[1], n_windows = .N, mean_lfc = round(mean(l), 3),
                 direction = ifelse(dir[1] > 0, "up", "down")), by = sd][, sd := NULL]
    fwrite(sub, file.path(cdir, "subdomains", paste0(nm, "_", t, ".bed")), sep = "\t", col.names = FALSE)
  }
  w2 <- R[[nm]][["2k"]]; w5 <- R[[nm]][["5k"]]
  cnt <- function(p, l) c(sum(!is.na(p) & p < PADJ & l > 0), sum(!is.na(p) & p < PADJ & l < 0))
  a <- cnt(dp, dl); b <- cnt(w2[[paste0("padj_", t)]], w2[[paste0("lfc_", t)]]); c5 <- cnt(w5[[paste0("padj_", t)]], w5[[paste0("lfc_", t)]])
  summ[[length(summ) + 1]] <- data.table(normalization = nm, test = t, domains_up = a[1], domains_down = a[2],
    win2k_up = b[1], win2k_down = b[2], win5k_up = c5[1], win5k_down = c5[2],
    subdomains_up = if (nrow(sub)) sum(sub$direction == "up") else 0L, subdomains_down = if (nrow(sub)) sum(sub$direction == "down") else 0L,
    subdomain_Mb = if (nrow(sub)) round(sum(sub$end - sub$start) / 1e6, 2) else 0)
}
summ_dt <- rbindlist(summ)
fwrite(summ_dt, file.path(cdir, "summary.tsv"), sep = "\t")
fwrite(dom_all, file.path(cdir, "domains_all_tests.tsv"), sep = "\t")
print(summ_dt)

# ---- comparison with the original method (whole domains, step 17; 5 kb windows, step 19)
sens <- file.path(dd, "sensitive", "regions_all_tests.tsv")
st <- if (file.exists(sens)) fread(sens) else NULL
comp <- list()
sig <- function(x) !is.na(x) & x < PADJ
for (t in tests) {
  tab <- file.path(dd, "deseq2", "tables", paste0(if (t == "KO_vs_Ctr_shared") "none" else t, ".tsv"))
  orig_dom <- if (t == "KO_vs_Ctr_shared" && !is.null(st)) st[sig(padj_shared), region]
              else if (file.exists(tab)) fread(tab)[sig(padj), region] else character(0)
  key <- sub("KO_vs_Ctr_", "", t)
  orig_win <- if (!is.null(st) && paste0("win_padj_", key) %in% names(st)) st[sig(get(paste0("win_padj_", key))), region] else character(0)
  bg <- dom_all[sig(get(paste0("padj_background_", t))), region]
  lib <- dom_all[sig(get(paste0("padj_libsize_", t))), region]
  mor <- dom_all[sig(get(paste0("padj_median_of_ratios_", t))), region]
  comp[[length(comp) + 1]] <- data.table(test = t, original_domain = length(orig_dom), original_win5k = length(orig_win),
    csaw_background = length(bg), csaw_libsize = length(lib), csaw_median_of_ratios = length(mor),
    background_and_original_domain = length(intersect(bg, orig_dom)), background_only = length(setdiff(bg, union(orig_dom, orig_win))),
    original_domain_only = length(setdiff(orig_dom, union(bg, union(lib, mor)))))
}
comp_dt <- rbindlist(comp)
fwrite(comp_dt, file.path(cdir, "comparison.tsv"), sep = "\t")
print(comp_dt)

htab <- function(d) {
  d <- as.data.frame(d); d[is.na(d)] <- ""
  paste0("<table><tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>",
         paste0(apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")), collapse = ""), "</table>")
}
writeLines(paste0('<!doctype html><html><head><meta charset="utf-8"><title>', mark, ' csaw-style</title>
<style>body{font-family:sans-serif;max-width:1300px;margin:24px auto;padding:0 16px}table{border-collapse:collapse;font-size:12px;margin-bottom:16px}
td,th{border:1px solid #ccc;padding:3px 6px;text-align:right}</style></head><body>
<h1>', mark, ': csaw-style differential analysis vs the original method</h1>
<p>Windows of 2 kb (', nrow(W[["2k"]]$m), ' kept) and 5 kb (', nrow(W[["5k"]]$m), ' kept) over the domains; kept if &ge; ', MIN_COUNT,
' fragments in &ge; ', MIN_LIBS, ' libraries and &ge; ', FC_BG, 'x the genome-wide background. Domain calls: Simes per window size, the two sizes
combined with Simes, BH across domains (padj &lt; ', PADJ, '). Sub-domains: adjacent significant 2 kb windows, same direction, &le; 30 kb.
up = higher in KO (KO vs Ctr) or in the later day (time).</p>
<h2>csaw-style counts</h2>', htab(summ_dt), '
<h2>Comparison with the original method (domains significant, padj &lt; ', PADJ, ')</h2>
<p>original_domain: whole-domain DESeq2 (step 17; shared test from step 19). original_win5k: 5 kb windows with median-of-ratios, Simes per domain
(step 19). background_only: csaw background calls missed by both original versions. original_domain_only: original domain calls missed by all three
csaw normalizations.</p>', htab(comp_dt), '
<h2>Size factors</h2>', htab(sf_tab), '
</body></html>'), file.path(cdir, "report.html"))
cat("done\n")
