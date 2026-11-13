#!/usr/bin/env Rscript
# Step 08: is the differential binding more robust with only the 6 earlier (non-batch-3) libraries?
#   Libraries: step 27's trio minus the batch 3 library of each day (D15 rep3/4, D18 rep3/4, D22 rep4/5; D18 rep4,
#   D15_CtrGFP_rep2, came with the batch 3 delivery but is labelled "earlier" in step 27, as agreed).
#   Same sites, counts and normalization as step 27 (DiffBind object, DESeq2 native, background 15 kb bins, size
#   factors from the 6 libraries); design ~ day (no batch term); contrasts D18 vs D15, D22 vs D18, D22 vs D15;
#   FDR < 0.05.
#   Held-out check: the 3 batch 3 libraries (D15 rep1, D18 rep1, D22 rep2) form an independent time course. For each
#   contrast, log2((b + 1) / (a + 1)) of their normalized signal (step 00, 9-library size factors) per site; share
#   of significant sites whose batch 3 change has the same sign as the 6-library log2FC, against all sites.
#   Alpha vs beta per day (as step 06, median signal per site): log2(alpha / beta marker sites) per library,
#   lm(~ day) on the 6 libraries (3 residual df).
# usage: Rscript 08_six_library.R <step 27 isl1_diff dir> <integration results dir> <markers.tsv> <out_dir>
# Output in <out_dir>: tables/diffbind_<contrast>.tsv, norm_factors.tsv, summary.tsv, overlap.tsv, heldout.tsv,
#   alpha_beta_tests.tsv

suppressPackageStartupMessages({ library(DiffBind); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
d <- args[1]; res <- args[2]; mfile <- args[3]; od <- args[4]
dir.create(file.path(od, "tables"), showWarnings = FALSE, recursive = TRUE)
TH <- 0.05
contrasts <- list(D18_vs_D15 = c("D18", "D15"), D22_vs_D18 = c("D22", "D18"), D22_vs_D15 = c("D22", "D15"))

smp <- fread(file.path(d, "samples.tsv"))
sites <- fread(file.path(d, "sites_used.bed"), header = FALSE,
               col.names = c("chr", "start", "end", "name", "score", "strand", "summit_day", "days_called", "n_days"))
name_of <- setNames(sites$name, paste(sites$chr, sites$start + 1, sep = ":"))
db <- readRDS(file.path(d, "diffbind.counted.rds"))
keep <- db$samples$SampleID %in% smp[batch == "earlier", sample]
stopifnot(sum(keep) == 6)
db6 <- dba(db, mask = keep)
db6$config$cores <- as.integer(Sys.getenv("DB_CORES", "4"))
cat("libraries:", paste(dba.show(db6)$ID, collapse = ", "), "\n")

dbn <- dba.normalize(db6, method = DBA_DESEQ2, normalize = DBA_NORM_NATIVE, library = DBA_LIBSIZE_BACKGROUND,
                     background = TRUE)
f <- dba.normalize(dbn, method = DBA_DESEQ2, bRetrieve = TRUE)
ids6 <- dba.show(dbn)$ID   # the masked libraries (dbn$samples keeps all 9 rows)
stopifnot(length(ids6) == 6, length(f$norm.factors) == 6, setequal(ids6, smp[batch == "earlier", sample]))
nf6 <- data.table(sample = ids6, lib_size = f$lib.sizes, norm_factor = f$norm.factors)
fwrite(nf6, file.path(od, "norm_factors.tsv"), sep = "\t", quote = FALSE)
dc <- NULL
for (ct in names(contrasts)) {
  g <- contrasts[[ct]]
  dc <- if (is.null(dc)) dba.contrast(dbn, design = "~ Condition", contrast = c("Condition", g[1], g[2]))
        else dba.contrast(dc, contrast = c("Condition", g[1], g[2]))
}
dc <- dba.analyze(dc, method = DBA_DESEQ2, bBlacklist = FALSE, bGreylist = FALSE)
print(dba.show(dc, bContrasts = TRUE))   # group sizes per contrast (2 vs 2)

tab6 <- list(); tab9 <- list()
for (i in seq_along(contrasts)) {
  ct <- names(contrasts)[i]
  r <- as.data.table(as.data.frame(dba.report(dc, contrast = i, th = 1, method = DBA_DESEQ2, bCounts = FALSE)))
  r[, name := name_of[paste(seqnames, start, sep = ":")]]
  stopifnot(!anyNA(r$name), nrow(r) == nrow(sites))
  setnames(r, c("Fold", "p.value"), c("log2FC", "pvalue"))
  r <- r[, .(name, seqnames, start, end, Conc, log2FC, pvalue, FDR)]
  setorder(r, FDR, pvalue, name)
  fwrite(r, file.path(od, "tables", sprintf("diffbind_%s.tsv", ct)), sep = "\t", quote = FALSE)
  tab6[[ct]] <- r
  tab9[[ct]] <- fread(file.path(d, "tables", sprintf("diffbind_background_%s.tsv", ct)))[, .(name, log2FC, FDR)]
}

summ <- rbindlist(lapply(names(contrasts), function(ct) rbindlist(lapply(c("6 earlier libraries", "9 libraries (step 27)"),
  function(set) {
    r <- if (set == "6 earlier libraries") tab6[[ct]] else tab9[[ct]]
    s <- r[FDR < TH]
    data.table(contrast = ct, libraries = set, significant = nrow(s), up = sum(s$log2FC > 0),
               down = sum(s$log2FC < 0), up_fc1 = sum(s$log2FC >= 1), down_fc1 = sum(s$log2FC <= -1))
  }))))
fwrite(summ, file.path(od, "summary.tsv"), sep = "\t", quote = FALSE)

ov <- rbindlist(lapply(names(contrasts), function(ct) {
  m <- merge(tab6[[ct]][, .(name, lfc6 = log2FC, fdr6 = FDR)], tab9[[ct]][, .(name, lfc9 = log2FC, fdr9 = FDR)], by = "name")
  s6 <- m$fdr6 < TH; s9 <- m$fdr9 < TH
  data.table(contrast = ct, sig_6 = sum(s6), sig_9 = sum(s9), sig_both = sum(s6 & s9),
             both_same_direction = sum(s6 & s9 & sign(m$lfc6) == sign(m$lfc9)), only_6 = sum(s6 & !s9),
             only_9 = sum(!s6 & s9), log2FC_pearson_r = round(cor(m$lfc6, m$lfc9), 3))
}))
fwrite(ov, file.path(od, "overlap.tsv"), sep = "\t", quote = FALSE)

# ---- held-out batch 3 time course (normalized with the 9-library size factors, step 00)
raw <- fread(file.path(res, "00_diffbind_counts", "counts_raw.tsv"))
nf9 <- fread(file.path(res, "00_diffbind_counts", "norm_factors.tsv"))
b3 <- smp[batch == "batch3"]
sig9 <- sweep(as.matrix(raw[, b3$sample, with = FALSE]), 2, nf9$norm_factor[match(b3$sample, nf9$sample)], "/")
rownames(sig9) <- raw$name; colnames(sig9) <- b3$day
ho <- rbindlist(lapply(names(contrasts), function(ct) {
  g <- contrasts[[ct]]
  b3lfc <- log2((sig9[, g[1]] + 1) / (sig9[, g[2]] + 1))
  rbindlist(lapply(list(list("6-library significant", tab6[[ct]][FDR < TH]), list("6-library, all sites", tab6[[ct]]),
                        list("9-library significant (includes batch 3)", tab9[[ct]][FDR < TH])), function(z) {
    r <- z[[2]]
    if (nrow(r) == 0) return(data.table(contrast = ct, set = z[[1]], n_sites = 0L))
    x <- b3lfc[r$name]
    data.table(contrast = ct, set = z[[1]], n_sites = nrow(r),
               pct_same_sign_in_batch3 = round(100 * mean(sign(x) == sign(r$log2FC)), 1),
               median_batch3_log2_ratio = round(median(x), 3),
               pearson_r_with_batch3 = if (nrow(r) > 2) round(cor(r$log2FC, x), 3) else NA_real_)
  }), fill = TRUE)
}), fill = TRUE)
fwrite(ho, file.path(od, "heldout.tsv"), sep = "\t", quote = FALSE)

# ---- alpha vs beta marker sites per day, 6 libraries
sgl <- fread(file.path(res, "01_great", "site_gene.tsv"), sep = "\t", quote = "")
mk <- fread(mfile, skip = "cell_type")
sets <- list(alpha = unique(sgl[gene_name %in% mk[cell_type == "alpha", gene], site]),
             beta = unique(sgl[gene_name %in% mk[cell_type == "beta", gene], site]))
s6 <- smp[batch == "earlier"]
m6 <- sweep(as.matrix(raw[, s6$sample, with = FALSE]), 2, nf6$norm_factor[match(s6$sample, nf6$sample)], "/")
rownames(m6) <- raw$name
ab <- rbindlist(lapply(c("median", "total"), function(metric) {
  summ_ <- function(s) if (metric == "total") colSums(m6[s, ]) else apply(m6[s, ], 2, median)
  y <- log2(summ_(sets$alpha) / summ_(sets$beta))
  dd <- data.frame(y = y, day = factor(s6$day, levels = c("D15", "D18", "D22")))
  fit <- lm(y ~ day, data = dd)
  W <- list(level_D15 = c(1, 0, 0), level_D18 = c(1, 1, 0), level_D22 = c(1, 0, 1),
            change_D18_vs_D15 = c(0, 1, 0), change_D22_vs_D18 = c(0, -1, 1), change_D22_vs_D15 = c(0, 0, 1))
  rbindlist(lapply(names(W), function(tn) {
    w <- W[[tn]]; est <- sum(w * coef(fit)); se <- sqrt(drop(t(w) %*% vcov(fit) %*% w))
    data.table(metric = metric, term = tn, estimate = round(est, 4), se = round(se, 4), df = fit$df.residual,
               p = signif(2 * pt(-abs(est / se), fit$df.residual), 3))
  }))
}))
ab[, `:=`(n_alpha_sites = length(sets$alpha), n_beta_sites = length(sets$beta))]
fwrite(ab, file.path(od, "alpha_beta_tests.tsv"), sep = "\t", quote = FALSE)

print(nf6); print(summ); print(ov); print(ho); print(ab)
