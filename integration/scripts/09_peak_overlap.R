#!/usr/bin/env Rscript
# Step 09: overlap of the per-day ISL1 peak calls (D15 / D18 / D22), checked against the count-based tests.
#   Presence   each of the 14,431 sites (step 26/27 sites after the IgG greylist) carries the days on which it lies in
#              that day's consensus peaks (macs3_ge2: MACS3 q 0.01 peaks in >= 2 of the day's 3 libraries, step 18).
#              The 7 combinations of days are the overlap classes.
#   Libraries  each site is also overlapped with each library's own MACS3 peaks: number of a day's libraries that call
#              it (0-3), and whether the batch 3 library does. Per library: peaks, sites called, library size and
#              share of reads in the 14,431 sites (DiffBind counts, step 00), to show how call rates follow the library.
#   Counts     per class: batch-corrected VST per day and the step 08 LRT / binding-pattern group; DiffBind
#              (background normalization, step 27 via step 03) log2FC and FDR per contrast.
#   Check      per contrast A -> B: sites called on B but not A, on A but not B, on both, on neither; for each, the
#              median DiffBind log2FC and the share with FDR < 0.05 up / down. If presence differences reflect binding
#              changes, the sites that appear should mostly go up in the counts.
# usage: Rscript 09_peak_overlap.R <step 27 isl1_diff dir> <integration results dir> <out_dir>
# Output in <out_dir>: sites_days.tsv (per site), classes.tsv (per class), class_days.tsv (per class and day),
#   presence_vs_counts.tsv (per contrast and presence category), libraries.tsv (per library)

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
cnr <- args[1]; res <- args[2]; od <- args[3]
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
smp[, peaks := ifelse(batch == "batch3",
                      file.path(root, "macs3", paste0(original_sample, "_peaks.filtered.narrowPeak")),
                      file.path(root, "ctr_matched", "macs3", paste0(sample, "_peaks.filtered.narrowPeak")))]
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

nf <- fread(file.path(res, "00_diffbind_counts", "norm_factors.tsv"))
cr <- fread(file.path(res, "00_diffbind_counts", "counts_raw.tsv"))
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
