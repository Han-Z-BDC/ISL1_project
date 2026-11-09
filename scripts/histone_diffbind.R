#!/usr/bin/env Rscript
# DiffBind cross-check for one histone mark (scripts/21; agreed 2026-11-09). Same 12 libraries as
# steps 16-17 (D15/D18/D22 x H3CTR/ISL1KO x rep1/rep2), no IgG in DiffBind (IgG was used for the
# peak calls), duplicates kept, paired-end fragments counted (summarizeOverlaps).
# Region sets:
#   diffbind_consensus  DiffBind default: per-library MACS3 peaks (step 03), in >= 2 libraries,
#                       re-centred on summits, 401 bp (summits = 200)
#   our_union           the step 16 union regions, also re-centred on summits, 401 bp
# Normalizations: libsize (DiffBind 3 default: full library size), rip_native (DESeq2 median of
#   ratios on reads in the regions), background (DESeq2 native on 15 kb genome bins)
# Designs: ~ Replicate + Factor (Factor = day_genotype): KO vs Ctr at each day, time within each
#   genotype; ~ Replicate + Treatment + Condition: shared KO vs Ctr. DESeq2, FDR < 0.05, no fold cut.
# Compared with our DESeq2 regions (step 19 padj-only tests, step 17 time tables): significant
# DiffBind sites overlapping our significant regions, and the reverse.
# Output in <diff_dir>/diffbind/: summary.tsv, comparison.tsv, sites/<set>_<norm>_<test>.tsv, report.html
# usage: Rscript histone_diffbind.R <results/diff/<mark>>

suppressPackageStartupMessages({
  library(DiffBind)
  library(data.table)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
dd <- args[1]
mark <- basename(dd)
od <- dirname(dirname(dd))                       # results/
out <- file.path(dd, "diffbind"); dir.create(file.path(out, "sites"), recursive = TRUE, showWarnings = FALSE)
TH <- 0.05                                       # FDR threshold
# DiffBind forks one counting worker per core, each holding a BAM's alignments: 12 workers ran out
# of 240 GB (jobs 35500667/8), so counting uses DB_CORES (default 4).
ncores <- as.integer(Sys.getenv("DB_CORES", "4"))
# Region sets to run (both by default); a re-run of one set keeps the other's site tables, and the
# summary and comparison are rebuilt from all site tables on disk.
sets <- strsplit(Sys.getenv("DB_SETS", "diffbind_consensus,our_union"), ",")[[1]]

# Per-library peaks and summit re-centring, from the environment (defaults = sharp marks).
# Broad marks: DB_PEAK_DIR=macs3_broad DB_PEAK_SUFFIX=_peaks.filtered.broadPeak DB_PEAK_FORMAT=bed
# DB_SUMMITS=0 (no re-centring: sites keep their full width).
peak_dir <- Sys.getenv("DB_PEAK_DIR", "macs3")
peak_suffix <- Sys.getenv("DB_PEAK_SUFFIX", "_peaks.filtered.narrowPeak")
peak_format <- Sys.getenv("DB_PEAK_FORMAT", "narrow")
summits <- as.integer(Sys.getenv("DB_SUMMITS", "200"))
# Summit re-centring for our_union separately (DB_UNION_SUMMITS=0: full-width union regions, as in
# our DESeq2 analysis); defaults to DB_SUMMITS.
union_summits <- as.integer(Sys.getenv("DB_UNION_SUMMITS", as.character(summits)))
cat(sprintf("peaks: %s/<sample>%s (%s); summits: %s\n", peak_dir, peak_suffix, peak_format, if (summits > 0) summits else "none"))

samples <- fread(file.path(dd, "samples.tsv"))
days <- intersect(c("D9", "D15", "D18", "D22"), unique(samples$day))
ss <- data.frame(SampleID = samples$sample, Tissue = samples$day, Factor = samples$group,
                 Condition = samples$genotype, Treatment = samples$day, Replicate = samples$replicate,
                 bamReads = samples$bam,
                 Peaks = file.path(od, peak_dir, paste0(samples$sample, peak_suffix)),
                 PeakCaller = peak_format, stringsAsFactors = FALSE)
stopifnot(all(file.exists(ss$Peaks)), all(file.exists(ss$bamReads)))

regions <- fread(file.path(dd, "regions.bed"), header = FALSE, col.names = c("chr", "start", "end", "region", "n_groups", "groups"))
union_gr <- GRanges(regions$chr, IRanges(regions$start + 1, regions$end))

ko_tests <- setNames(lapply(days, function(d) c(paste0(d, "_ISL1KO"), paste0(d, "_H3CTR"))), paste0("KO_vs_Ctr_", days))
time_tests <- list()
for (g in c("H3CTR", "ISL1KO")) for (ab in list(c("D18", "D15"), c("D22", "D18"), c("D22", "D15")))
  if (all(ab %in% days)) time_tests[[sprintf("time_%s_%s_vs_%s", g, ab[1], ab[2])]] <- paste0(ab, "_", g)
group_tests <- c(ko_tests, time_tests)

norms <- list(
  libsize    = list(normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL, background = FALSE),
  rip_native = list(normalize = DBA_NORM_NATIVE, library = DBA_LIBSIZE_PEAKREADS, background = FALSE),
  background = list(normalize = DBA_NORM_NATIVE, library = DBA_LIBSIZE_BACKGROUND, background = TRUE))

count_set <- function(set) {
  db <- dba(sampleSheet = ss, minOverlap = 2)
  db$config$cores <- ncores
  db$config$singleEnd <- FALSE
  # peaks must be omitted (not NULL) for DiffBind's own consensus: peaks = NULL means "re-use counts"
  sm <- if (set == "our_union") union_summits else summits
  a <- list(db, summits = if (sm > 0) sm else FALSE, bUseSummarizeOverlaps = TRUE,
            bRemoveDuplicates = FALSE, minCount = 0)
  if (set == "our_union") a$peaks <- union_gr
  db <- do.call(dba.count, a)
  cat(sprintf("%s: %d sites counted\n", set, nrow(dba.peakset(db, bRetrieve = TRUE, DataType = DBA_DATA_FRAME))))
  db
}

summ <- list(); sig_sites <- list()
for (set in sets) {
  db <- count_set(set)
  saveRDS(db, file.path(out, paste0(set, ".counted.rds")))
  for (nm in names(norms)) {
    n <- norms[[nm]]
    dbn <- dba.normalize(db, method = DBA_DESEQ2, normalize = n$normalize, library = n$library, background = n$background)
    # per-group contrasts
    d1 <- NULL
    for (t in names(group_tests)) {
      g <- group_tests[[t]]
      d1 <- if (is.null(d1)) dba.contrast(dbn, design = "~ Replicate + Factor", contrast = c("Factor", g[1], g[2]))
            else dba.contrast(d1, contrast = c("Factor", g[1], g[2]))
    }
    d1 <- dba.analyze(d1, method = DBA_DESEQ2, bBlacklist = FALSE, bGreylist = FALSE)
    d2 <- dba.contrast(dbn, design = "~ Replicate + Treatment + Condition", contrast = c("Condition", "ISL1KO", "H3CTR"))
    d2 <- dba.analyze(d2, method = DBA_DESEQ2, bBlacklist = FALSE, bGreylist = FALSE)
    reps <- c(setNames(lapply(seq_along(group_tests), function(i) list(db = d1, i = i)), names(group_tests)),
              list(KO_vs_Ctr_shared = list(db = d2, i = 1)))
    for (t in names(reps)) {
      r <- dba.report(reps[[t]]$db, contrast = reps[[t]]$i, th = 1, method = DBA_DESEQ2, bCounts = FALSE)
      dt <- as.data.table(as.data.frame(r))
      fwrite(dt, file.path(out, "sites", paste0(set, "_", nm, "_", t, ".tsv")), sep = "\t")
      s <- dt[FDR < TH]
      summ[[length(summ) + 1]] <- data.table(sites = set, normalization = nm, test = t, tested = nrow(dt),
                                             up = sum(s$Fold > 0), down = sum(s$Fold < 0))
      sig_sites[[paste(set, nm, t)]] <- GRanges(s$seqnames, IRanges(s$start, s$end))
      cat(sprintf("%s %s %s: %d up, %d down (FDR < %.2f)\n", set, nm, t, sum(s$Fold > 0), sum(s$Fold < 0), TH))
    }
  }
}
# Rebuild from every site table on disk (covers sets run in earlier jobs).
summ <- list(); sig_sites <- list()
for (f in list.files(file.path(out, "sites"), pattern = "\\.tsv$", full.names = TRUE)) {
  m <- regmatches(basename(f), regexec("^(diffbind_consensus|our_union)_(libsize|rip_native|background)_(.+)\\.tsv$", basename(f)))[[1]]
  if (length(m) != 4) next
  dt <- fread(f)
  s <- dt[FDR < TH & seqnames %in% paste0("chr", c(1:22, "X", "Y"))]
  summ[[length(summ) + 1]] <- data.table(sites = m[2], normalization = m[3], test = m[4], tested = nrow(dt),
                                         up = sum(s$Fold > 0), down = sum(s$Fold < 0))
  sig_sites[[paste(m[2], m[3], m[4])]] <- GRanges(s$seqnames, IRanges(s$start, s$end))
}
summ_dt <- rbindlist(summ)[order(sites, normalization, test)]
fwrite(summ_dt, file.path(out, "summary.tsv"), sep = "\t")

# ---- comparison with our DESeq2 region-level calls (padj only)
ours <- list()
st <- fread(file.path(dd, "sensitive", "regions_all_tests.tsv"))
for (t in c(names(ko_tests), "KO_vs_Ctr_shared")) {
  k <- sub("KO_vs_Ctr_", "", t)
  ours[[t]] <- st[!is.na(get(paste0("padj_", k))) & get(paste0("padj_", k)) < TH, region]
}
for (t in names(time_tests)) {
  f <- file.path(dd, "deseq2", "tables", paste0(t, ".tsv"))
  ours[[t]] <- if (file.exists(f)) fread(f)[!is.na(padj) & padj < TH, region] else character(0)
}
reg_gr <- GRanges(regions$chr, IRanges(regions$start + 1, regions$end), region = regions$region)
comp <- list()
for (k in names(sig_sites)) {
  p <- strsplit(k, " ")[[1]]; t <- p[3]
  o <- reg_gr[reg_gr$region %in% ours[[t]]]
  s <- sig_sites[[k]]
  comp[[length(comp) + 1]] <- data.table(sites = p[1], normalization = p[2], test = t, diffbind_sig = length(s), ours_sig = length(o),
    diffbind_in_ours = sum(overlapsAny(s, o)), ours_in_diffbind = sum(overlapsAny(o, s)))
}
comp_dt <- rbindlist(comp)
fwrite(comp_dt, file.path(out, "comparison.tsv"), sep = "\t")
print(summ_dt); print(comp_dt)

htab <- function(d) {
  d <- as.data.frame(d); d[is.na(d)] <- ""
  paste0("<table><tr>", paste0("<th>", names(d), "</th>", collapse = ""), "</tr>",
         paste0(apply(d, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>")), collapse = ""), "</table>")
}
writeLines(paste0('<!doctype html><html><head><meta charset="utf-8"><title>', mark, ' DiffBind</title>
<style>body{font-family:sans-serif;max-width:1200px;margin:24px auto;padding:0 16px}table{border-collapse:collapse;font-size:12px;margin-bottom:16px}
td,th{border:1px solid #ccc;padding:3px 6px;text-align:right}</style></head><body>
<h1>', mark, ': DiffBind (', as.character(packageVersion("DiffBind")), ') cross-check</h1>
<p>DESeq2, FDR &lt; ', TH, ', no fold-change cut; summits = ', if (summits > 0) summits else "none", ' (DiffBind consensus), ', if (union_summits > 0) union_summits else "none", ' (our union). up = higher in KO, or in the later day.
diffbind_consensus: per-library MACS3 peaks in &ge; 2 libraries; our_union: the step 16 union regions, re-centred on summits.</p>
<h2>Counts</h2>', htab(summ_dt), '
<h2>Overlap with our DESeq2 region calls (padj &lt; ', TH, ')</h2>', htab(comp_dt), '
</body></html>'), file.path(out, "report.html"))
cat("done\n")
