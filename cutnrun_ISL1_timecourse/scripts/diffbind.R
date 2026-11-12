#!/usr/bin/env Rscript
# Step 07 analysis: ISL1 differential binding over time. Run by scripts/07_diffbind.sb, which
# writes <dir>/libraries.tsv and the featureCounts tables.
#   Libraries  control ISL1, D15 / D18 / D22 rep1-3 (9)
#   Sites      step 06 sites/sites.bed, minus an IgG greylist (GreyListChIP on each trio library's
#              matched IgG, p 0.999, merged over the 9 IgGs; DiffBind's default settings)
#   Design     ~ batch + day, batch = batch 3 vs earlier libraries
#   Contrasts  D18 vs D15, D22 vs D18, D22 vs D15 (log2FC > 0 = higher on the later day)
#   Main       DiffBind (DESeq2): background (DESeq2 native on 15 kb bins) and rle (DESeq2 native
#              on reads in the sites); duplicates kept, paired-end fragments (summarizeOverlaps),
#              no control subtraction
#   Check      DESeq2 on featureCounts counts of the same sites, size factors from 15 kb bins
#              (background) or from the sites (rle); unshrunk log2FC, as DiffBind reports
# FDR < 0.05, with and without |log2FC| >= 1.
# Output in <dir>: greylist.bed, sites_used.bed, samples.tsv, norm_factors.tsv,
#   tables/{diffbind,deseq2}_<norm>_<contrast>.tsv, beds/diffbind_<norm>_<contrast>_{up,down}.bed,
#   summary.tsv, agreement.tsv, pca.tsv, pca_variance.tsv, counts_sites.tsv
# usage: Rscript diffbind.R <results/diff>

suppressPackageStartupMessages({
  library(DiffBind)
  library(DESeq2)
  library(GreyListChIP)
  library(Rsamtools)
  library(GenomicRanges)
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
d <- args[1]
for (s in c("tables", "beds")) dir.create(file.path(d, s), showWarnings = FALSE)
TH <- 0.05
FC1 <- 1
ncores <- as.integer(Sys.getenv("DB_CORES", "4"))
CANON <- paste0("chr", c(1:22, "X", "Y"))
contrasts <- list(D18_vs_D15 = c("D18", "D15"), D22_vs_D18 = c("D22", "D18"), D22_vs_D15 = c("D22", "D15"))
set.seed(1)

# ---- libraries
trio <- fread(file.path(d, "libraries.tsv"), header = FALSE,
              col.names = c("sample", "day", "replicate", "batch", "original_sample", "isl1_bam", "igg_bam"))
stopifnot(all(trio$batch %in% c("batch3", "earlier")))
trio[, day := factor(day, levels = c("D15", "D18", "D22"))]
stopifnot(nrow(trio) == 9, all(file.exists(trio$isl1_bam)), all(file.exists(trio$igg_bam)))
fwrite(trio[, .(sample, day, replicate, batch, original_sample, isl1_bam, igg_bam)], file.path(d, "samples.tsv"), sep = "\t")

# ---- sites and IgG greylist
sites <- fread(file.path(d, "sites", "sites.bed"), header = FALSE,
               col.names = c("chr", "start", "end", "name", "score", "strand", "summit_day", "days_called", "n_days"))
sites_gr <- GRanges(sites$chr, IRanges(sites$start + 1, sites$end), name = sites$name)
hdr <- scanBamHeader(trio$igg_bam[1])[[1]]$targets
hdr <- hdr[names(hdr) %in% CANON]
ki <- Seqinfo(names(hdr), as.integer(hdr))
grey <- GRanges()
for (b in unique(trio$igg_bam)) {
  gl <- new("GreyList", karyotype = ki)
  gl <- countReads(gl, b)
  gl <- calcThreshold(gl, reps = 100, sampleSize = 30000, p = 0.999, cores = 1)
  gl <- makeGreyList(gl, maxGap = 16384)
  cat(sprintf("greylist %s: %d regions, %.2f Mb\n", basename(b), length(gl@regions), sum(width(gl@regions)) / 1e6))
  grey <- c(grey, granges(gl@regions))
}
grey <- reduce(grey)
fwrite(data.table(chr = as.character(seqnames(grey)), start = start(grey) - 1, end = end(grey)),
       file.path(d, "greylist.bed"), sep = "\t", col.names = FALSE)
hit <- overlapsAny(sites_gr, grey)
cat(sprintf("greylist: %d regions, %.2f Mb; %d of %d sites removed\n", length(grey), sum(width(grey)) / 1e6, sum(hit), length(hit)))
sites <- sites[!hit]
sites_gr <- sites_gr[!hit]
fwrite(sites, file.path(d, "sites_used.bed"), sep = "\t", col.names = FALSE)
fwrite(data.table(item = c("greylist regions", "greylist Mb", "sites removed by greylist", "sites tested"),
                  n = c(length(grey), round(sum(width(grey)) / 1e6, 2), sum(hit), nrow(sites))),
       file.path(d, "greylist_summary.tsv"), sep = "\t")

# ---- DiffBind
ss <- data.frame(SampleID = trio$sample, Tissue = "ISL1", Factor = as.character(trio$day),
                 Condition = as.character(trio$day), Treatment = trio$batch, Replicate = trio$replicate,
                 bamReads = trio$isl1_bam, Peaks = file.path(d, "sites_used.bed"), PeakCaller = "bed",
                 stringsAsFactors = FALSE)
db <- dba(sampleSheet = ss, minOverlap = 1)
db$config$cores <- ncores
db$config$singleEnd <- FALSE
db <- dba.count(db, peaks = sites_gr, summits = FALSE, bUseSummarizeOverlaps = TRUE, bRemoveDuplicates = FALSE,
                bSubControl = FALSE, minCount = 0, filter = 0)
pk <- as.data.table(dba.peakset(db, bRetrieve = TRUE, DataType = DBA_DATA_FRAME))
cat(sprintf("DiffBind counted %d sites (expected %d)\n", nrow(pk), length(sites_gr)))
stopifnot(nrow(pk) == length(sites_gr))
saveRDS(db, file.path(d, "diffbind.counted.rds"))
key <- function(chr, start) paste(chr, start, sep = ":")
name_of <- setNames(sites$name, key(sites$chr, sites$start + 1))

norms <- list(
  background = list(normalize = DBA_NORM_NATIVE, library = DBA_LIBSIZE_BACKGROUND, background = TRUE),
  rle        = list(normalize = DBA_NORM_NATIVE, library = DBA_LIBSIZE_PEAKREADS, background = FALSE))

nf <- list(); summ <- list(); sig <- list()
add_summary <- function(method, norm, ct, dt, fdr, lfc) {
  s <- dt[!is.na(get(fdr)) & get(fdr) < TH]
  summ[[length(summ) + 1]] <<- data.table(method = method, normalization = norm, contrast = ct, tested = nrow(dt),
    up = sum(s[[lfc]] > 0), down = sum(s[[lfc]] < 0), up_fc1 = sum(s[[lfc]] >= FC1), down_fc1 = sum(s[[lfc]] <= -FC1))
  sig[[paste(method, norm, ct)]] <<- list(up = s[get(lfc) > 0, name], down = s[get(lfc) < 0, name])
}
for (nm in names(norms)) {
  n <- norms[[nm]]
  dbn <- dba.normalize(db, method = DBA_DESEQ2, normalize = n$normalize, library = n$library, background = n$background)
  f <- dba.normalize(dbn, method = DBA_DESEQ2, bRetrieve = TRUE)
  nf[[length(nf) + 1]] <- data.table(method = "diffbind", normalization = nm, sample = trio$sample,
                                     lib_size = f$lib.sizes, norm_factor = f$norm.factors)
  dc <- NULL
  for (ct in names(contrasts)) {
    g <- contrasts[[ct]]
    dc <- if (is.null(dc)) dba.contrast(dbn, design = "~ Treatment + Condition", contrast = c("Condition", g[1], g[2]))
          else dba.contrast(dc, contrast = c("Condition", g[1], g[2]))
  }
  dc <- dba.analyze(dc, method = DBA_DESEQ2, bBlacklist = FALSE, bGreylist = FALSE)
  for (i in seq_along(contrasts)) {
    ct <- names(contrasts)[i]
    r <- as.data.table(as.data.frame(dba.report(dc, contrast = i, th = 1, method = DBA_DESEQ2, bCounts = FALSE)))
    r[, name := name_of[key(as.character(seqnames), start)]]
    stopifnot(!anyNA(r$name))
    setnames(r, c("Fold", "p.value"), c("log2FC", "pvalue"))
    r <- r[, c("name", "seqnames", "start", "end", "Conc", grep("^Conc_", names(r), value = TRUE), "log2FC", "pvalue", "FDR"), with = FALSE]
    setorder(r, FDR, pvalue, name)
    fwrite(r, file.path(d, "tables", sprintf("diffbind_%s_%s.tsv", nm, ct)), sep = "\t")
    add_summary("diffbind", nm, ct, r, "FDR", "log2FC")
    for (dir in c("up", "down")) {
      s <- r[FDR < TH & (if (dir == "up") log2FC > 0 else log2FC < 0)]
      setorder(s, seqnames, start)
      fwrite(s[, .(seqnames, start - 1, end, name, round(log2FC, 3), ".")], sep = "\t", col.names = FALSE,
             file.path(d, "beds", sprintf("diffbind_%s_%s_%s.bed", nm, ct, dir)))
    }
    cat(sprintf("diffbind %s %s: %d up, %d down\n", nm, ct, length(sig[[paste("diffbind", nm, ct)]]$up), length(sig[[paste("diffbind", nm, ct)]]$down)))
  }
}

# ---- DESeq2 check on featureCounts counts
read_fc <- function(f) {
  x <- fread(f, skip = "Geneid")
  m <- as.matrix(x[, 7:ncol(x)]); rownames(m) <- x$Geneid; colnames(m) <- trio$sample
  stopifnot(identical(basename(names(x)[7:ncol(x)]), basename(trio$isl1_bam)))
  m
}
cs <- read_fc(file.path(d, "featurecounts", "sites.counts.txt"))
cs <- cs[sites$name, , drop = FALSE]
cb <- read_fc(file.path(d, "featurecounts", "bins15k.counts.txt"))
fwrite(data.table(name = rownames(cs), cs), file.path(d, "counts_sites.tsv"), sep = "\t")
cd <- data.frame(row.names = trio$sample, day = trio$day, batch = factor(trio$batch, levels = c("earlier", "batch3")))
sf <- list(background = estimateSizeFactorsForMatrix(cb[rowSums(cb) > 0, ]),
           rle = estimateSizeFactorsForMatrix(cs))
for (nm in names(sf)) {
  nf[[length(nf) + 1]] <- data.table(method = "deseq2", normalization = nm, sample = trio$sample,
                                     lib_size = colSums(cs), norm_factor = sf[[nm]])
  dds <- DESeqDataSetFromMatrix(cs, cd, design = ~ batch + day)
  sizeFactors(dds) <- sf[[nm]]
  dds <- DESeq(dds, quiet = TRUE)
  for (ct in names(contrasts)) {
    g <- contrasts[[ct]]
    res <- results(dds, contrast = c("day", g[1], g[2]), alpha = TH)
    r <- data.table(name = rownames(res), baseMean = res$baseMean, log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
                    pvalue = res$pvalue, padj = res$padj)
    setorder(r, padj, pvalue, name, na.last = TRUE)
    fwrite(r, file.path(d, "tables", sprintf("deseq2_%s_%s.tsv", nm, ct)), sep = "\t", na = "")
    add_summary("deseq2", nm, ct, r, "padj", "log2FC")
  }
}
fwrite(rbindlist(nf), file.path(d, "norm_factors.tsv"), sep = "\t")
summ_dt <- rbindlist(summ)
fwrite(summ_dt, file.path(d, "summary.tsv"), sep = "\t")

# ---- agreement (same sites, same direction)
agr <- list()
for (ct in names(contrasts)) {
  pairs <- list(c("diffbind background", "diffbind rle"), c("diffbind background", "deseq2 background"),
                c("diffbind rle", "deseq2 rle"), c("deseq2 background", "deseq2 rle"))
  for (p in pairs) {
    a <- sig[[paste(p[1], ct)]]; b <- sig[[paste(p[2], ct)]]
    na <- length(a$up) + length(a$down); nb <- length(b$up) + length(b$down)
    both <- length(intersect(a$up, b$up)) + length(intersect(a$down, b$down))
    agr[[length(agr) + 1]] <- data.table(contrast = ct, a = p[1], b = p[2], sig_a = na, sig_b = nb, both_same_direction = both,
      opposite = length(intersect(a$up, b$down)) + length(intersect(a$down, b$up)))
  }
}
fwrite(rbindlist(agr), file.path(d, "agreement.tsv"), sep = "\t")

# ---- PCA (vst of the site counts, rle size factors; before and after removing the batch term)
dds <- DESeqDataSetFromMatrix(cs, cd, design = ~ batch + day)
sizeFactors(dds) <- sf$rle
v <- assay(vst(dds, blind = TRUE))
top <- head(order(apply(v, 1, var), decreasing = TRUE), 5000)
pcs <- list(); pvar <- list()
for (type in c("vst", "vst_batch_removed")) {
  m <- if (type == "vst") v[top, ] else removeBatchEffect(v, batch = cd$batch, design = model.matrix(~ day, cd))[top, ]
  p <- prcomp(t(m), center = TRUE, scale. = FALSE)
  pcs[[type]] <- data.table(type = type, sample = trio$sample, day = as.character(trio$day), batch = trio$batch,
                            PC1 = p$x[, 1], PC2 = p$x[, 2], PC3 = p$x[, 3])
  pvar[[type]] <- data.table(type = type, PC = paste0("PC", 1:3), pct = round(100 * p$sdev[1:3]^2 / sum(p$sdev^2), 1))
}
fwrite(rbindlist(pcs), file.path(d, "pca.tsv"), sep = "\t")
fwrite(rbindlist(pvar), file.path(d, "pca_variance.tsv"), sep = "\t")

print(summ_dt)
print(rbindlist(agr))
cat(sprintf("DiffBind %s, DESeq2 %s\ndone\n", packageVersion("DiffBind"), packageVersion("DESeq2")))
