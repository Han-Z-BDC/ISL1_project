#!/usr/bin/env Rscript
# DESeq2 differential ISL1 binding on one consensus region set (10_consensus_counts.sb).
#   design  ~ batch + group   (batch = Admera delivery; group = day x genotype, reference D12_Ctr)
#   size factors: DESeq2 median-of-ratios on the consensus counts
#   contrasts: every pair of Ctr stages; KO+rescue vs Ctr at each day that has both
# Writes to <consensus_dir>/deseq2/:
#   dds.rds, samples.tsv (with size factors and fragments counted), pca.tsv (+ pca_variance.tsv,
#   top 500 most variable regions of the VST), de_<A>_vs_<B>.tsv (all regions, annotated, by padj),
#   de_summary.tsv (regions up / down at padj < 0.05 and |log2FC| >= 1, median log2FC)
# With KO+rescue groups and <qc_dir>, the KO+rescue vs Ctr contrasts are also run with library-size
# normalization (size factors = total filtered fragments, from <qc_dir>/<sample>.flagstat.txt
# read1, scaled to geometric mean 1), so a global loss of binding is not normalized away:
#   deseq2/libsize/{samples.tsv, de_<A>_vs_<B>.tsv, de_summary.tsv}
# usage: Rscript deseq2_consensus.R <consensus_dir> [<qc_dir>]

suppressPackageStartupMessages({
  library(DESeq2)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2)
cdir <- args[1]
qc_dir <- if (length(args) == 2) args[2] else NA
out <- file.path(cdir, "deseq2")
dir.create(out, showWarnings = FALSE)
PADJ <- 0.05; LFC <- 1

counts <- fread(file.path(cdir, "counts.tsv"))
regions <- fread(file.path(cdir, "regions.bed"), header = FALSE,
                 col.names = c("chr", "start", "end", "region", "groups", "n_groups"))
samples <- fread(file.path(cdir, "samples.tsv"))
m <- as.matrix(counts[, -1]); rownames(m) <- counts$region
stopifnot(identical(colnames(m), samples$sample))

groups <- intersect(c("D12_Ctr", "D15_Ctr", "D22_Ctr", "D12_KO", "D15_KO"), unique(samples$group))
cd <- data.frame(row.names = samples$sample, group = factor(samples$group, levels = groups),
                 batch = factor(samples$batch))
dds <- DESeqDataSetFromMatrix(m, cd, design = ~ batch + group)
dds <- DESeq(dds, quiet = TRUE)
saveRDS(dds, file.path(out, "dds.rds"))

samples[, size_factor := sizeFactors(dds)[sample]]
samples[, fragments_counted := colSums(m)[sample]]
fwrite(samples, file.path(out, "samples.tsv"), sep = "\t")

vsd <- vst(dds, blind = TRUE)
v <- assay(vsd)
top <- head(order(apply(v, 1, var), decreasing = TRUE), 500)
pc <- prcomp(t(v[top, ]))
pv <- 100 * pc$sdev^2 / sum(pc$sdev^2)
fwrite(data.table(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3]),
       file.path(out, "pca.tsv"), sep = "\t")
fwrite(data.table(pc = paste0("PC", seq_along(pv)), percent = pv), file.path(out, "pca_variance.tsv"), sep = "\t")

# Annotate each region at its centre: nearest gene and distance to its TSS.
gr <- GRanges(regions$chr, IRanges(as.integer((regions$start + regions$end) / 2) + 1, width = 1), region = regions$region)
ann <- as.data.table(as.data.frame(annotatePeak(gr, tssRegion = c(-3000, 3000),
                                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                                annoDb = "org.Hs.eg.db", verbose = FALSE)))
ann <- ann[, .(region, annotation, SYMBOL, distanceToTSS)]
regions <- merge(regions, ann, by = "region", all.x = TRUE)

contrasts <- list()
ctr <- grep("_Ctr$", groups, value = TRUE)
for (i in seq_along(ctr)) for (j in seq_along(ctr)) if (i > j) contrasts[[length(contrasts) + 1]] <- c(ctr[i], ctr[j])
for (ko in grep("_KO$", groups, value = TRUE)) {
  c0 <- sub("_KO$", "_Ctr", ko)
  if (c0 %in% groups) contrasts[[length(contrasts) + 1]] <- c(ko, c0)
}

run_contrasts <- function(dds, contrasts, odir) {
  dir.create(odir, showWarnings = FALSE, recursive = TRUE)
  summ <- list()
  for (cc in contrasts) {
    r <- results(dds, contrast = c("group", cc[1], cc[2]), alpha = PADJ)
    d <- data.table(region = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange,
                    lfcSE = r$lfcSE, pvalue = r$pvalue, padj = r$padj)
    d <- merge(regions, d, by = "region")
    setorder(d, padj, na.last = TRUE)
    name <- paste0(cc[1], "_vs_", cc[2])
    fwrite(d, file.path(odir, paste0("de_", name, ".tsv")), sep = "\t")
    sig <- d[!is.na(padj) & padj < PADJ & abs(log2FC) >= LFC]
    summ[[length(summ) + 1]] <- data.table(contrast = name, a = cc[1], b = cc[2], tested = sum(!is.na(d$padj)),
                                           up = sum(sig$log2FC > 0), down = sum(sig$log2FC < 0),
                                           median_log2FC = median(d$log2FC, na.rm = TRUE))
    cat(sprintf("%s %s %s: %d up, %d down (padj < %.2f, |log2FC| >= %g), median log2FC %.2f\n", basename(cdir),
                basename(odir), name, sum(sig$log2FC > 0), sum(sig$log2FC < 0), PADJ, LFC, median(d$log2FC, na.rm = TRUE)))
  }
  fwrite(rbindlist(summ), file.path(odir, "de_summary.tsv"), sep = "\t")
}

run_contrasts(dds, contrasts, out)

# KO+rescue vs Ctr again with library-size normalization: same model and counts, size factors
# from total filtered fragments instead of median-of-ratios on the regions.
ko_contrasts <- Filter(function(cc) grepl("_KO$", cc[1]), contrasts)
if (length(ko_contrasts) && !is.na(qc_dir)) {
  total <- sapply(samples$sample, function(s) {
    x <- grep(" read1$", readLines(file.path(qc_dir, paste0(s, ".flagstat.txt"))), value = TRUE)
    as.numeric(sub(" .*", "", x))
  })
  sf <- total / exp(mean(log(total)))
  dds_lib <- DESeqDataSetFromMatrix(m, cd, design = ~ batch + group)
  sizeFactors(dds_lib) <- sf
  dds_lib <- DESeq(dds_lib, quiet = TRUE)
  lib_dir <- file.path(out, "libsize")
  dir.create(lib_dir, showWarnings = FALSE)
  saveRDS(dds_lib, file.path(lib_dir, "dds.rds"))
  fwrite(data.table(sample = samples$sample, group = samples$group, total_fragments = total, size_factor = sf),
         file.path(lib_dir, "samples.tsv"), sep = "\t")
  run_contrasts(dds_lib, ko_contrasts, lib_dir)
}
