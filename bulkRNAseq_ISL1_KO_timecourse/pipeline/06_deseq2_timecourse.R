#!/usr/bin/env Rscript
# DESeq2, ISL1 KO time course (YJZ5LS): genotype (Ctr, KO) x day (D15, D18, D22), n = 2.
#
# Model: ~ genotype + day + genotype:day, fitted as the equivalent 6-level `group`
# (Ctr_D15 ... KO_D22) so every comparison is a named contrast. Day is a factor, not a
# number. No replicate term: rep1/rep2 are not differentiation batches.
#
# Tests (all on one fit of all 12 samples, so dispersions are shared):
#   Wald  KO vs Ctr at D15, D18, D22                          -> KO_vs_Ctr_<day>.tsv
#   Wald  day vs day within each genotype                      -> <geno>_<dayB>_vs_<dayA>.tsv
#   LRT   full vs ~ genotype + day  (KO effect changes w/ day) -> LRT_interaction.tsv
#   LRT   full vs ~ day             (KO effect at any day)     -> LRT_genotype.tsv
# Wald tables: log2FoldChange/padj from the MLE fit, log2FC_shrunk from apeglm (for ranking
# and plots; significance is not taken from the shrunk values).
#
# Input:  results/counts/featureCounts.gene.txt (UMI-deduplicated, featureCounts -s 1)
#         results/counts/gene_annotation.tsv    (gene_id, gene_name, gene_type; from the GTF)
#         sample_sheet.csv
# Output: results/DE_DESeq2/
# Usage:  Rscript pipeline/06_deseq2_timecourse.R <folder root>   (via 06_deseq2_timecourse.sb)

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(data.table)
})
set.seed(12345)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
count_file <- file.path(root, "results", "counts", "featureCounts.gene.txt")
annot_file <- file.path(root, "results", "counts", "gene_annotation.tsv")
sheet_file <- file.path(root, "sample_sheet.csv")
out_dir    <- file.path(root, "results", "DE_DESeq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ALPHA <- 0.05

# --- Counts and metadata -------------------------------------------------------
fc <- fread(count_file, skip = "Geneid", sep = "\t", check.names = FALSE)
counts <- as.matrix(fc[, -(1:6)])
rownames(counts) <- fc$Geneid
colnames(counts) <- sub("\\.dedup\\.bam$", "", basename(colnames(counts)))

meta <- read.csv(sheet_file, stringsAsFactors = FALSE)
stopifnot(setequal(meta$sample, colnames(counts)))
meta$genotype <- factor(meta$genotype, levels = c("Ctr", "KO"))
meta$day      <- factor(meta$day, levels = c("D15", "D18", "D22"))
meta$group    <- factor(paste(meta$genotype, meta$day, sep = "_"),
                        levels = as.vector(t(outer(levels(meta$genotype), levels(meta$day), paste, sep = "_"))))
rownames(meta) <- meta$sample
counts <- counts[, meta$sample]
print(table(meta$genotype, meta$day))

annot <- fread(annot_file, sep = "\t")
setkey(annot, gene_id)

# Keep genes with >= 10 counts in >= 2 samples (2 = smallest group size).
keep <- rowSums(counts >= 10) >= 2
message(sprintf("Genes kept: %d of %d", sum(keep), length(keep)))
counts <- counts[keep, ]

# --- Fit -------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ group)
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = TRUE)

add_annot <- function(df) {
  a <- annot[J(df$gene_id), .(gene_name, gene_type)]
  cbind(df[, "gene_id", drop = FALSE], a, df[, setdiff(names(df), "gene_id"), drop = FALSE])
}
sort_res <- function(df) {
  # Deterministic: padj, then pvalue, then gene_id (NA padj last).
  df[order(is.na(df$padj), df$padj, df$pvalue, df$gene_id, method = "radix"), ]
}
write_tsv <- function(df, name) {
  fwrite(df, file.path(out_dir, name), sep = "\t", na = "NA")
}

summary_rows <- list()
summarise <- function(name, df, has_lfc = TRUE) {
  sig <- !is.na(df$padj) & df$padj < ALPHA
  summary_rows[[name]] <<- data.frame(
    test = name, n_tested = sum(!is.na(df$padj)), n_padj_lt_0.05 = sum(sig),
    n_up = if (has_lfc) sum(sig & df$log2FoldChange > 0) else NA,
    n_down = if (has_lfc) sum(sig & df$log2FoldChange < 0) else NA)
}

# --- Wald contrasts ------------------------------------------------------------------
# apeglm needs a model coefficient, so for each contrast the group factor is releveled to
# the reference group and the Wald test re-run (same dispersions, same fit).
wald_contrast <- function(num, den) {
  d <- dds
  d$group <- relevel(d$group, ref = den)
  d <- nbinomWaldTest(d, quiet = TRUE)
  coef <- paste0("group_", num, "_vs_", den)
  stopifnot(coef %in% resultsNames(d))
  res <- results(d, name = coef, alpha = ALPHA)
  shr <- lfcShrink(d, coef = coef, res = res, type = "apeglm", quiet = TRUE)
  df <- data.frame(gene_id = rownames(res), baseMean = res$baseMean,
                   log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                   log2FC_shrunk = shr$log2FoldChange,
                   stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  sort_res(add_annot(df))
}

contrasts <- list(
  KO_vs_Ctr_D15 = c("KO_D15", "Ctr_D15"),
  KO_vs_Ctr_D18 = c("KO_D18", "Ctr_D18"),
  KO_vs_Ctr_D22 = c("KO_D22", "Ctr_D22"),
  Ctr_D18_vs_D15 = c("Ctr_D18", "Ctr_D15"),
  Ctr_D22_vs_D15 = c("Ctr_D22", "Ctr_D15"),
  Ctr_D22_vs_D18 = c("Ctr_D22", "Ctr_D18"),
  KO_D18_vs_D15 = c("KO_D18", "KO_D15"),
  KO_D22_vs_D15 = c("KO_D22", "KO_D15"),
  KO_D22_vs_D18 = c("KO_D22", "KO_D18")
)
for (nm in names(contrasts)) {
  message("Wald: ", nm)
  df <- wald_contrast(contrasts[[nm]][1], contrasts[[nm]][2])
  write_tsv(df, paste0(nm, ".tsv"))
  summarise(nm, df)
}

# --- Likelihood-ratio tests ------------------------------------------------------------
lrt <- function(reduced) {
  d <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ genotype + day + genotype:day)
  d <- DESeq(d, test = "LRT", reduced = reduced, quiet = TRUE)
  res <- results(d, alpha = ALPHA)
  df <- data.frame(gene_id = rownames(res), baseMean = res$baseMean,
                   stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  sort_res(add_annot(df))
}
message("LRT: interaction")
df <- lrt(~ genotype + day); write_tsv(df, "LRT_interaction.tsv"); summarise("LRT_interaction", df, FALSE)
message("LRT: genotype")
df <- lrt(~ day);            write_tsv(df, "LRT_genotype.tsv");    summarise("LRT_genotype", df, FALSE)

# --- Tables for QC and the report --------------------------------------------------------
write_tsv(do.call(rbind, summary_rows), "summary.tsv")

nc <- data.frame(gene_id = rownames(norm_counts), round(norm_counts, 3), check.names = FALSE)
write_tsv(add_annot(nc), "normalized_counts.tsv")
vm <- data.frame(gene_id = rownames(vsd), round(assay(vsd), 4), check.names = FALSE)
write_tsv(add_annot(vm), "vst.tsv")

# PCA on the 500 most variable VST genes (ties broken by gene_id).
v <- assay(vsd)
rv <- apply(v, 1, var)
top <- order(-rv, rownames(v), method = "radix")[seq_len(min(500, nrow(v)))]
pca <- prcomp(t(v[top, ]))
pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pc <- data.frame(sample = rownames(pca$x), pca$x[, 1:4], meta[rownames(pca$x), c("genotype", "day", "replicate")])
write_tsv(pc, "pca_top500.tsv")
write_tsv(data.frame(PC = paste0("PC", seq_along(pct)), pct_variance = round(pct, 2)), "pca_variance.tsv")

cm <- cor(v, method = "pearson")
write_tsv(data.frame(sample = rownames(cm), round(cm, 5), check.names = FALSE), "sample_correlation_vst.tsv")
write_tsv(data.frame(sample = colnames(dds), size_factor = round(sizeFactors(dds), 5)), "size_factors.tsv")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(do.call(rbind, summary_rows), row.names = FALSE)
message("Done: ", out_dir)
