#!/usr/bin/env Rscript
# limma-voom cross-check of the DESeq2 time-course analysis (06_deseq2_timecourse.R).
#
# Same input counts, same gene filter (>= 10 counts in >= 2 samples), same 6-level group
# model (~ 0 + group, Ctr_D15 ... KO_D22), same contrasts, so results can be matched gene by
# gene in 08_compare_DE_methods.R. TMM normalisation (edgeR), voom, lmFit, eBayes.
#
# Tests (file names match results/DE_DESeq2/):
#   moderated t  KO vs Ctr at D15, D18, D22; day vs day within each genotype
#   F-test       LRT_interaction: (KO-Ctr)_D18 - (KO-Ctr)_D15 and (KO-Ctr)_D22 - (KO-Ctr)_D15 = 0
#                (the counterpart of the DESeq2 LRT full vs ~ genotype + day)
#   F-test       LRT_genotype: KO-Ctr at D15, D18 and D22 all = 0 (DESeq2 LRT full vs ~ day)
# Columns: gene_id, gene_name, gene_type, AveExpr (log2 CPM), log2FoldChange, stat (t or F),
#          pvalue, padj (BH). F-test tables have no log2FoldChange.
#
# Input:  results/counts/featureCounts.gene.txt, results/counts/gene_annotation.tsv, sample_sheet.csv
# Output: results/DE_limma/
# Usage:  Rscript pipeline/07_limma_voom_timecourse.R <folder root>   (via 07_limma_voom_timecourse.sb)

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
count_file <- file.path(root, "results", "counts", "featureCounts.gene.txt")
annot_file <- file.path(root, "results", "counts", "gene_annotation.tsv")
sheet_file <- file.path(root, "sample_sheet.csv")
out_dir    <- file.path(root, "results", "DE_limma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05

# --- Counts and metadata (identical to 06) ------------------------------------------
fc <- fread(count_file, skip = "Geneid", sep = "\t", check.names = FALSE)
counts <- as.matrix(fc[, -(1:6)])
rownames(counts) <- fc$Geneid
colnames(counts) <- sub("\\.dedup\\.bam$", "", basename(colnames(counts)))

meta <- read.csv(sheet_file, stringsAsFactors = FALSE)
stopifnot(setequal(meta$sample, colnames(counts)))
meta$genotype <- factor(meta$genotype, levels = c("Ctr", "KO"))
meta$day      <- factor(meta$day, levels = c("D15", "D18", "D22"))
grp_levels    <- as.vector(t(outer(levels(meta$genotype), levels(meta$day), paste, sep = "_")))
meta$group    <- factor(paste(meta$genotype, meta$day, sep = "_"), levels = grp_levels)
rownames(meta) <- meta$sample
counts <- counts[, meta$sample]

annot <- fread(annot_file, sep = "\t")
setkey(annot, gene_id)

keep <- rowSums(counts >= 10) >= 2
message(sprintf("Genes kept: %d of %d", sum(keep), length(keep)))
counts <- counts[keep, ]

# --- Fit ----------------------------------------------------------------------------
dge <- DGEList(counts = counts, group = meta$group)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

cm <- makeContrasts(
  KO_vs_Ctr_D15  = KO_D15 - Ctr_D15,
  KO_vs_Ctr_D18  = KO_D18 - Ctr_D18,
  KO_vs_Ctr_D22  = KO_D22 - Ctr_D22,
  Ctr_D18_vs_D15 = Ctr_D18 - Ctr_D15,
  Ctr_D22_vs_D15 = Ctr_D22 - Ctr_D15,
  Ctr_D22_vs_D18 = Ctr_D22 - Ctr_D18,
  KO_D18_vs_D15  = KO_D18 - KO_D15,
  KO_D22_vs_D15  = KO_D22 - KO_D15,
  KO_D22_vs_D18  = KO_D22 - KO_D18,
  INT_D18        = (KO_D18 - Ctr_D18) - (KO_D15 - Ctr_D15),
  INT_D22        = (KO_D22 - Ctr_D22) - (KO_D15 - Ctr_D15),
  levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))

add_annot <- function(df) {
  a <- annot[J(df$gene_id), .(gene_name, gene_type)]
  cbind(df[, "gene_id", drop = FALSE], a, df[, setdiff(names(df), "gene_id"), drop = FALSE])
}
sort_res <- function(df) df[order(is.na(df$padj), df$padj, df$pvalue, df$gene_id, method = "radix"), ]
write_tsv <- function(df, name) fwrite(df, file.path(out_dir, name), sep = "\t", na = "NA", quote = FALSE)

summary_rows <- list()
summarise <- function(name, df, has_lfc = TRUE) {
  sig <- !is.na(df$padj) & df$padj < ALPHA
  summary_rows[[name]] <<- data.frame(
    test = name, n_tested = sum(!is.na(df$padj)), n_padj_lt_0.05 = sum(sig),
    n_up = if (has_lfc) sum(sig & df$log2FoldChange > 0) else NA,
    n_down = if (has_lfc) sum(sig & df$log2FoldChange < 0) else NA)
}

for (nm in colnames(cm)[1:9]) {
  tt <- topTable(fit2, coef = nm, number = Inf, sort.by = "none")
  df <- data.frame(gene_id = rownames(tt), AveExpr = tt$AveExpr, log2FoldChange = tt$logFC,
                   stat = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val)
  df <- sort_res(add_annot(df))
  write_tsv(df, paste0(nm, ".tsv"))
  summarise(nm, df)
}

ftest <- function(coefs, name) {
  tt <- topTable(fit2, coef = coefs, number = Inf, sort.by = "none")
  df <- data.frame(gene_id = rownames(tt), AveExpr = tt$AveExpr, stat = tt$F,
                   pvalue = tt$P.Value, padj = tt$adj.P.Val)
  df <- sort_res(add_annot(df))
  write_tsv(df, paste0(name, ".tsv"))
  summarise(name, df, FALSE)
}
ftest(c("INT_D18", "INT_D22"), "LRT_interaction")
ftest(c("KO_vs_Ctr_D15", "KO_vs_Ctr_D18", "KO_vs_Ctr_D22"), "LRT_genotype")

write_tsv(do.call(rbind, summary_rows), "summary.tsv")
write_tsv(data.frame(sample = colnames(dge), lib_size = dge$samples$lib.size,
                     norm_factor = round(dge$samples$norm.factors, 5)), "tmm_factors.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(do.call(rbind, summary_rows), row.names = FALSE)
message("Done: ", out_dir)
