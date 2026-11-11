#!/usr/bin/env Rscript
# Compare DESeq2 (results/DE_DESeq2/) and limma-voom (results/DE_limma/) test by test.
#
# For each of the 11 tests, genes are matched by gene_id (both use the same gene filter).
# Significant = padj < 0.05 in that method.
#   n_both_tested, n_sig_DESeq2, n_sig_limma, n_sig_both, jaccard (both / either),
#   pct_DESeq2_confirmed (both / DESeq2 sig), pct_limma_confirmed (both / limma sig),
#   Pairwise tests also: pearson_log2FC (all genes), spearman_log2FC, direction_agree_pct
#   (same sign among genes significant in both). F/LRT tests: spearman of -log10 p.
# Output: results/DE_comparison/summary.tsv and one merged table per test,
#         compare_<test>.tsv (gene, both methods' statistics, `call`: both / DESeq2_only /
#         limma_only / neither).
# Usage:  Rscript pipeline/08_compare_DE_methods.R <folder root>   (via 08_compare_DE_methods.sb)

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
d_dir <- file.path(root, "results", "DE_DESeq2")
l_dir <- file.path(root, "results", "DE_limma")
out_dir <- file.path(root, "results", "DE_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05

tests <- c("KO_vs_Ctr_D15", "KO_vs_Ctr_D18", "KO_vs_Ctr_D22",
           "Ctr_D18_vs_D15", "Ctr_D22_vs_D15", "Ctr_D22_vs_D18",
           "KO_D18_vs_D15", "KO_D22_vs_D15", "KO_D22_vs_D18",
           "LRT_interaction", "LRT_genotype")

rows <- list()
for (t in tests) {
  d <- fread(file.path(d_dir, paste0(t, ".tsv")), sep = "\t")
  l <- fread(file.path(l_dir, paste0(t, ".tsv")), sep = "\t")
  pairwise <- "log2FoldChange" %in% names(l)
  dk <- c("gene_id", "gene_name", "baseMean", if (pairwise) c("log2FoldChange", "log2FC_shrunk"), "pvalue", "padj")
  lk <- c("gene_id", "AveExpr", if (pairwise) "log2FoldChange", "pvalue", "padj")
  m <- merge(d[, ..dk], l[, ..lk], by = "gene_id", suffixes = c("_DESeq2", "_limma"))
  setnames(m, c("baseMean", "AveExpr"), c("baseMean_DESeq2", "AveExpr_limma"))
  if (pairwise) setnames(m, "log2FC_shrunk", "log2FC_shrunk_DESeq2")
  sd_ <- !is.na(m$padj_DESeq2) & m$padj_DESeq2 < ALPHA
  sl_ <- !is.na(m$padj_limma) & m$padj_limma < ALPHA
  m[, call := fifelse(sd_ & sl_, "both", fifelse(sd_, "DESeq2_only", fifelse(sl_, "limma_only", "neither")))]
  m <- m[order(match(call, c("both", "DESeq2_only", "limma_only", "neither")), padj_DESeq2, gene_id, method = "radix")]
  fwrite(m, file.path(out_dir, paste0("compare_", t, ".tsv")), sep = "\t", na = "NA", quote = FALSE)

  sd_ <- !is.na(m$padj_DESeq2) & m$padj_DESeq2 < ALPHA
  sl_ <- !is.na(m$padj_limma) & m$padj_limma < ALPHA
  both <- sum(sd_ & sl_)
  r <- data.frame(test = t, n_both_tested = sum(!is.na(m$padj_DESeq2) & !is.na(m$padj_limma)),
                  n_sig_DESeq2 = sum(sd_), n_sig_limma = sum(sl_), n_sig_both = both,
                  jaccard = round(both / max(1, sum(sd_ | sl_)), 3),
                  pct_DESeq2_confirmed = round(100 * both / max(1, sum(sd_)), 1),
                  pct_limma_confirmed = round(100 * both / max(1, sum(sl_)), 1))
  if (pairwise) {
    ok <- is.finite(m$log2FoldChange_DESeq2) & is.finite(m$log2FoldChange_limma)
    r$pearson_log2FC  <- round(cor(m$log2FoldChange_DESeq2[ok], m$log2FoldChange_limma[ok]), 3)
    r$spearman_log2FC <- round(cor(m$log2FoldChange_DESeq2[ok], m$log2FoldChange_limma[ok], method = "spearman"), 3)
    bb <- sd_ & sl_
    r$direction_agree_pct <- round(100 * mean(sign(m$log2FoldChange_DESeq2[bb]) == sign(m$log2FoldChange_limma[bb])), 1)
    r$spearman_log10p <- NA
  } else {
    r$pearson_log2FC <- NA; r$spearman_log2FC <- NA; r$direction_agree_pct <- NA
    ok <- is.finite(m$pvalue_DESeq2) & is.finite(m$pvalue_limma)
    r$spearman_log10p <- round(cor(-log10(m$pvalue_DESeq2[ok]), -log10(m$pvalue_limma[ok]), method = "spearman"), 3)
  }
  rows[[t]] <- r
}
s <- do.call(rbind, rows)
fwrite(s, file.path(out_dir, "summary.tsv"), sep = "\t", na = "NA", quote = FALSE)
print(s, row.names = FALSE)
message("Done: ", out_dir)
