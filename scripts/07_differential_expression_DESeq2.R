#!/usr/bin/env Rscript
# Differential expression, ISL1 rescue bulk RNA-seq -- DESeq2 version.
#
# Same design and same contrasts as scripts/06_differential_expression.R
# (limma-voom), run through a different statistical engine (negative-binomial
# GLM vs linear model + empirical Bayes on log-CPM) so results can be
# cross-checked -- see scripts/08_compare_DE_methods.R.
#
# Design note: not a full factorial -- Ctr cells only ever get GFP, only KO
# cells get the ISL1-rescue construct -- so we use one 6-level `group` factor
# (day x group3) rather than separate genotype/construct main effects. DESeq2's
# results(contrast=c("group","A","B")) computes any pairwise level difference
# correctly regardless of which level is the reference, so no releveling is
# needed for the 6 primary pairwise contrasts.
#
# Input:  results/rsem/<sample>.genes.results  (same as the limma script)
# Output: results/DE_DESeq2/<contrast>.tsv, results/DE_DESeq2/volcano_<contrast>.pdf
#
# Run in the lab's `r_env` (DESeq2 + apeglm added there specifically for this).
# Usage: Rscript scripts/07_differential_expression_DESeq2.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

# --- Paths ---------------------------------------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
repo_dir     <- normalizePath(file.path(get_script_dir(), ".."))
sample_sheet <- file.path(repo_dir, "scripts", "sample_sheet.csv")
rsem_dir     <- file.path(repo_dir, "results", "rsem")
out_dir      <- file.path(repo_dir, "results", "DE_DESeq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load sample metadata + build the counts matrix -----------------------
meta <- read.csv(sample_sheet, stringsAsFactors = FALSE)

read_rsem_counts <- function(sample) {
  f <- file.path(rsem_dir, paste0(sample, ".genes.results"))
  d <- read.delim(f, stringsAsFactors = FALSE)
  setNames(d$expected_count, d$gene_id)
}
counts_list <- lapply(meta$sample, read_rsem_counts)
stopifnot(length(unique(sapply(counts_list, length))) == 1)
counts <- do.call(cbind, counts_list)
colnames(counts) <- meta$sample
rownames(counts) <- names(counts_list[[1]])
counts <- round(counts)  # DESeq2 requires integer counts

meta$group <- factor(paste(meta$day, meta$group3, sep = "_"))
rownames(meta) <- meta$sample

# --- DESeqDataSet, filter low counts, run DESeq() --------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]  # same spirit as filterByExpr, simpler cutoff
message(sprintf("Testing %d genes after low-count filtering", nrow(dds)))
dds <- DESeq(dds)

# --- PCA sanity check on VST-transformed counts, before trusting any result -
vsd <- vst(dds, blind = TRUE)
pdf(file.path(out_dir, "pca_vst.pdf"), width = 6, height = 5)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
plot(pcaData$PC1, pcaData$PC2, col = as.integer(factor(meta$day)) + 1, pch = 16,
     xlab = "PC1", ylab = "PC2", main = "PCA on VST counts (color = day)")
text(pcaData$PC1, pcaData$PC2, labels = meta$group, pos = 3, cex = 0.6)
dev.off()

# --- Primary pairwise contrasts: direct level-vs-level, any reference ------
pairwise <- list(
  KOeffect_D12     = c("D12_KO_GFP",  "D12_Ctr_GFP"),
  KOeffect_D15     = c("D15_KO_GFP",  "D15_Ctr_GFP"),
  RescueEffect_D12 = c("D12_KO_ISL1", "D12_KO_GFP"),
  RescueEffect_D15 = c("D15_KO_ISL1", "D15_KO_GFP"),
  RescueVsCtrl_D12 = c("D12_KO_ISL1", "D12_Ctr_GFP"),
  RescueVsCtrl_D15 = c("D15_KO_ISL1", "D15_Ctr_GFP")
)

write_result <- function(res, cn) {
  # shrink LFCs (apeglm) for more reliable effect sizes with n=2/group
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[order(res_df$pvalue), ]
  write.table(res_df, file.path(out_dir, paste0(cn, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(out_dir, paste0("volcano_", cn, ".pdf")), width = 5, height = 5)
  sig <- !is.na(res_df$padj) & res_df$padj < 0.05
  plot(res_df$log2FoldChange, -log10(res_df$pvalue),
       col = ifelse(sig, "red", "grey60"), pch = 16, cex = 0.6,
       xlab = "log2 fold change", ylab = "-log10(p-value)", main = cn)
  abline(v = c(-1, 1), lty = 2, col = "grey40")
  dev.off()

  message(sprintf("%-18s  %4d genes FDR<0.05  (of %d tested)",
                   cn, sum(sig, na.rm = TRUE), nrow(res_df)))
}

for (cn in names(pairwise)) {
  lvls <- pairwise[[cn]]
  res <- results(dds, contrast = c("group", lvls[1], lvls[2]), alpha = 0.05)
  write_result(res, cn)
}

# --- Day-difference-of-rescue-effect (interaction), exploratory only -------
# (D15_ISL1 - D15_GFP) - (D12_ISL1 - D12_GFP), built from named model
# coefficients since this isn't a simple two-level contrast. Reference level
# is whatever R's default factor ordering picked; we match coefficient names
# by grep so this doesn't silently break if that ordering changes.
rn <- resultsNames(dds)
find_coef <- function(level) rn[grepl(paste0("group.*", level, "$"), rn)]
c_d15_isl1 <- find_coef("D15_KO_ISL1"); c_d15_gfp <- find_coef("D15_KO_GFP")
c_d12_isl1 <- find_coef("D12_KO_ISL1"); c_d12_gfp <- find_coef("D12_KO_GFP")

if (length(c_d15_isl1) == 1 && length(c_d15_gfp) == 1 &&
    length(c_d12_isl1) == 1 && length(c_d12_gfp) == 1) {
  res_int <- results(dds,
    contrast = list(c(c_d15_isl1, c_d12_gfp), c(c_d15_gfp, c_d12_isl1)),
    alpha = 0.05)
  write_result(res_int, "RescueEffect_DayDiff")
} else {
  # One of the 4 levels is the reference level itself (coef name won't exist
  # for it) -- e.g. if D12_Ctr_GFP is the reference, all 4 needed coefs DO
  # exist since none of them IS D12_Ctr_GFP, so this branch should not
  # normally trigger. Flagging rather than silently skipping if it does.
  warning("Could not resolve all 4 coefficients for RescueEffect_DayDiff; check resultsNames(dds): ",
          paste(rn, collapse = ", "))
}

message("Done. Results in ", out_dir)

# NOTE: as with the limma script, RescueEffect_DayDiff has only n=2/group per
# day -- treat as exploratory, not confirmatory.
