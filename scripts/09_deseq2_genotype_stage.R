#!/usr/bin/env Rscript
# DESeq2, ISL1 rescue bulk RNA-seq -- genotype x stage only, rescue EXCLUDED.
#
# Restricted to Ctr-GFP and KO-GFP samples (8 of the 12) -- no KO-ISL1-rescue
# samples at all. With rescue out of the picture this IS a full, clean 2x2
# factorial (genotype: Ctr/KO x stage: D12/D15), unlike the 3-group design in
# script 07. Purpose: identify genes regulated BY ISL1 (KO vs Ctr = genotype
# effect) at each stage, for direct integration with stage-matched CUT&RUN
# peaks in script 10 -- KO-up genes are ISL1-repressed targets, KO-down genes
# are ISL1-activated targets.
#
# Input:  results/rsem/*.genes.results (same RSEM output as scripts 06/07,
#         strandedness-corrected)
# Output: results/DE_genotype_stage/{KOeffect_D12,KOeffect_D15,
#         StageEffect_Ctr,StageEffect_KO}.tsv
#
# Run in `r_env`. Usage: Rscript scripts/09_deseq2_genotype_stage.R

suppressPackageStartupMessages(library(DESeq2))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
repo_dir     <- normalizePath(file.path(get_script_dir(), ".."))
sample_sheet <- file.path(repo_dir, "scripts", "sample_sheet.csv")
rsem_dir     <- file.path(repo_dir, "results", "rsem")
out_dir      <- file.path(repo_dir, "results", "DE_genotype_stage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load metadata, restrict to Ctr_GFP / KO_GFP only (no rescue) ----------
meta <- read.csv(sample_sheet, stringsAsFactors = FALSE)
meta <- meta[meta$group3 %in% c("Ctr_GFP", "KO_GFP"), ]
meta$genotype <- ifelse(meta$group3 == "Ctr_GFP", "Ctr", "KO")
meta$group <- factor(paste(meta$day, meta$genotype, sep = "_"))
rownames(meta) <- meta$sample
message(sprintf("Using %d samples (rescue excluded):", nrow(meta)))
print(table(meta$day, meta$genotype))

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
counts <- round(counts)

# --- DESeqDataSet, filter, fit ---------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
message(sprintf("Testing %d genes after low-count filtering", nrow(dds)))
dds <- DESeq(dds)

pdf(file.path(out_dir, "pca_vst.pdf"), width = 6, height = 5)
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
plot(pcaData$PC1, pcaData$PC2, col = as.integer(factor(meta$genotype)) + 1, pch = 16,
     xlab = "PC1", ylab = "PC2", main = "PCA on VST counts (color = genotype)")
text(pcaData$PC1, pcaData$PC2, labels = meta$group, pos = 3, cex = 0.6)
dev.off()

# --- Contrasts: genotype effect (KOeffect) and stage effect ----------------
# All four are simple two-level differences within the single `group` factor
# -- results(contrast=c("group",A,B)) is correct regardless of reference level.
pairwise <- list(
  KOeffect_D12    = c("D12_KO", "D12_Ctr"),   # ISL1's regulatory targets at D12
  KOeffect_D15    = c("D15_KO", "D15_Ctr"),   # ISL1's regulatory targets at D15
  StageEffect_Ctr = c("D15_Ctr", "D12_Ctr"),  # normal D12->D15 progression, ISL1 intact
  StageEffect_KO  = c("D15_KO", "D12_KO")     # D12->D15 progression without ISL1
)

write_result <- function(res, cn) {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[order(res_df$pvalue), ]
  write.table(res_df, file.path(out_dir, paste0(cn, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  sig <- !is.na(res_df$padj) & res_df$padj < 0.05
  message(sprintf("%-16s  %4d genes FDR<0.05  (of %d tested)", cn, sum(sig, na.rm = TRUE), nrow(res_df)))
}

for (cn in names(pairwise)) {
  lvls <- pairwise[[cn]]
  res <- results(dds, contrast = c("group", lvls[1], lvls[2]), alpha = 0.05)
  write_result(res, cn)
}

message("Done. Results in ", out_dir)
