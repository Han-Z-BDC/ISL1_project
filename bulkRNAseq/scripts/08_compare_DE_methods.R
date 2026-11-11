#!/usr/bin/env Rscript
# Compare limma-voom (scripts/06) vs DESeq2 (scripts/07) results, per contrast.
#
# For each contrast run by both methods: merge on gene_id, report the
# correlation of the two logFC estimates, a 2-way Venn of genes called
# significant (FDR<0.05) by each method, and how many of each method's hits
# are confirmed by the other. Large disagreement on a contrast is a flag to
# look at that contrast's PCA/volcano plots more closely, not necessarily a
# sign either tool is "wrong" -- with n=2/group, individual gene calls near
# the significance boundary are expected to be sensitive to modeling choices.
#
# Run in `r_env` (needs limma for vennDiagram + DESeq2/limma both already
# present there after scripts/06 and scripts/07's setup).
# Usage: Rscript scripts/08_compare_DE_methods.R

suppressPackageStartupMessages(library(limma))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
repo_dir    <- normalizePath(file.path(get_script_dir(), ".."))
limma_dir   <- file.path(repo_dir, "results", "DE")
deseq2_dir  <- file.path(repo_dir, "results", "DE_DESeq2")
out_dir     <- file.path(repo_dir, "results", "DE_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrasts <- c("KOeffect_D12", "KOeffect_D15",
               "RescueEffect_D12", "RescueEffect_D15",
               "RescueVsCtrl_D12", "RescueVsCtrl_D15",
               "RescueEffect_DayDiff")

summary_rows <- list()

for (cn in contrasts) {
  f_limma  <- file.path(limma_dir,  paste0(cn, ".tsv"))
  f_deseq2 <- file.path(deseq2_dir, paste0(cn, ".tsv"))
  if (!file.exists(f_limma) || !file.exists(f_deseq2)) {
    message("Skipping ", cn, " (missing output from one method -- run scripts 06 and 07 first)")
    next
  }

  lm  <- read.delim(f_limma,  stringsAsFactors = FALSE)   # cols: logFC, adj.P.Val, gene_id, ...
  ds  <- read.delim(f_deseq2, stringsAsFactors = FALSE)   # cols: log2FoldChange, padj, gene_id, ...

  m <- merge(lm[, c("gene_id", "logFC", "adj.P.Val")],
             ds[, c("gene_id", "log2FoldChange", "padj")],
             by = "gene_id")
  m$padj <- ifelse(is.na(m$padj), 1, m$padj)  # DESeq2 NAs padj for low-count/outlier genes

  r <- cor(m$logFC, m$log2FoldChange, use = "complete.obs")

  sig_limma  <- m$adj.P.Val < 0.05
  sig_deseq2 <- m$padj      < 0.05
  both       <- sig_limma & sig_deseq2

  # 2-way Venn via limma::vennDiagram, reusing limma (already a dependency)
  ind <- cbind(limma = sig_limma, DESeq2 = sig_deseq2)
  pdf(file.path(out_dir, paste0("venn_", cn, ".pdf")), width = 4, height = 4)
  vennDiagram(vennCounts(ind), main = cn)
  dev.off()

  write.table(m, file.path(out_dir, paste0("merged_", cn, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summary_rows[[cn]] <- data.frame(
    contrast       = cn,
    logFC_cor      = round(r, 3),
    n_limma_sig    = sum(sig_limma),
    n_deseq2_sig   = sum(sig_deseq2),
    n_both_sig     = sum(both),
    frac_limma_confirmed  = round(sum(both) / max(sum(sig_limma), 1), 3),
    frac_deseq2_confirmed = round(sum(both) / max(sum(sig_deseq2), 1), 3)
  )

  message(sprintf("%-20s  logFC r=%.3f  limma=%d DESeq2=%d both=%d",
                   cn, r, sum(sig_limma), sum(sig_deseq2), sum(both)))
}

summary_tbl <- do.call(rbind, summary_rows)
write.table(summary_tbl, file.path(out_dir, "method_comparison_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Done. Summary + per-contrast merged tables + Venns in ", out_dir)
