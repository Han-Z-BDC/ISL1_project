#!/usr/bin/env Rscript
# Differential expression, ISL1 rescue bulk RNA-seq.
#
# Design: NOT a full factorial -- Ctr cells only ever get the GFP construct,
# only KO cells get the ISL1-rescue construct, so genotype x construct isn't
# separately estimable. Instead we use one 3-level factor (group3: Ctr_GFP,
# KO_GFP, KO_ISL1) crossed with day (D12, D15), fit as a single no-intercept
# 6-level group factor, and pull out exactly the contrasts we care about.
#
# Input:  results/rsem/<sample>.genes.results  (RSEM gene-level output,
#         one flat file per sample, prefix = sample name from sample_sheet.csv)
# Output: results/DE/<contrast>.tsv, results/DE/volcano_<contrast>.pdf,
#         results/DE/pca_logCPM.pdf
#
# Run in the lab's `r_env` (has edgeR + limma already; DESeq2 is not
# installed there -- this uses limma-voom instead, which needs no new pkgs).
#
# Usage: Rscript scripts/06_differential_expression.R

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

# --- Paths ---------------------------------------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()  # fallback if sourced interactively from scripts/
}
repo_dir     <- normalizePath(file.path(get_script_dir(), ".."))
sample_sheet <- file.path(repo_dir, "scripts", "sample_sheet.csv")
rsem_dir     <- file.path(repo_dir, "results", "rsem")
out_dir      <- file.path(repo_dir, "results", "DE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load sample metadata + build the counts matrix -----------------------
meta <- read.csv(sample_sheet, stringsAsFactors = FALSE)
stopifnot(all(c("sample", "day", "group3") %in% names(meta)))

read_rsem_counts <- function(sample) {
  f <- file.path(rsem_dir, paste0(sample, ".genes.results"))
  d <- read.delim(f, stringsAsFactors = FALSE)
  setNames(d$expected_count, d$gene_id)
}

counts_list <- lapply(meta$sample, read_rsem_counts)
stopifnot(length(unique(sapply(counts_list, length))) == 1)  # all same gene set
counts <- do.call(cbind, counts_list)
colnames(counts) <- meta$sample
rownames(counts) <- names(counts_list[[1]])
counts <- round(counts)  # RSEM expected_count is fractional; edgeR/limma want integers

# --- Group factor: 6 levels = day x group3 --------------------------------
meta$group <- factor(paste(meta$day, meta$group3, sep = "_"))
group <- meta$group
message("Groups:\n", paste(capture.output(table(group)), collapse = "\n"))

# --- Build DGEList, filter, normalize -------------------------------------
dge <- DGEList(counts = counts, group = group)
keep <- filterByExpr(dge, group = group)
message(sprintf("Keeping %d / %d genes after filterByExpr", sum(keep), length(keep)))
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# --- PCA / MDS sanity check on logCPM, BEFORE trusting any DE result ------
logcpm <- cpm(dge, log = TRUE)
pdf(file.path(out_dir, "pca_logCPM.pdf"), width = 6, height = 5)
plotMDS(logcpm, labels = meta$group, col = as.integer(factor(meta$day)) + 1,
        main = "MDS on logCPM (color = day, label = group)")
dev.off()

# --- Design: no-intercept, one column per group ---------------------------
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

# --- Contrasts of biological interest --------------------------------------
# Adjust these level names if your day/group3 labels differ from D12/D15,
# Ctr_GFP/KO_GFP/KO_ISL1.
cm <- makeContrasts(
  KOeffect_D12          = D12_KO_GFP  - D12_Ctr_GFP,
  KOeffect_D15          = D15_KO_GFP  - D15_Ctr_GFP,
  RescueEffect_D12      = D12_KO_ISL1 - D12_KO_GFP,
  RescueEffect_D15      = D15_KO_ISL1 - D15_KO_GFP,
  RescueVsCtrl_D12      = D12_KO_ISL1 - D12_Ctr_GFP,
  RescueVsCtrl_D15      = D15_KO_ISL1 - D15_Ctr_GFP,
  RescueEffect_DayDiff  = (D15_KO_ISL1 - D15_KO_GFP) - (D12_KO_ISL1 - D12_KO_GFP),
  levels = design
)

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

# --- Write results + volcano plot per contrast ------------------------------
for (cn in colnames(cm)) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  tt$gene_id <- rownames(tt)
  write.table(tt, file.path(out_dir, paste0(cn, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(out_dir, paste0("volcano_", cn, ".pdf")), width = 5, height = 5)
  sig <- tt$adj.P.Val < 0.05
  plot(tt$logFC, -log10(tt$P.Value),
       col = ifelse(sig, "red", "grey60"), pch = 16, cex = 0.6,
       xlab = "log2 fold change", ylab = "-log10(p-value)", main = cn)
  abline(v = c(-1, 1), lty = 2, col = "grey40")
  dev.off()

  message(sprintf("%-22s  %4d genes FDR<0.05  (of %d tested)",
                   cn, sum(sig, na.rm = TRUE), nrow(tt)))
}

message("Done. Results in ", out_dir)

# NOTE on RescueEffect_DayDiff: only n=2/group per day, so treat this
# interaction contrast as exploratory, not confirmatory -- report the
# per-day RescueEffect_D12 / RescueEffect_D15 contrasts as the primary result
# and use this one only to flag candidates worth a qualitative look.
