# Step 7.2 (R part): edgeR quasi-likelihood test, KO vs ctr, on per-sample pseudobulk counts.
# usage: Rscript 7.2_pseudobulk_edger.R <outdir>   (reads pb_counts.csv, pb_meta.csv; writes edger_*.csv)
suppressPackageStartupMessages(library(edgeR))
out <- commandArgs(trailingOnly = TRUE)[1]
counts <- as.matrix(read.csv(file.path(out, "pb_counts.csv"), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(out, "pb_meta.csv"), row.names = 1)
meta <- meta[colnames(counts), , drop = FALSE]
meta$genotype <- factor(meta$genotype, levels = c("ctr", "KO"))

y <- DGEList(counts, group = meta$genotype)
keep <- filterByExpr(y, group = meta$genotype)
cat("genes kept by filterByExpr:", sum(keep), "of", length(keep), "\n")
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ genotype, data = meta)
y <- estimateDisp(y, design)
cat("common BCV:", round(sqrt(y$common.dispersion), 3), "\n")
fit <- glmQLFit(y, design, robust = TRUE)
res <- topTags(glmQLFTest(fit, coef = "genotypeKO"), n = Inf, sort.by = "PValue")$table
res <- data.frame(gene = rownames(res), res, check.names = FALSE)
write.csv(res, file.path(out, "edger_KO_vs_ctr.csv"), row.names = FALSE)
write.csv(cpm(y, log = TRUE, prior.count = 1), file.path(out, "edger_logCPM.csv"))
cat("FDR < 0.05:", sum(res$FDR < 0.05), " (up", sum(res$FDR < 0.05 & res$logFC > 0),
    ", down", sum(res$FDR < 0.05 & res$logFC < 0), ")\n")
