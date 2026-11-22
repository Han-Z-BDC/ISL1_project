# Step 7.1n: joint DESeq2 of the SC_EC1 / SC_EC2 pseudobulks of D22 and D40 (after 7.1m, which writes the
# pseudobulks). Request (user, 2026-11-22): a comprehensive comparison with the variables cell type (EC1 / EC2), stage
# (D22 / D40) and genotype (control / KO); mainly EC1 vs EC2, then control vs KO within each cell type.
#
# Pseudobulks: SoupX counts summed per sample x label (7.1m), SC_EC1 and SC_EC2 of the 8 samples (2 days x 2 genotypes
# x 2 replicates) = 16. Genes with >= 10 counts in >= 2 pseudobulks. Day = stage = differentiation (D22 diff 37, D40
# diff 34, separate 10x runs), so day effects include batch; the two days are read as two independent experiments.
#   Model A, cell type (EC2 vs EC1), paired within sample (DESeq2 vignette design for a within-individual condition
#     with individuals nested in groups): ~ group + group:rep + group:celltype, group = day_genotype. Contrasts:
#     EC2 vs EC1 within each group (D22_ctr, D22_KO, D40_ctr, D40_KO); average over the 4 groups; interactions =
#     difference of the EC2 vs EC1 effect KO - control (averaged over days) and D40 - D22 (averaged over genotypes).
#   Model B, genotype within cell type (KO vs control, unpaired: genotype differs between samples): ~ 0 + cell,
#     cell = celltype_genotype_day (8 groups). Contrasts: KO vs control in EC1 at D22, D40, average of the days; the
#     same in EC2; interactions = KO effect in EC2 - KO effect in EC1 (average of the days) and KO effect at D40 - D22
#     (average of the cell types).
#   Each contrast: Wald test (numeric contrast), BH padj, MLE log2FC and ashr-shrunk log2FC; call = padj < 0.05 and
#   |log2FC (MLE)| >= 0.5 (as 7.1m / 7.3l). ISL1 is flagged (the KO allele is still transcribed, higher in KO).
#   Gene classes (model B, average-of-days calls): KO shared (EC1 and EC2, same direction), EC1 only, EC2 only;
#   per day calls kept. Core EC1 / EC2 genes (model A): average call and the same sign in all 4 groups.
# Writes ${OUT_DIR}/day_analysis/ec_deseq2/ (reruns overwrite):
#   samples.csv, normalized_counts.csv, vst.csv
#   results_A.csv / results_B.csv   one row per gene: per contrast log2FC, log2FC_shrunk, pvalue, padj, call; TF flag
#   genes_celltype.csv              core EC2-high / EC1-high genes (model A)
#   genes_KO.csv                    KO classes (model B)
#   summary.csv                     calls per contrast (up / down)
# Usage: Rscript 7.1n_ec_deseq2.R   (DESeq2 from the DiffBind env; see 7.1n_ec_deseq2.sb)
suppressPackageStartupMessages({ library(DESeq2); library(data.table) })
da <- file.path(Sys.getenv("OUT_DIR"), "day_analysis")
src <- file.path(da, "ec_similarity"); od <- file.path(da, "ec_deseq2"); dir.create(od, showWarnings = FALSE)
DAYS <- c("D22", "D40"); EC <- c("SC_EC1", "SC_EC2"); LFC <- 0.5; ALPHA <- 0.05
TF <- readLines(file.path(Sys.getenv("REPO_DIR"), "study_config", "annotation", "TF_names_v_1.01.txt"))

# ---- pseudobulks
cnt <- list(); smp <- list()
for (d in DAYS) {
  x <- fread(cmd = sprintf("zcat %s", file.path(src, sprintf("pb_counts_%s.csv.gz", d))))   # no R.utils in this env
  m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]
  s <- fread(file.path(src, sprintf("pb_samples_%s.csv", d)))[match(colnames(m), pb)]
  k <- s$label %in% EC
  m <- m[, k]; s <- s[k]
  s[, `:=`(day = d, rep = sub(".*_", "", sample), celltype = sub("SC_", "", label))]
  s[, id := sprintf("%s_%s_%s_%s", celltype, genotype, d, rep)]
  colnames(m) <- s$id
  cnt[[d]] <- m; smp[[d]] <- s
}
g <- Reduce(intersect, lapply(cnt, rownames))
m <- do.call(cbind, lapply(cnt, function(x) x[g, ]))
s <- rbindlist(smp)
stopifnot(identical(colnames(m), s$id), nrow(s) == 16)
s[, `:=`(genotype = factor(genotype, c("ctr", "KO")), day = factor(day, DAYS), rep = factor(rep),
         celltype = factor(celltype, c("EC1", "EC2")))]
s[, group := factor(paste(day, genotype, sep = "_"), c("D22_ctr", "D22_KO", "D40_ctr", "D40_KO"))]
s[, cell := factor(paste(celltype, genotype, day, sep = "_"))]
fwrite(s[, .(id, celltype, genotype, day, rep, sample, cells)], file.path(od, "samples.csv"))
print(s[, .(id, cells)])
keep <- rowSums(m >= 10) >= 2
m <- m[keep, ]
cat("genes kept:", nrow(m), "of", length(keep), "\n")
cd <- as.data.frame(s); rownames(cd) <- s$id

run_contrasts <- function(dds, cons) {
  rn <- resultsNames(dds)
  out <- data.table(gene = rownames(dds))
  for (cn in names(cons)) {
    w <- cons[[cn]]; v <- setNames(numeric(length(rn)), rn); v[names(w)] <- w
    r <- results(dds, contrast = v, alpha = ALPHA)
    sh <- lfcShrink(dds, contrast = v, res = r, type = "ashr", quiet = TRUE)
    out[, (paste0(cn, "__log2FC")) := r$log2FoldChange]
    out[, (paste0(cn, "__log2FC_shrunk")) := sh$log2FoldChange]
    out[, (paste0(cn, "__pvalue")) := r$pvalue]
    out[, (paste0(cn, "__padj")) := r$padj]
    out[, (paste0(cn, "__call")) := !is.na(r$padj) & r$padj < ALPHA & abs(r$log2FoldChange) >= LFC]
    cat(sprintf("  %-28s called %5d (up %d / down %d)\n", cn, sum(out[[paste0(cn, "__call")]]),
                sum(out[[paste0(cn, "__call")]] & r$log2FoldChange > 0, na.rm = TRUE),
                sum(out[[paste0(cn, "__call")]] & r$log2FoldChange < 0, na.rm = TRUE)))
  }
  out[, TF := gene %in% TF][, ISL1_flag := gene == "ISL1"]
  out
}
summ <- function(res, model) {
  cs <- sub("__call$", "", grep("__call$", names(res), value = TRUE))
  rbindlist(lapply(cs, function(cn) data.table(model = model, contrast = cn, called = sum(res[[paste0(cn, "__call")]]),
    up = sum(res[[paste0(cn, "__call")]] & res[[paste0(cn, "__log2FC")]] > 0),
    down = sum(res[[paste0(cn, "__call")]] & res[[paste0(cn, "__log2FC")]] < 0))))
}

# ---- model A: EC2 vs EC1, paired within sample
mmA <- model.matrix(~ group + group:rep + group:celltype, cd)
mmA <- mmA[, colSums(mmA != 0) > 0]
stopifnot(qr(mmA)$rank == ncol(mmA))
ddsA <- DESeq(DESeqDataSetFromMatrix(m, cd, design = mmA), quiet = TRUE)
rnA <- resultsNames(ddsA); ct <- grep("celltypeEC2$", rnA, value = TRUE)
cat("model A coefficients for EC2 vs EC1:", paste(ct, collapse = ", "), "\n")
nm <- function(gp) grep(sprintf("^group%s.celltypeEC2$", gp), ct, value = TRUE)   # ":" may become "."
consA <- list(
  EC2_vs_EC1_D22_ctr = setNames(1, nm("D22_ctr")), EC2_vs_EC1_D22_KO = setNames(1, nm("D22_KO")),
  EC2_vs_EC1_D40_ctr = setNames(1, nm("D40_ctr")), EC2_vs_EC1_D40_KO = setNames(1, nm("D40_KO")),
  EC2_vs_EC1_average = setNames(rep(0.25, 4), ct),
  EC2_vs_EC1_diff_KO_minus_ctr = setNames(c(-0.5, 0.5, -0.5, 0.5), c(nm("D22_ctr"), nm("D22_KO"), nm("D40_ctr"), nm("D40_KO"))),
  EC2_vs_EC1_diff_D40_minus_D22 = setNames(c(-0.5, -0.5, 0.5, 0.5), c(nm("D22_ctr"), nm("D22_KO"), nm("D40_ctr"), nm("D40_KO"))))
cat("model A\n"); resA <- run_contrasts(ddsA, consA)
fwrite(resA, file.path(od, "results_A.csv"))

# ---- model B: KO vs control within each cell type
ddsB <- DESeq(DESeqDataSetFromMatrix(m, cd, design = ~ 0 + cell), quiet = TRUE)
cc <- function(x) paste0("cell", x)
consB <- list(
  KO_vs_ctr_EC1_D22 = setNames(c(1, -1), cc(c("EC1_KO_D22", "EC1_ctr_D22"))),
  KO_vs_ctr_EC1_D40 = setNames(c(1, -1), cc(c("EC1_KO_D40", "EC1_ctr_D40"))),
  KO_vs_ctr_EC1_average = setNames(c(0.5, -0.5, 0.5, -0.5), cc(c("EC1_KO_D22", "EC1_ctr_D22", "EC1_KO_D40", "EC1_ctr_D40"))),
  KO_vs_ctr_EC2_D22 = setNames(c(1, -1), cc(c("EC2_KO_D22", "EC2_ctr_D22"))),
  KO_vs_ctr_EC2_D40 = setNames(c(1, -1), cc(c("EC2_KO_D40", "EC2_ctr_D40"))),
  KO_vs_ctr_EC2_average = setNames(c(0.5, -0.5, 0.5, -0.5), cc(c("EC2_KO_D22", "EC2_ctr_D22", "EC2_KO_D40", "EC2_ctr_D40"))),
  KO_effect_diff_EC2_minus_EC1 = setNames(c(0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, 0.5),
    cc(c("EC2_KO_D22", "EC2_ctr_D22", "EC2_KO_D40", "EC2_ctr_D40", "EC1_KO_D22", "EC1_ctr_D22", "EC1_KO_D40", "EC1_ctr_D40"))),
  KO_effect_diff_D40_minus_D22 = setNames(c(0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, 0.5),
    cc(c("EC1_KO_D40", "EC1_ctr_D40", "EC2_KO_D40", "EC2_ctr_D40", "EC1_KO_D22", "EC1_ctr_D22", "EC2_KO_D22", "EC2_ctr_D22"))))
stopifnot(all(unlist(lapply(consB, names)) %in% resultsNames(ddsB)))
cat("model B\n"); resB <- run_contrasts(ddsB, consB)
fwrite(resB, file.path(od, "results_B.csv"))

# ---- gene classes
grp4 <- c("EC2_vs_EC1_D22_ctr", "EC2_vs_EC1_D22_KO", "EC2_vs_EC1_D40_ctr", "EC2_vs_EC1_D40_KO")
sg <- sign(as.matrix(resA[, paste0(grp4, "__log2FC"), with = FALSE]))
resA[, same_sign_4 := apply(sg, 1, function(z) !anyNA(z) & length(unique(z)) == 1)]
gA <- resA[EC2_vs_EC1_average__call == TRUE & same_sign_4 == TRUE,
           .(gene, class = fifelse(EC2_vs_EC1_average__log2FC > 0, "EC2-high", "EC1-high"),
             log2FC = EC2_vs_EC1_average__log2FC, log2FC_shrunk = EC2_vs_EC1_average__log2FC_shrunk,
             padj = EC2_vs_EC1_average__padj, genotype_dependent = EC2_vs_EC1_diff_KO_minus_ctr__call,
             day_dependent = EC2_vs_EC1_diff_D40_minus_D22__call, TF)][order(class, -abs(log2FC))]
fwrite(gA, file.path(od, "genes_celltype.csv"))
cat("core EC1 / EC2 genes:\n"); print(gA[, .N, by = class])

b1 <- resB$KO_vs_ctr_EC1_average__call; b2 <- resB$KO_vs_ctr_EC2_average__call
f1 <- resB$KO_vs_ctr_EC1_average__log2FC; f2 <- resB$KO_vs_ctr_EC2_average__log2FC
resB[, KO_class := fcase(b1 & b2 & sign(f1) == sign(f2), "shared",
                         b1 & b2, "both, opposite direction",
                         b1, "EC1 only", b2, "EC2 only", default = "none")]
resB[, KO_direction := fifelse(KO_class == "shared", fifelse(f1 > 0, "up in KO", "down in KO"),
                       fifelse(KO_class == "EC1 only", fifelse(f1 > 0, "up in KO", "down in KO"),
                       fifelse(KO_class == "EC2 only", fifelse(f2 > 0, "up in KO", "down in KO"), NA_character_)))]
resB[, both_days_EC1 := KO_vs_ctr_EC1_D22__call & KO_vs_ctr_EC1_D40__call &
                        sign(KO_vs_ctr_EC1_D22__log2FC) == sign(KO_vs_ctr_EC1_D40__log2FC)]
resB[, both_days_EC2 := KO_vs_ctr_EC2_D22__call & KO_vs_ctr_EC2_D40__call &
                        sign(KO_vs_ctr_EC2_D22__log2FC) == sign(KO_vs_ctr_EC2_D40__log2FC)]
gB <- resB[KO_class != "none", .(gene, KO_class, KO_direction,
  EC1_log2FC = KO_vs_ctr_EC1_average__log2FC, EC1_padj = KO_vs_ctr_EC1_average__padj,
  EC2_log2FC = KO_vs_ctr_EC2_average__log2FC, EC2_padj = KO_vs_ctr_EC2_average__padj,
  both_days_EC1, both_days_EC2, celltype_dependent = KO_effect_diff_EC2_minus_EC1__call,
  day_dependent = KO_effect_diff_D40_minus_D22__call, TF, ISL1_flag)][order(KO_class, KO_direction, EC1_padj)]
fwrite(gB, file.path(od, "genes_KO.csv"))
cat("KO classes:\n"); print(resB[, .N, by = .(KO_class, KO_direction)][order(KO_class)])
fwrite(rbind(summ(resA, "A"), summ(resB, "B")), file.path(od, "summary.csv"))

# ---- normalized counts / VST for plots
fwrite(data.table(gene = rownames(ddsB), counts(ddsB, normalized = TRUE)), file.path(od, "normalized_counts.csv"))
fwrite(data.table(gene = rownames(ddsB), assay(vst(ddsB, blind = FALSE))), file.path(od, "vst.csv"))
cat("wrote", od, "\n")
