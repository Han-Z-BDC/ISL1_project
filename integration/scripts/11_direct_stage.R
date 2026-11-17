#!/usr/bin/env Rscript
# Step 11: ISL1 direct targets with a stage-specific KO effect, clustered by the timing of the effect; run once for
# the direct-activated and once for the direct-repressed genes.
#   Genes           direct-<label> in step 02: activated = lower in KO (padj < 0.05) on >= 1 of D15 / D18 / D22, never
#                   higher; repressed = higher in KO, never lower; >= 1 ISL1 site in the GREAT domain; ISL1 excluded.
#   Stage-specific  LRT_interaction padj < 0.05 (the KO effect differs between stages). The others: KO effect not
#                   shown to change between stages. The per-day significance pattern (sig_days) is descriptive.
#   Profile         shrunk KO log2FC at D15, D18, D22, divided by the gene's largest |log2FC| (timing, not size).
#   Clusters        k-means on the profiles of the stage-specific genes, k = 2-8 (seed 1, 25 starts, 100 iterations);
#                   clusters numbered by size within each k (ties by the order of the first gene). The user picks k.
#   Also            max_effect_day: the day with the largest shrunk KO log2FC in the label's direction.
#   group           the k = 4 cluster (user's choice) for stage-specific genes, "constant" for the others (KO effect
#                   not shown to change between stages); all_days = significant on D15, D18 and D22 (descriptive).
# usage: Rscript 11_direct_stage.R <integration results dir> <activated|repressed> <out_dir>
# Output in <out_dir>: genes.tsv (one row per gene of the label; k2-k8 for the stage-specific ones, all_days, group),
#   elbow.tsv, centers.tsv (per k and cluster: genes, mean scaled profile, median shrunk log2FC), patterns.tsv,
#   groups.tsv (the 4 clusters and the constant group: genes, median shrunk log2FC, share significant on all days)

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3, args[2] %in% c("activated", "repressed"))
res <- args[1]; lab <- args[2]; od <- args[3]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
days <- c("D15", "D18", "D22"); K <- 2:8

g <- fread(file.path(res, "02_direct_targets", "gene_table.tsv"), sep = "\t", na.strings = "NA", quote = "")
g <- g[direct_label == paste0("direct-", lab)]
g[, stage_specific := !is.na(padj_interaction) & padj_interaction < 0.05]
lfc <- as.matrix(g[, paste0("ko_lfc_", days), with = FALSE])
stopifnot(!anyNA(lfc))
mx <- apply(abs(lfc), 1, max)
prof <- lfc / mx
for (j in seq_along(days)) g[, (paste0("scaled_", days[j])) := round(prof[, j], 4)]
g[, max_effect_day := days[apply(lfc, 1, if (lab == "activated") which.min else which.max)]]
setorder(g, gene_id)   # fixed input order for k-means
prof <- as.matrix(g[, paste0("scaled_", days), with = FALSE])

run_km <- function(idx, pre) {
  m <- prof[idx, , drop = FALSE]
  elbow <- list(); centers <- list()
  for (k in K) {
    set.seed(1)
    f <- kmeans(m, centers = k, nstart = 25, iter.max = 100)
    sz <- tabulate(f$cluster, k)
    first <- sapply(seq_len(k), function(c) min(which(f$cluster == c)))
    ord <- order(-sz, first)
    cl <- match(f$cluster, ord)
    g[, (paste0(pre, k)) := NA_integer_]
    g[idx, (paste0(pre, k)) := cl]
    elbow[[length(elbow) + 1]] <- data.table(k = k, tot_withinss = round(f$tot.withinss, 3),
                                             pct_between = round(100 * f$betweenss / f$totss, 1))
    for (c in seq_len(k)) {
      i <- idx[cl == c]
      centers[[length(centers) + 1]] <- data.table(k = k, cluster = c, genes = length(i),
        scaled_D15 = round(mean(prof[i, 1]), 3), scaled_D18 = round(mean(prof[i, 2]), 3), scaled_D22 = round(mean(prof[i, 3]), 3),
        lfc_D15 = round(median(g$ko_lfc_D15[i]), 3), lfc_D18 = round(median(g$ko_lfc_D18[i]), 3),
        lfc_D22 = round(median(g$ko_lfc_D22[i]), 3))
    }
  }
  list(elbow = rbindlist(elbow), centers = rbindlist(centers))
}
g[, all_days := sig_days == "D15,D18,D22"]
ss <- which(g$stage_specific)
A <- run_km(ss, "k")
elbow <- A$elbow; centers <- A$centers
g[, group := fifelse(stage_specific, paste0("cluster", k4), "constant")]
grp <- g[, .(genes = .N, lfc_D15 = round(median(ko_lfc_D15), 3), lfc_D18 = round(median(ko_lfc_D18), 3),
             lfc_D22 = round(median(ko_lfc_D22), 3), all_days = sum(all_days),
             pct_all_days = round(100 * mean(all_days), 1)), by = group]
setorder(grp, group)

pat <- g[, .N, by = .(sig_days, stage_specific)]
pat <- dcast(pat, sig_days ~ stage_specific, value.var = "N", fill = 0)
setnames(pat, c("FALSE", "TRUE"), c("not_stage_specific", "stage_specific"))
pat[, sig_days := fifelse(sig_days == "", "none (interaction only)", sig_days)]
pat[, total := not_stage_specific + stage_specific]
setorder(pat, -total, sig_days)

cols <- c("gene_id", "gene_name", "gene_type", "stage_specific", "padj_interaction", "sig_days",
          paste0("ko_lfc_", days), paste0("ko_padj_", days), paste0("scaled_", days), "max_effect_day",
          paste0("k", K), "all_days", "group", "n_sites", "bound_near50")
out <- g[, cols, with = FALSE]
setorder(out, -stage_specific, gene_name, gene_id)
fwrite(out, file.path(od, "genes.tsv"), sep = "\t", na = "NA", quote = FALSE)
fwrite(elbow, file.path(od, "elbow.tsv"), sep = "\t", quote = FALSE)
fwrite(centers, file.path(od, "centers.tsv"), sep = "\t", quote = FALSE)
fwrite(pat, file.path(od, "patterns.tsv"), sep = "\t", quote = FALSE)
fwrite(grp, file.path(od, "groups.tsv"), sep = "\t", quote = FALSE)

message(sprintf("%d direct-%s genes, %d stage-specific", nrow(g), lab, length(ss)))
print(pat); print(elbow); print(centers)
print(g[ss, .N, by = max_effect_day][order(max_effect_day)])
print(grp)
