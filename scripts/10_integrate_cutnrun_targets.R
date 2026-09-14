#!/usr/bin/env Rscript
# Integrate stage-matched, endogenous (non-rescue) ISL1 CUT&RUN peaks with the
# genotype-effect (KOeffect) bulk RNA-seq DE results from script 09, to call
# direct ISL1 target genes -- activated or repressed -- at D12 and D15.
#
# CUT&RUN input: HOMER-annotated peaks for the Ctrl-GFP (endogenous, no
# rescue plasmid) ISL1 CUT&RUN samples, already computed in the companion
# CUT&RUN project:
#   /pl/active/ZhuLab_BDC/Han_analysis/ISL1_cut_N_run/analysis/d12_vs_d15_ctrl/
#     D12-Ctr-GFP-3d-1_ISL1.annotation.txt
#     D15-Ctr-GFP-3d-1_ISL1.annotation.txt
# These are genuine endogenous ISL1 binding (Ctrl cells, no OE-ISL1 plasmid)
# -- the rescue-condition CUT&RUN samples are deliberately not used here.
#
# Target-gene definition: HOMER's nearest-gene assignment (Nearest Ensembl),
# restricted to |Distance to TSS| <= 10kb (promoter-proximal direct binding;
# adjust DIST_CUTOFF below if a different window is wanted -- e.g. this does
# NOT use the project's ABC-model distal enhancer-gene links, only proximal
# nearest-TSS assignment).
#
# Direction logic: genotype-effect log2FC is KO vs Ctr (script 09). ISL1 is
# absent in KO, so:
#   log2FC > 0 in KO (gene goes UP when ISL1 is lost)  -> ISL1 normally REPRESSES it
#   log2FC < 0 in KO (gene goes DOWN when ISL1 is lost) -> ISL1 normally ACTIVATES it
# A "direct target" here = has a proximal CUT&RUN peak AND is significantly
# DE (FDR<0.05) in that stage's KOeffect contrast.
#
# Output: results/cutnrun_integration/D{12,15}_ISL1_{activated,repressed}_targets.tsv
#         results/cutnrun_integration/summary.tsv
#
# Run in `r_env`. Usage: Rscript scripts/10_integrate_cutnrun_targets.R

suppressPackageStartupMessages(library(dplyr))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
repo_dir <- normalizePath(file.path(get_script_dir(), ".."))
de_dir   <- file.path(repo_dir, "results", "DE_genotype_stage")
out_dir  <- file.path(repo_dir, "results", "cutnrun_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CUTRUN_DIR <- "/pl/active/ZhuLab_BDC/Han_analysis/ISL1_cut_N_run/analysis/d12_vs_d15_ctrl"
DIST_CUTOFF <- 10000  # bp from TSS, promoter-proximal direct-binding cutoff

peak_files <- list(
  D12 = file.path(CUTRUN_DIR, "D12-Ctr-GFP-3d-1_ISL1.annotation.txt"),
  D15 = file.path(CUTRUN_DIR, "D15-Ctr-GFP-3d-1_ISL1.annotation.txt")
)
de_files <- list(
  D12 = file.path(de_dir, "KOeffect_D12.tsv"),
  D15 = file.path(de_dir, "KOeffect_D15.tsv")
)

# HOMER annotatePeaks.pl output: read the real header (row 1) but the first
# column name is a long cmd-string -- rename positionally instead of relying
# on exact text.
read_homer_targets <- function(f, dist_cutoff) {
  d <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  names(d)[1] <- "PeakID"
  names(d)[names(d) == "Distance to TSS"] <- "DistanceToTSS"
  names(d)[names(d) == "Nearest Ensembl"] <- "gene_id"
  names(d)[names(d) == "Gene Name"] <- "gene_name"
  d <- d[!is.na(d$DistanceToTSS) & abs(d$DistanceToTSS) <= dist_cutoff & d$gene_id != "", ]
  # a gene can have multiple proximal peaks -- record the count, keep the
  # strongest (max Peak Score) peak's stats per gene
  d <- d %>% group_by(gene_id) %>% mutate(n_proximal_peaks = n()) %>% ungroup()
  d %>%
    group_by(gene_id) %>%
    slice_max(`Peak Score`, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(gene_id, gene_name, n_proximal_peaks,
              best_peak_score = `Peak Score`, best_peak_distance_to_tss = DistanceToTSS,
              best_peak_annotation = Annotation)
}

summary_rows <- list()

for (day in c("D12", "D15")) {
  message("=== ", day, " ===")
  targets <- read_homer_targets(peak_files[[day]], DIST_CUTOFF)
  message(sprintf("%d genes with a proximal (<=%dbp) ISL1 CUT&RUN peak", nrow(targets), DIST_CUTOFF))

  de <- read.delim(de_files[[day]], stringsAsFactors = FALSE)
  de_sig <- de[!is.na(de$padj) & de$padj < 0.05, ]
  message(sprintf("%d significant KOeffect genes (FDR<0.05)", nrow(de_sig)))

  merged <- merge(targets, de_sig, by = "gene_id")
  merged$direction <- ifelse(merged$log2FoldChange > 0, "ISL1_repressed", "ISL1_activated")
  merged <- merged[order(merged$padj), ]
  message(sprintf("%d direct ISL1 targets (peak + significant DE): %d activated, %d repressed",
                   nrow(merged), sum(merged$direction == "ISL1_activated"), sum(merged$direction == "ISL1_repressed")))

  cols <- c("gene_id", "gene_name", "log2FoldChange", "padj", "direction",
            "best_peak_distance_to_tss", "best_peak_score", "best_peak_annotation")
  write.table(merged[merged$direction == "ISL1_activated", cols],
              file.path(out_dir, paste0(day, "_ISL1_activated_targets.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(merged[merged$direction == "ISL1_repressed", cols],
              file.path(out_dir, paste0(day, "_ISL1_repressed_targets.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summary_rows[[day]] <- data.frame(
    day = day,
    n_genes_with_proximal_peak = nrow(targets),
    n_sig_KOeffect_genes = nrow(de_sig),
    n_direct_targets = nrow(merged),
    n_activated = sum(merged$direction == "ISL1_activated"),
    n_repressed = sum(merged$direction == "ISL1_repressed")
  )
}

summary_tbl <- do.call(rbind, summary_rows)
write.table(summary_tbl, file.path(out_dir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_tbl)
message("Done. Results in ", out_dir)
