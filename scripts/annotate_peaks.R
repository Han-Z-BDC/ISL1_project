#!/usr/bin/env Rscript
# Genomic annotation and target genes for every final peak set (06), with ChIPseeker.
#   peaks: <final_peaks>/<set>/<group>.narrowPeak (column 8 = -log10 p from IDR, column 10 = summit)
#   genes: TxDb.Hsapiens.UCSC.hg38.knownGene, symbols from org.Hs.eg.db; promoter = TSS +/- 3 kb
# Writes to <out_dir>:
#   <set>/<group>.annotated.tsv     one row per peak: annotation, nearest gene, distance to TSS
#   <set>/<group>.target_genes.tsv  one row per nearest gene: peaks, best -log10 p, closest distance,
#                                   whether any peak is in its promoter
#   annotation_summary.tsv          peaks per collapsed category, per set x group
#   tss_distance_summary.tsv        peaks per distance-to-TSS bin, per set x group
#   islet_panel.tsv                 distance from each panel gene's TSS to the nearest summit
# usage: Rscript annotate_peaks.R <final_peaks_dir> <out_dir>

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
fp_dir <- args[1]; out_dir <- args[2]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
canon <- paste0("chr", c(1:22, "X", "Y"))

PANEL <- c("INS", "GCG", "SST", "PPY", "GHRL", "IAPP", "CHGA", "PDX1", "NKX6-1", "NKX2-2",
           "NEUROD1", "NEUROG3", "PAX6", "MAFA", "MAFB", "ARX", "ISL1", "LMX1B", "FOXA2", "INSM1")

collapse_cat <- function(a) {
  fcase(startsWith(a, "Promoter (<=1kb)"), "Promoter (<=1 kb)",
        startsWith(a, "Promoter"), "Promoter (1-3 kb)",
        startsWith(a, "5' UTR") | startsWith(a, "3' UTR") | startsWith(a, "Exon"), "UTR / exon",
        startsWith(a, "Intron"), "Intron",
        startsWith(a, "Downstream"), "Downstream (<=300 bp)",
        startsWith(a, "Distal"), "Distal intergenic",
        default = "Other")
}
BREAKS <- c(-Inf, -1e5, -1e4, -1e3, 0, 1e3, 1e4, 1e5, Inf)
BIN_LABELS <- c("< -100 kb", "-100 to -10 kb", "-10 to -1 kb", "-1 to 0 kb",
                "0 to 1 kb", "1 to 10 kb", "10 to 100 kb", "> 100 kb")

# Panel gene TSSs (strand-aware), canonical chromosomes, one per gene (longest span).
g <- genes(txdb, single.strand.genes.only = TRUE)
g <- g[seqnames(g) %in% canon]
sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(g), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
g$symbol <- unname(sym)
g <- g[!is.na(g$symbol) & g$symbol %in% PANEL]
g <- g[order(-width(g))]; g <- g[!duplicated(g$symbol)]
panel_tss <- resize(g, 1, fix = "start")

ann_rows <- list(); tss_rows <- list(); panel_rows <- list()
for (set in c("conservative", "optimal")) {
  dir.create(file.path(out_dir, set), recursive = TRUE, showWarnings = FALSE)
  for (f in sort(Sys.glob(file.path(fp_dir, set, "*.narrowPeak")))) {
    group <- sub("\\.narrowPeak$", "", basename(f))
    np <- fread(f, header = FALSE, col.names = c("chr", "start", "end", "name", "score", "strand",
                                                  "signal", "neglog10p", "neglog10q", "summit"))
    # Annotate at the summit, so "nearest gene" and the TSS distance refer to the binding site.
    summit_gr <- GRanges(np$chr, IRanges(np$start + np$summit + 1, width = 1), name = np$name, neglog10p = np$neglog10p)
    a <- as.data.table(as.data.frame(annotatePeak(summit_gr, tssRegion = c(-3000, 3000), TxDb = txdb,
                                                  annoDb = "org.Hs.eg.db", verbose = FALSE)))
    a[, category := collapse_cat(annotation)]
    a <- merge(np[, .(name, chr, start, end)], a[, .(name, summit = start, neglog10p, annotation, category,
                                                     geneId, SYMBOL, GENENAME, distanceToTSS)], by = "name")
    setorder(a, -neglog10p)
    fwrite(a, file.path(out_dir, set, paste0(group, ".annotated.tsv")), sep = "\t")

    tg <- a[!is.na(SYMBOL), .(n_peaks = .N, best_neglog10p = max(neglog10p),
                              min_abs_dist_to_tss = min(abs(distanceToTSS)),
                              promoter_bound = any(abs(distanceToTSS) <= 3000)), by = SYMBOL]
    setorder(tg, -best_neglog10p)
    fwrite(tg, file.path(out_dir, set, paste0(group, ".target_genes.tsv")), sep = "\t")

    ann_rows[[length(ann_rows) + 1]] <- a[, .(n = .N), by = category][, `:=`(set = set, group = group)]
    a[, bin := cut(distanceToTSS, BREAKS, labels = BIN_LABELS, right = FALSE)]
    tss_rows[[length(tss_rows) + 1]] <- a[, .(n = .N), by = bin][, `:=`(set = set, group = group)]

    hit <- distanceToNearest(panel_tss, summit_gr, ignore.strand = TRUE)
    d <- rep(NA_real_, length(panel_tss)); d[queryHits(hit)] <- mcols(hit)$distance
    panel_rows[[length(panel_rows) + 1]] <- data.table(set = set, group = group, gene = panel_tss$symbol, dist_bp = d)
    cat(sprintf("%s %s: %d peaks, %d nearest genes, %d promoter-bound\n",
                set, group, nrow(a), nrow(tg), sum(tg$promoter_bound)))
  }
}
fwrite(rbindlist(ann_rows)[, .(set, group, category, n)], file.path(out_dir, "annotation_summary.tsv"), sep = "\t")
fwrite(rbindlist(tss_rows)[, .(set, group, bin, n)], file.path(out_dir, "tss_distance_summary.tsv"), sep = "\t")
fwrite(rbindlist(panel_rows), file.path(out_dir, "islet_panel.tsv"), sep = "\t")
missing <- setdiff(PANEL, panel_tss$symbol)
if (length(missing)) cat("panel genes not found in TxDb:", paste(missing, collapse = ", "), "\n")
