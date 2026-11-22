# Step 13: ISL1 target genes from TSS overlap and ABC enhancer-gene links (user, 2026-11-22).
#
# ISL1 sites   the macs3_ge2 consensus regions of D15, D18 and D22 (step 05), pooled and merged where they overlap
#              (bookended regions merged too); each site keeps the days it was called on and the days it is in the
#              3-of-3 set.
# TSS link     the site lies within +/- TSS_WIN bp of any annotated transcript start of a gene (GENES_GTF, the 10x
#              GRCh38 2020-A / GENCODE v32 annotation the scRNA-seq used, all transcripts, all gene types) -> that gene.
# ABC link     the site overlaps (>= 1 bp) an ABC element linked to a gene (ABC >= 0.02, own-promoter links included)
#              in any of the ABC_GROUPS of the diff-multiome ABC run (ENP2, ENP3, early_SC_alpha; PP2 left out by the
#              user) -> that gene. ABC genes are named through their Ensembl ID with the GTF names (the ABC gene list
#              uses its own symbols); IDs missing from the GTF keep the ABC symbol.
# A site can link to several genes; a gene is an ISL1 target when it has >= 1 TSS or ABC link.
# usage: Rscript abc_targets.R <consensus dir> <genes.gtf> <ABC results dir> <groups, comma-separated> <TSS window bp>
#        <out dir>
# Output in <out dir>:
#   sites.tsv    site, chr, start, end, days, days_3of3, score (max over days), n_genes, genes
#   links.tsv    one row per site x gene: source (TSS / ABC / TSS+ABC), TSS distance (nearest TSS of the gene, from
#                the site edge; 0 = overlaps), ABC: best score, groups, element classes, own promoter
#   genes.tsv    one row per linked gene: sites (all / TSS / ABC-only), days, TSS-bound, best ABC, groups
#   summary.tsv  counts; unmatched.tsv  ABC genes whose Ensembl ID is not in the GTF
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) == 6)
cdir <- a[1]; gtf_f <- a[2]; abc <- a[3]; groups <- strsplit(a[4], ",")[[1]]; WIN <- as.integer(a[5]); od <- a[6]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
DAYS <- c("D15", "D18", "D22")
canon <- function(d) d[grepl("^chr([0-9]+|X|Y)$", chr)]

# ---- ISL1 sites: pooled day regions, merged
rd_bed <- function(day, set) {
  d <- fread(file.path(cdir, sprintf("%s_Ctr_ISL1.%s.bed", day, set)), header = FALSE,
             col.names = c("chr", "start", "end", "name", "n_reps", "reps", "score", "n_batch3", "n_earlier"))
  d[, `:=`(day = day, set = set)]
}
reg <- rbindlist(lapply(DAYS, rd_bed, set = "macs3_ge2"))
r3 <- rbindlist(lapply(DAYS, rd_bed, set = "macs3_3of3"))
setorder(reg, chr, start)
reg[, grp := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = chr]   # bed half-open: bookended merge (start == end)
sites <- reg[, .(start = min(start), end = max(end), days = paste(intersect(DAYS, day), collapse = ","),
                 score = max(score)), by = .(chr, grp)]
sites <- canon(sites)[order(chr, start)]
sites[, site := sprintf("isl1_%05d", .I)][, grp := NULL]
# overlaps on 1-based closed coordinates (foverlaps counts shared end points; BED is 0-based half-open)
cl <- function(d) { d <- copy(d); d[, start := start + 1L]; d }
S1 <- cl(sites); setkey(S1, chr, start, end)
o3 <- foverlaps(cl(r3)[, .(chr, start, end, day)], S1, nomatch = NULL)
sites[o3[, .(d3 = paste(intersect(DAYS, day), collapse = ",")), by = site], days_3of3 := i.d3, on = "site"]
sites[is.na(days_3of3), days_3of3 := ""]
cat("ISL1 sites:", nrow(sites), "(from", nrow(reg), "day regions)\n"); print(sites[, .N, by = days][order(-N)])

# ---- TSSs (all transcripts) from the GTF
g <- fread(cmd = sprintf("awk -F'\\t' '$3==\"transcript\"' %s", gtf_f), header = FALSE, sep = "\t", quote = "")
g <- g[, .(chr = V1, tss = ifelse(V7 == "+", V4, V5), attr = V9)]
ex <- function(x, k) sub(sprintf('.*%s "([^"]+)".*', k), "\\1", x)
g[, `:=`(gene_id = ex(attr, "gene_id"), gene_name = ex(attr, "gene_name"), gene_type = ex(attr, "gene_type"))]
g[, attr := NULL]
g <- canon(unique(g))
genes_gtf <- unique(g[, .(gene_id, gene_name, gene_type)])
cat("GTF:", nrow(g), "TSSs of", nrow(genes_gtf), "genes\n")
tw <- g[, .(chr, start = tss - WIN, end = tss + WIN, tss, gene_id, gene_name, gene_type)]   # 1-based, closed
setkey(tw, chr, start, end)
ot <- foverlaps(S1[, .(chr, start, end, site)], tw, nomatch = NULL)
ot[, dist := pmax(0L, i.start - tss, tss - i.end)]
tss_links <- ot[, .(tss_dist = min(dist)), by = .(site, gene_id, gene_name, gene_type)]

# ---- ABC links
abc_l <- rbindlist(lapply(groups, function(k) {
  f <- file.path(abc, k, "Predictions", "EnhancerPredictionsFull_threshold0.02_self_promoter.tsv")
  d <- fread(f, select = c("chr", "start", "end", "class", "TargetGene", "TargetGeneEnsembl_ID", "isSelfPromoter",
                           "ABC.Score"))
  d[, group := k]
}))
cat("ABC links:", nrow(abc_l), "\n"); print(abc_l[, .N, by = group])
abc_l[, gene_id := sub("\\..*$", "", TargetGeneEnsembl_ID)]
abc_l[genes_gtf, `:=`(gene_name = i.gene_name, gene_type = i.gene_type), on = "gene_id"]
um <- unique(abc_l[is.na(gene_name), .(gene_id, abc_symbol = TargetGene)])
fwrite(um, file.path(od, "unmatched.tsv"), sep = "\t")
cat("ABC genes without a GTF match (kept under the ABC symbol):", nrow(um), "\n")
abc_l[is.na(gene_name), `:=`(gene_name = TargetGene, gene_type = "ABC_only")]
abc_l <- cl(canon(abc_l))
setkey(abc_l, chr, start, end)
oa <- foverlaps(S1[, .(chr, start, end, site)], abc_l, nomatch = NULL)
abc_links <- oa[, .(abc_best = max(ABC.Score), abc_groups = paste(intersect(groups, group), collapse = ","),
                    abc_class = paste(sort(unique(class)), collapse = ","), abc_own_promoter = any(isSelfPromoter)),
                by = .(site, gene_id, gene_name, gene_type)]

# ---- combine
links <- merge(tss_links, abc_links, by = c("site", "gene_id", "gene_name", "gene_type"), all = TRUE)
links[, source := fifelse(!is.na(tss_dist) & !is.na(abc_best), "TSS+ABC", fifelse(!is.na(tss_dist), "TSS", "ABC"))]
links <- merge(sites[, .(site, chr, start, end, days, days_3of3, score)], links, by = "site")
setcolorder(links, c("site", "chr", "start", "end", "days", "days_3of3", "score", "gene_name", "gene_id", "gene_type",
                     "source", "tss_dist", "abc_best", "abc_groups", "abc_class", "abc_own_promoter"))
setorder(links, chr, start, gene_name)
fwrite(links, file.path(od, "links.tsv"), sep = "\t", na = "NA")

sg <- links[, .(n_genes = uniqueN(gene_name), genes = paste(unique(gene_name), collapse = ",")), by = site]
sites[sg, `:=`(n_genes = i.n_genes, genes = i.genes), on = "site"][is.na(n_genes), `:=`(n_genes = 0L, genes = "")]
fwrite(sites[, .(site, chr, start, end, days, days_3of3, score, n_genes, genes)], file.path(od, "sites.tsv"),
       sep = "\t")

dd <- function(x) paste(intersect(DAYS, unlist(strsplit(x, ","))), collapse = ",")
gt <- links[, .(gene_id = gene_id[1], gene_type = gene_type[1], n_sites = uniqueN(site),
                n_sites_tss = uniqueN(site[source != "ABC"]), n_sites_abc_only = uniqueN(site[source == "ABC"]),
                tss_bound = any(source != "ABC"), abc_linked = any(!is.na(abc_best)),
                abc_best = suppressWarnings(max(abc_best, na.rm = TRUE)),
                days = dd(paste(days, collapse = ",")), days_tss = dd(paste(days[source != "ABC"], collapse = ",")),
                best_site_score = max(score)), by = gene_name]
gt[!is.finite(abc_best), abc_best := NA]
gg <- links[!is.na(abc_groups), .(g2 = paste(intersect(groups, unlist(strsplit(abc_groups, ","))), collapse = ",")),
            by = gene_name]
gt[, abc_groups := ""][gg, abc_groups := i.g2, on = "gene_name"]
gt[, link := fifelse(tss_bound & n_sites_abc_only > 0, "TSS + enhancer", fifelse(tss_bound, "TSS only",
              "enhancer only (ABC)"))]
setorder(gt, gene_name)
fwrite(gt, file.path(od, "genes.tsv"), sep = "\t", na = "NA")

# ---- summary
abc_universe <- unique(rbindlist(lapply(groups, function(k)
  fread(file.path(abc, k, "processed_genes_file.bed"), header = FALSE, select = 7, col.names = "gene_id"))))
abc_universe[, gene_id := sub("\\..*$", "", gene_id)]
sm <- rbindlist(list(
  data.table(item = "ISL1 sites", n = nrow(sites)),
  data.table(item = paste("sites per day set:", sites[, .N, by = days]$days), n = sites[, .N, by = days]$N),
  data.table(item = "sites with >= 1 gene", n = sites[n_genes > 0, .N]),
  data.table(item = "sites with a TSS link", n = links[source != "ABC", uniqueN(site)]),
  data.table(item = "sites with an ABC link", n = links[!is.na(abc_best), uniqueN(site)]),
  data.table(item = "sites with an ABC link only (no TSS link to any gene)",
             n = length(setdiff(links[!is.na(abc_best), unique(site)], links[source != "ABC", unique(site)]))),
  data.table(item = "sites without a gene", n = sites[n_genes == 0, .N]),
  data.table(item = paste("site-gene links:", links[, .N, by = source]$source), n = links[, .N, by = source]$N),
  data.table(item = "target genes", n = nrow(gt)),
  data.table(item = paste("target genes:", gt[, .N, by = link]$link), n = gt[, .N, by = link]$N),
  data.table(item = "target genes, protein coding", n = gt[gene_type == "protein_coding", .N]),
  data.table(item = "GTF genes", n = nrow(genes_gtf)),
  data.table(item = "GTF protein-coding genes", n = genes_gtf[gene_type == "protein_coding", .N]),
  data.table(item = "ABC genes (union of groups)", n = nrow(abc_universe)),
  data.table(item = paste("ABC-linked genes supported by", groups),
             n = sapply(groups, function(k) gt[grepl(k, abc_groups, fixed = TRUE), .N]))))
fwrite(sm, file.path(od, "summary.tsv"), sep = "\t")
print(sm)
cat("wrote", od, "\n")
