#!/usr/bin/env Rscript
# Step 08: group all ISL1 sites by their binding pattern across D15 / D18 / D22.
#   Counts     step 27's DiffBind fragment counts and background size factors (step 00), 9 libraries.
#   Varying    DESeq2 LRT, ~ batch + day against ~ batch, size factors fixed to the DiffBind ones. Sites with
#              LRT p < 0.05 are clustered (a filter against noise, not a significance claim); the others form the
#              "stable" group (cluster 0).
#   Profile    DESeq2 VST (blind = FALSE), batch 3 effect removed (limma removeBatchEffect, day kept in the
#              design), mean per day, z-scored per site across the 3 days.
#   Clusters   k-means on the z-scored profiles of the varying sites, k = 2-10 (seed 1, 25 starts, 100 iterations);
#              clusters numbered by size within each k. The report uses one k (chosen from the elbow and from
#              whether extra clusters only split one shape).
#   Held-out   the clustering is repeated on the 6 earlier libraries only (LRT ~ day on 6 libraries, p < 0.05; VST
#              day means of those libraries). For each varying site, the z-scored profile of its 3 batch 3 libraries
#              (D15 rep1, D18 rep1, D22 rep2; VST of the 9-library fit) is compared with the 6-library cluster
#              centres: share of sites whose batch 3 profile is closest (Pearson r) to their own cluster's centre,
#              against the share expected by chance (sum of squared cluster fractions).
#   Links      per cluster: sites gained / lost in step 27 (step 03 site_stats), sites linked (GREAT, step 01) to
#              each bulk KO-effect group (step 02) and to the alpha / beta markers; Fisher test of each cluster
#              against the stable sites for each group.
#   Genes      gene level, for the report's k (K_GENE = 4): linked genes tested in the bulk RNA-seq (step 02); a gene is
#              in a site group if at least one of its GREAT-linked sites is (a gene can be in several groups); the
#              reference is the genes whose sites are all stable. Per site group and KO-effect group: Fisher test,
#              and logistic regression in_site_group ~ in_KO_group + log2(1 + n_sites) + log10(domain_kb) (genes
#              with more sites more often have a varying site, as in step 02), BH over all tests. Done twice: genes
#              linked by GREAT (main), and every bulk-tested gene with its TSS within 50 kb of the site summit (check;
#              covariate log2(1 + n_sites) only).
# usage: Rscript 08_site_clusters.R <step 27 isl1_diff dir> <integration results dir> <markers.tsv> <out_dir>
# Output in <out_dir>: site_clusters.tsv, elbow.tsv, centers.tsv, heldout.tsv (per site), heldout_clusters.tsv
#   (mean batch 3 profile per 6-library cluster), links.tsv, vst_libraries.tsv, gene_site_groups.tsv (per link and
#   gene: KO group, marker, sites per site group), gene_links.tsv (per link, site group and KO group)

suppressPackageStartupMessages({ library(DESeq2); library(limma); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
cnr <- args[1]; res <- args[2]; mfile <- args[3]; od <- args[4]
dir.create(od, showWarnings = FALSE, recursive = TRUE)
days <- c("D15", "D18", "D22"); P_VAR <- 0.05; KS <- 2:10; K_GENE <- 4

cs <- fread(file.path(res, "00_diffbind_counts", "counts_raw.tsv"))
m <- as.matrix(cs[, -1]); rownames(m) <- cs$name
smp <- fread(file.path(cnr, "samples.tsv"))
stopifnot(identical(colnames(m), smp$sample))
nf <- fread(file.path(res, "00_diffbind_counts", "norm_factors.tsv"))
sf <- setNames(nf$norm_factor[match(smp$sample, nf$sample)], smp$sample)
stopifnot(!anyNA(sf))
cd <- data.frame(day = factor(smp$day, levels = days), batch = factor(smp$batch, levels = c("earlier", "batch3")),
                 row.names = smp$sample)

# ---- LRT and VST, 9 libraries
dds <- DESeqDataSetFromMatrix(m, cd, design = ~ batch + day)
sizeFactors(dds) <- sf
dds <- DESeq(dds, test = "LRT", reduced = ~ batch, quiet = TRUE)
lrt <- results(dds)
vs <- assay(vst(dds, blind = FALSE))
vc <- removeBatchEffect(vs, batch = cd$batch, design = model.matrix(~ day, cd))
dm <- sapply(days, function(d) rowMeans(vc[, cd$day == d, drop = FALSE]))
zsc <- function(x) { s <- apply(x, 1, sd); (x - rowMeans(x)) / ifelse(s > 0, s, 1) }
z <- zsc(dm)
varying <- !is.na(lrt$pvalue) & lrt$pvalue < P_VAR
cat(sprintf("sites %d; LRT p < %.2f: %d; FDR < 0.1: %d\n", nrow(m), P_VAR, sum(varying),
            sum(lrt$padj < 0.1, na.rm = TRUE)))

km_by_size <- function(x, k) {
  set.seed(1)
  fit <- kmeans(x, centers = k, nstart = 25, iter.max = 100)
  sz <- tabulate(fit$cluster, k)
  ord <- order(-sz, seq_len(k))                       # by size, ties by first-found order
  list(cluster = match(fit$cluster, ord), centers = fit$centers[ord, , drop = FALSE],
       tot_within = fit$tot.withinss, between_share = fit$betweenss / fit$totss)
}

st <- data.table(site = rownames(m), lrt_stat = round(lrt$stat, 3), lrt_p = signif(lrt$pvalue, 4),
                 lrt_fdr = signif(lrt$padj, 4), varying = varying)
for (d in days) st[, (paste0("vst_", d)) := round(dm[, d], 4)]
for (d in days) st[, (paste0("z_", d)) := round(z[, d], 4)]
elbow <- list(); centers <- list(); fits <- list()
for (k in KS) {
  f <- km_by_size(z[varying, ], k); fits[[as.character(k)]] <- f
  cl <- integer(nrow(m)); cl[varying] <- f$cluster
  st[, (paste0("k", k)) := cl]
  elbow[[length(elbow) + 1]] <- data.table(k = k, tot_within_ss = round(f$tot_within, 2),
                                           between_share = round(f$between_share, 4))
  centers[[length(centers) + 1]] <- data.table(k = k, cluster = seq_len(k), n_sites = tabulate(f$cluster, k),
    z_D15 = round(f$centers[, 1], 3), z_D18 = round(f$centers[, 2], 3), z_D22 = round(f$centers[, 3], 3),
    vst_D15 = round(sapply(seq_len(k), function(c) median(dm[varying, 1][f$cluster == c])), 3),
    vst_D18 = round(sapply(seq_len(k), function(c) median(dm[varying, 2][f$cluster == c])), 3),
    vst_D22 = round(sapply(seq_len(k), function(c) median(dm[varying, 3][f$cluster == c])), 3))
}
centers <- rbindlist(centers)
centers <- rbind(centers, data.table(k = NA_integer_, cluster = 0L, n_sites = sum(!varying),
  z_D15 = NA_real_, z_D18 = NA_real_, z_D22 = NA_real_,
  vst_D15 = round(median(dm[!varying, 1]), 3), vst_D18 = round(median(dm[!varying, 2]), 3),
  vst_D22 = round(median(dm[!varying, 3]), 3)))
fwrite(st, file.path(od, "site_clusters.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(elbow), file.path(od, "elbow.tsv"), sep = "\t", quote = FALSE)
fwrite(centers, file.path(od, "centers.tsv"), sep = "\t", quote = FALSE)

# ---- held-out check: cluster with the 6 earlier libraries, test on the 3 batch 3 libraries
e6 <- cd$batch == "earlier"
d6 <- DESeqDataSetFromMatrix(m[, e6], droplevels(cd[e6, ]), design = ~ day)
sizeFactors(d6) <- sf[e6]
d6 <- DESeq(d6, test = "LRT", reduced = ~ 1, quiet = TRUE)
p6 <- results(d6)$pvalue
v6 <- !is.na(p6) & p6 < P_VAR
vs6 <- assay(vst(d6, blind = FALSE))
z6 <- zsc(sapply(days, function(d) rowMeans(vs6[, cd$day[e6] == d, drop = FALSE])))
b3 <- smp[batch == "batch3"]
zb3 <- zsc(vs[, b3$sample[match(days, b3$day)]])     # VST of the 9-library fit, batch 3 libraries, D15/D18/D22
hoc <- list()
ho <- rbindlist(lapply(KS, function(k) {
  f <- km_by_size(z6[v6, ], k)
  # per cluster: mean batch 3 profile over the cluster's sites against the 6-library centre
  zb <- zb3[v6, , drop = FALSE]
  hoc[[length(hoc) + 1]] <<- rbindlist(lapply(seq_len(k), function(c) {
    mb <- colMeans(zb[f$cluster == c, , drop = FALSE], na.rm = TRUE)
    data.table(k = k, cluster = c, n_sites = sum(f$cluster == c),
               centre_D15 = round(f$centers[c, 1], 3), centre_D18 = round(f$centers[c, 2], 3),
               centre_D22 = round(f$centers[c, 3], 3), batch3_D15 = round(mb[1], 3), batch3_D18 = round(mb[2], 3),
               batch3_D22 = round(mb[3], 3), pearson_r = round(cor(f$centers[c, ], mb), 3))
  }))
  r <- cor(t(zb3[v6, ]), t(f$centers))                 # sites x centres
  near <- max.col(r, ties.method = "first")
  frac <- tabulate(f$cluster, k) / sum(v6)
  # agreement with the 9-library clusters for sites varying in both
  both <- v6 & varying
  data.table(k = k, sites_6lib = sum(v6), pct_batch3_nearest_own = round(100 * mean(near == f$cluster, na.rm = TRUE), 1),
             pct_expected_by_chance = round(100 * sum(frac^2), 1),
             sites_varying_both = sum(both),
             pct_same_as_9lib_best_match = {
               x <- table(f$cluster[both[v6]], st[[paste0("k", k)]][both])
               round(100 * sum(apply(x, 1, max)) / sum(x), 1) })
}))
fwrite(ho, file.path(od, "heldout.tsv"), sep = "\t", quote = FALSE)
hoc <- rbindlist(hoc)
fwrite(hoc, file.path(od, "heldout_clusters.tsv"), sep = "\t", quote = FALSE)
# all sites, per day: median batch 3 z-profile (the shared day trend of the held-out libraries)
cat("median batch 3 z per day, all sites:", round(apply(zb3, 2, median, na.rm = TRUE), 3), "\n")
cat("median 6-library z per day, all sites:", round(apply(z6, 2, median, na.rm = TRUE), 3), "\n")
# per library VST (9-library fit), for the report heatmap
fwrite(data.table(site = rownames(vs), round(vs, 4)), file.path(od, "vst_libraries.tsv"), sep = "\t", quote = FALSE)

# ---- links to gained / lost sites, KO-effect groups and markers
ss <- fread(file.path(res, "03_binding_change", "site_stats.tsv"))[, .(site, gained, lost)]
sgl <- fread(file.path(res, "01_great", "site_gene.tsv"), sep = "\t", quote = "")
gt <- fread(file.path(res, "02_direct_targets", "gene_table.tsv"), sep = "\t", quote = "", na.strings = "NA")
mk <- fread(mfile, skip = "cell_type")
sgl[gt, group := i.group, on = "gene_id"]
feat <- list(gained_step27 = ss[gained == TRUE, site], lost_step27 = ss[lost == TRUE, site])
for (g in c(paste0("cluster", 1:6), "stays", "other_KO", "no_KO_effect"))
  feat[[paste0("KO_", g)]] <- unique(sgl[group == g, site])
feat$alpha_markers <- unique(sgl[gene_name %in% mk[cell_type == "alpha", gene], site])
feat$beta_markers <- unique(sgl[gene_name %in% mk[cell_type == "beta", gene], site])
links <- rbindlist(lapply(KS, function(k) {
  cl <- st[[paste0("k", k)]]
  rbindlist(lapply(names(feat), function(fn) {
    has <- st$site %in% feat[[fn]]
    rbindlist(lapply(0:k, function(c) {
      a <- sum(has & cl == c); n <- sum(cl == c)
      ft <- if (c == 0) NULL else fisher.test(matrix(c(a, n - a, sum(has & cl == 0), sum(cl == 0) - sum(has & cl == 0)), 2))
      data.table(k = k, feature = fn, cluster = c, n_sites = n, n_with = a, pct_with = round(100 * a / n, 2),
                 odds_ratio_vs_stable = if (is.null(ft)) NA_real_ else signif(unname(ft$estimate), 3),
                 p_vs_stable = if (is.null(ft)) NA_real_ else signif(ft$p.value, 3))
    }))
  }))
}))
links[!is.na(p_vs_stable), fdr_vs_stable := signif(p.adjust(p_vs_stable, "BH"), 3), by = k]
fwrite(links, file.path(od, "links.tsv"), sep = "\t", quote = FALSE)

# ---- gene level: site groups of each gene's sites against its KO-effect group
#   link "GREAT": step 01 domains; link "within_50kb" (check): every bulk-tested gene whose TSS is within 50 kb of the
#   site summit (a site can link to several genes; catches sites GREAT gives to a neighbouring gene, e.g. near PAX6)
kc <- paste0("k", K_GENE)
reg <- fread(file.path(res, "01_great", "gene_regions.bed"), header = FALSE,
             col.names = c("chr", "start", "end", "gene_id", "gene_name", "strand", "tss", "basal_start",
                           "basal_end", "gene_type", "tested"))
sb <- fread(file.path(cnr, "sites_used.bed"), header = FALSE, select = c(1:4),
            col.names = c("chr", "start", "end", "site"))
sb[, `:=`(s0 = start + (end - start) %/% 2, s1 = start + (end - start) %/% 2)]
tw <- reg[gene_id %in% gt$gene_id, .(chr, s0 = pmax(tss - 50000L, 0L), s1 = tss + 50000L, gene_id)]
setkey(tw, chr, s0, s1)
w50 <- unique(foverlaps(sb[, .(chr, s0, s1, site)], tw, nomatch = NULL)[, .(site, gene_id)])
links_by <- list(GREAT = sgl[tested_in_bulk == 1, .(site, gene_id)], within_50kb = w50)
kgroups <- c(paste0("cluster", 1:6), "stays", "other_KO", "no_KO_effect")
gg_all <- list(); gl_all <- list()
for (ln in names(links_by)) {
  gs <- merge(links_by[[ln]], st[, c("site", kc), with = FALSE], by = "site")
  setnames(gs, kc, "sgroup")
  gs[, sgroup := factor(sgroup, levels = 0:K_GENE)]
  gw <- dcast(gs, gene_id ~ sgroup, fun.aggregate = length, value.var = "site", drop = FALSE)
  setnames(gw, as.character(0:K_GENE), paste0("sites_group", 0:K_GENE))
  gg <- merge(gt[, .(gene_id, gene_name, group, domain_kb)], gw, by = "gene_id")
  gg[, n_sites := rowSums(.SD), .SDcols = patterns("^sites_group")]
  gg <- gg[n_sites > 0]
  stopifnot(!anyDuplicated(gg$gene_id))
  if (ln == "GREAT") stopifnot(all(gg$n_sites == gt[match(gg$gene_id, gene_id), n_sites]))
  gg[mk, marker := i.cell_type, on = c(gene_name = "gene")]
  gg[, `:=`(link = ln, stable_only = sites_group0 == n_sites)]
  # covariates: number of linked sites; for GREAT also the domain size (the 50 kb window is the same for all genes)
  cov <- if (ln == "GREAT") "+ log2(1 + n_sites) + log10(domain_kb)" else "+ log2(1 + n_sites)"
  gl_all[[ln]] <- rbindlist(lapply(seq_len(K_GENE), function(c) {
    sub <- gg[get(paste0("sites_group", c)) > 0 | stable_only]
    sub[, y := as.integer(get(paste0("sites_group", c)) > 0)]
    rbindlist(lapply(kgroups, function(kg) {
      sub[, x := as.integer(group == kg)]
      ft <- fisher.test(table(factor(sub$x, 1:0), factor(sub$y, 1:0)))
      fit <- glm(as.formula(paste("y ~ x", cov)), data = sub, family = binomial)
      cf <- summary(fit)$coefficients["x", ]
      data.table(link = ln, site_group = c, ko_group = kg, genes = sub[x == 1 & y == 1, .N],
                 genes_in_site_group = sum(sub$y), pct_in_site_group = round(100 * sub[x == 1 & y == 1, .N] / sum(sub$y), 2),
                 pct_stable_only = round(100 * sub[x == 1 & y == 0, .N] / sum(sub$y == 0), 2),
                 odds_ratio = signif(unname(ft$estimate), 3), p = signif(ft$p.value, 3),
                 adj_odds_ratio = signif(exp(cf[["Estimate"]]), 3), adj_p = signif(cf[["Pr(>|z|)"]], 3))
    }))
  }))
  gg_all[[ln]] <- gg
  cat(sprintf("%s: linked genes %d; with only stable sites %d\n", ln, nrow(gg), sum(gg$stable_only)))
}
gg <- rbindlist(gg_all, use.names = TRUE)
setcolorder(gg, c("link", "gene_id", "gene_name", "group", "marker", "n_sites"))
setorder(gg, link, group, gene_name)
fwrite(gg, file.path(od, "gene_site_groups.tsv"), sep = "\t", quote = FALSE)
gl <- rbindlist(gl_all)
gl[, `:=`(fdr = signif(p.adjust(p, "BH"), 3), adj_fdr = signif(p.adjust(adj_p, "BH"), 3)), by = link]
fwrite(gl, file.path(od, "gene_links.tsv"), sep = "\t", quote = FALSE)
print(gl[adj_fdr < 0.1 | fdr < 0.1])
print(gg[gene_name %in% c("INS", "PAX6", "ARX", "PAX4", "MAFA", "MAFB", "PDX1", "NKX6-1")])

print(rbindlist(elbow)); print(ho); print(centers[k %in% c(4, 5, 6, NA)])
