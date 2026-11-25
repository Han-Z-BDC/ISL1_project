# Interactive D15 gene browser (user, 2026-11-24: "an html report where I can select genes (all genes from the C1-C10
# direct ISL1 targets) to plot in these two methods (trajectory based and feature plots)"). Runs locally (R).
# Genes: the 7.3h C1-C10 trajectory genes with ISL1 binding at any stage (integration step 12 genes.tsv, bound_any:
#   a CUT&RUN site at the TSS or ABC-linked, D15 / D18 / D22).
# Per gene, the page draws the two figures of fig_d15_gene_curves.R and fig_feature_dot.R (D15_endocrine) in the
#   browser from data computed here:
#   trajectory: Slingshot lineage cells with curve weight >= 0.5, 10 equal-width pseudotime bins per lineage (beta /
#     alpha / EC), curve-weighted mean log1p per genotype, NA where < 15 cells, 3-bin moving average; KO stop =
#     KO 95th percentile of pseudotime on beta / alpha. Identical numbers to fig_d15_gene_curves.R.
#   feature: per-cell log1p on the D15 endocrine UMAP (5,994 cells), capped at the 99th percentile of the expressing
#     cells (min 0.1), stored as 8-bit levels of the cap (any expressing cell >= level 1).
# Inputs (git-ignored): reports/trajectory/D15/slingshot_pc30_{cells,lineages}.csv, reports/labels/D15/
#   endocrine_subclusters.csv, reports/figure_data/D15_genes_targets.csv.gz (fig_export_genes.py D15 with the gene
#   list; copied under this name), reports/isl1_binding_C1_C4/step12_genes.tsv (integration step 12 genes.tsv),
#   analysis/scripts/gene_browser_template.html
# Output (git-ignored): reports/isl1_target_gene_browser.html (one self-contained file); reports/gene_browser_artifact/
#   (isl1_gene_browser.html + gene_browser_data.json, the version published as a claude.ai Artifact)
# Usage (from scRNAseq/): Rscript analysis/scripts/fig_d15_gene_browser.R
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
TR <- "reports/trajectory/D15"
EXPR <- Sys.getenv("GB_EXPR", "reports/figure_data/D15_genes_targets.csv.gz")
OUT <- "reports/isl1_target_gene_browser.html"
NBIN <- 10; MIN_CELLS <- 15; W_MIN <- 0.5; SEED <- 1234
LIN <- c(pre_beta2 = "beta", pre_alpha = "alpha", pre_EC = "EC")
GENO <- c("ctr", "KO")
pal <- fread("study_config/annotation/label_colours.tsv"); pal <- setNames(pal$colour, pal$label)

# ---- genes
gi <- fread("reports/isl1_binding_C1_C4/step12_genes.tsv")[grepl("^C[0-9]+$", group) & bound_any == TRUE]
e <- fread(cmd = sprintf("gzcat %s", EXPR))
miss <- setdiff(gi$gene, names(e)); if (length(miss)) cat("not in the expression export, dropped:", miss, "\n")
gi <- gi[gene %in% names(e)]
gi[, cl := as.integer(sub("C", "", group))]; setorder(gi, cl, gene)
G <- gi$gene
cat(length(G), "genes\n")

# ---- cells: slingshot table = the 5,994 D15 endocrine cells
sl <- fread(file.path(TR, "slingshot_pc30_cells.csv")); setnames(sl, 1, "barcode")
ln <- fread(file.path(TR, "slingshot_pc30_lineages.csv"))
lend <- setNames(sapply(strsplit(ln$clusters, " -> "), function(p) p[length(p)]), ln$lineage)
s <- fread("reports/labels/D15/endocrine_subclusters.csv"); setnames(s, 1, "barcode")
stopifnot(setequal(s$barcode, sl$barcode))
E <- as.matrix(e[match(sl$barcode, barcode), ..G]); stopifnot(!anyNA(E))

# ---- trajectory curves (as fig_d15_gene_curves.R, all genes at once)
curves <- list(); mids <- list(); stops <- list()
for (end in names(LIN)) {
  L <- names(lend)[lend == end]
  pt <- sl[[paste0("pt_", L)]]; w <- sl[[paste0("w_", L)]]; w[is.na(w)] <- 0
  on <- w >= W_MIN & !is.na(pt)
  br <- seq(min(pt[on]), max(pt[on]), length.out = NBIN + 1)
  b <- pmin(pmax(findInterval(pt, br, rightmost.closed = TRUE), 1), NBIN)
  mids[[LIN[[end]]]] <- round((br[-1] + br[-(NBIN + 1)]) / 2, 2)
  if (end != "pre_EC") stops[[LIN[[end]]]] <- round(unname(quantile(pt[on & sl$genotype == "KO"], 0.95)), 2)
  for (gt in GENO) {
    V <- matrix(NA_real_, NBIN, length(G))
    for (i in seq_len(NBIN)) {
      m <- on & sl$genotype == gt & b == i
      if (sum(m) >= MIN_CELLS) V[i, ] <- colSums(E[m, , drop = FALSE] * w[m]) / sum(w[m])
    }
    S <- V; for (i in seq_len(NBIN)) { k <- max(1, i - 1):min(NBIN, i + 1)
      S[i, ] <- ifelse(is.na(V[i, ]), NA, colMeans(V[k, , drop = FALSE], na.rm = TRUE)) }
    curves[[sprintf("%s_%s", LIN[[end]], gt)]] <- S
  }
}

# ---- feature plot cells (as fig_feature_dot.R D15_endocrine: same shuffle; UMAP of the endocrine subset)
set.seed(SEED); ord <- sample(nrow(sl))   # draw order only matters for ties; the page sorts by value
ux <- s$UMAP1[match(sl$barcode, s$barcode)]; uy <- s$UMAP2[match(sl$barcode, s$barcode)]
cells <- list(x = round(ux[ord], 3), y = round(uy[ord], 3), ko = as.integer(sl$genotype[ord] == "KO"))
NBYTE <- ceiling(nrow(sl) / 8)
enc <- function(v) {   # bitmap of expressing cells (LSB first) + uint8 levels of those cells, base64
  cap <- max(quantile(v[v > 0], 0.99), 0.1)
  nz <- v > 0; q <- pmax(1L, as.integer(round(pmin(v[nz], cap) / cap * 255)))
  bits <- packBits(c(nz, rep(FALSE, NBYTE * 8 - length(v))), type = "raw")
  list(cap = signif(cap, 3), n = sum(nz), b = base64_enc(c(bits, as.raw(q))))
}

genes <- lapply(seq_along(G), function(j) {
  f <- enc(E[ord, j])
  cv <- lapply(curves, function(S) ifelse(is.na(S[, j]), NA, round(S[, j], 4)))
  list(g = G[j], c = gi$cl[j],
       b = paste0(as.integer(gi$bound_D15[j]), as.integer(gi$bound_D18[j]), as.integer(gi$bound_D22[j])),
       t = paste0(as.integer(gi$tss_D15[j]), as.integer(gi$tss_D18[j]), as.integer(gi$tss_D22[j])),
       s = gi$n_sites[j], cap = f$cap, n = f$n, e = f$b, cv = cv)
})
data <- list(cells = cells, mids = mids, stops = stops, genes = genes,
             lcol = list(beta = unname(pal["pre_beta2"]), alpha = unname(pal["pre_alpha"]), EC = unname(pal["pre_EC"])),
             nctr = sum(sl$genotype == "ctr"), nko = sum(sl$genotype == "KO"), date = format(Sys.Date()))
js <- toJSON(data, auto_unbox = TRUE, na = "null", digits = NA)
tpl <- paste(readLines("analysis/scripts/gene_browser_template.html", warn = FALSE), collapse = "\n")
parts <- strsplit(tpl, "/*__DATA__*/", fixed = TRUE)[[1]]; stopifnot(length(parts) == 2)
html <- paste0(parts[1], "const DATA = ", js, ";", parts[2])
writeLines(html, OUT)
cat(sprintf("wrote %s (%.1f MB)\n", OUT, file.size(OUT) / 1e6))
# Artifact version (user, 2026-11-24: "publish the gene browser as an artifact"): page without data + the data as
# gene_browser_data.json beside it (fetched by the page)
AD <- "reports/gene_browser_artifact"; dir.create(AD, showWarnings = FALSE)
writeLines(paste0(parts[1], "const DATA = null;", parts[2]), file.path(AD, "isl1_gene_browser.html"))
writeLines(js, file.path(AD, "gene_browser_data.json"))
cat(sprintf("wrote %s/ (page + %.1f MB data)\n", AD, file.size(file.path(AD, "gene_browser_data.json")) / 1e6))
