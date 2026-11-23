# Step 7.1o (enrichment): over-representation of the D22 EC1 beta-axis quadrant categories in GO and pathway gene sets
# (after 7.1o_ec_beta_model.R; user, 2026-11-23: GO BP, MF and CC plus pathways, KEGG, Reactome and Hallmark, in one
# plot). Replaces the GO-BP-only 7.1o_ec_beta_go.R.
# Gene sets: MSigDB via msigdbr (isl1_downstream env, installed 2026-11-23): GO:BP, GO:MF, GO:CC (C5), KEGG Medicus and
#   KEGG Legacy, Reactome, WikiPathways (C2 canonical pathways), Hallmark (H), human symbols. msigdbr downloads the
#   MSigDB release at first use, so the sets are written once to ${OUT_DIR}/reference/msigdb_<version>.csv.gz (e.g. msigdb_2026.1.Hs.csv.gz) and
#   read from there afterwards (no internet needed then).
# Sets tested (quadrants.csv, EC1): beta-high up in KO EC, beta-high down in KO EC, EC-high down in KO EC, EC-high up
#   in KO EC. fgsea::fora (hypergeometric), universe = all genes tested in the 7.1o model (de_wide.csv), gene sets
#   with 15-500 universe genes, BH padj per category x collection.
# Writes ${OUT_DIR}/day_analysis/ec_beta_D22/enrich_quadrants_EC1.csv (category, collection, term, overlap, set size,
#   term size, pval, padj, genes) and prints the top terms.
# Run: sbatch analysis/scripts/7.1o_ec_beta_enrich.sb   (isl1_downstream env; first run on acompile for the download)
suppressPackageStartupMessages({ library(data.table); library(fgsea) })
out <- Sys.getenv("OUT_DIR"); od <- file.path(out, "day_analysis", "ec_beta_D22")
COLL <- list("GO BP" = c("C5", "GO:BP"), "GO MF" = c("C5", "GO:MF"), "GO CC" = c("C5", "GO:CC"),
             "KEGG Medicus" = c("C2", "CP:KEGG_MEDICUS"), "KEGG Legacy" = c("C2", "CP:KEGG_LEGACY"),
             "Reactome" = c("C2", "CP:REACTOME"), "WikiPathways" = c("C2", "CP:WIKIPATHWAYS"), "Hallmark" = c("H", ""))
ref <- file.path(out, "reference"); dir.create(ref, showWarnings = FALSE, recursive = TRUE)
gsf <- Sys.glob(file.path(ref, "msigdb_*.Hs.csv.gz"))
if (!length(gsf)) {
  suppressPackageStartupMessages(library(msigdbr))
  ver <- as.character(msigdbr_collections()$db_version[1])
  g <- rbindlist(lapply(names(COLL), function(cn) {
    c2 <- COLL[[cn]]
    m <- if (c2[2] == "") msigdbr(species = "Homo sapiens", collection = c2[1]) else
      msigdbr(species = "Homo sapiens", collection = c2[1], subcollection = c2[2])
    unique(data.table(collection = cn, term = m$gs_name, gene = m$gene_symbol))
  }))
  gsf <- file.path(ref, sprintf("msigdb_%s.csv.gz", ver)); fwrite(g, gsf)
  cat("wrote", gsf, "\n")
} else g <- fread(gsf[1])
cat("gene sets:", gsf[1], "\n"); print(g[, .(sets = uniqueN(term)), by = collection])

uni <- fread(file.path(od, "de_wide.csv"), select = "gene")$gene
q <- fread(file.path(od, "quadrants.csv"))[ec == "EC1"]
CATS <- c("beta-high, up in KO EC", "beta-high, down in KO EC (check ambient)", "EC-high, down in KO EC",
          "EC-high, up in KO EC")
sets <- setNames(lapply(CATS, function(cn) q[quadrant == cn, gene]), CATS)
print(sapply(sets, length))
g <- g[gene %in% uni]
res <- rbindlist(lapply(names(COLL), function(cn) {
  pw <- split(g[collection == cn, gene], g[collection == cn, term])
  pw <- pw[lengths(pw) >= 15 & lengths(pw) <= 500]
  rbindlist(lapply(CATS, function(sn) {
    r <- as.data.table(fora(pathways = pw, genes = intersect(sets[[sn]], uni), universe = uni, minSize = 15, maxSize = 500))
    if (!nrow(r)) return(NULL)
    r[, `:=`(category = sn, collection = cn, set_size = length(intersect(sets[[sn]], uni)), terms_tested = length(pw),
             genes = sapply(overlapGenes, paste, collapse = ","))][, overlapGenes := NULL]
  }))
}))
setnames(res, c("pathway", "size"), c("term", "term_size"))
setcolorder(res, c("category", "collection", "term", "overlap", "set_size", "term_size", "pval", "padj", "terms_tested", "genes"))
res <- res[order(category, collection, pval)]
fwrite(res, file.path(od, "enrich_quadrants_EC1.csv"))
print(dcast(res[padj < 0.05, .N, by = .(category, collection)], collection ~ category, value.var = "N", fill = 0))
for (sn in CATS) for (cn in names(COLL)) {
  x <- res[category == sn & collection == cn & padj < 0.05]
  if (nrow(x)) {
    cat(sprintf("\n== %s | %s: %d terms\n", sn, cn, nrow(x)))
    print(head(x[, .(term = substr(term, 1, 60), overlap, padj = signif(padj, 2), genes = substr(genes, 1, 60))], 5))
  }
}
cat("wrote", file.path(od, "enrich_quadrants_EC1.csv"), "\n")
