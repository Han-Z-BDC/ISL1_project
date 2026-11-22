# Step 7.1n (GO): GO biological-process over-representation of the 7.1n gene sets (after 7.1n_ec_deseq2.R).
# Offline: GO annotation from org.Hs.eg.db (GOALL, ontology BP, i.e. genes annotated to a term or its descendants),
# term names from GO.db; fgsea::fora (hypergeometric), universe = genes tested in 7.1n; terms with 15-500 universe
# genes; BH padj over the terms tested per set.
# Sets: core EC2-high, core EC1-high (genes_celltype.csv); KO up / down, shared and per cell type (genes_KO.csv).
# Writes ${OUT_DIR}/day_analysis/ec_deseq2/go_bp.csv (set, term, GO id, overlap, set size, term size, pval, padj,
#   genes) and prints the top terms per set.
# Usage: Rscript 7.1n_ec_go.R   (isl1_downstream env; see 7.1n_ec_deseq2.sb)
suppressPackageStartupMessages({ library(data.table); library(fgsea); library(AnnotationDbi); library(org.Hs.eg.db)
                                 library(GO.db) })
od <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "ec_deseq2")
uni <- fread(file.path(od, "results_B.csv"), select = "gene")$gene
gA <- fread(file.path(od, "genes_celltype.csv")); gB <- fread(file.path(od, "genes_KO.csv"))
sets <- list(
  "EC2-high (core)" = gA[class == "EC2-high", gene], "EC1-high (core)" = gA[class == "EC1-high", gene],
  "KO up, shared EC1 + EC2" = gB[KO_class == "shared" & KO_direction == "up in KO", gene],
  "KO down, shared EC1 + EC2" = gB[KO_class == "shared" & KO_direction == "down in KO", gene],
  "KO up, EC1 only" = gB[KO_class == "EC1 only" & KO_direction == "up in KO", gene],
  "KO down, EC1 only" = gB[KO_class == "EC1 only" & KO_direction == "down in KO", gene],
  "KO up, EC2 only" = gB[KO_class == "EC2 only" & KO_direction == "up in KO", gene],
  "KO down, EC2 only" = gB[KO_class == "EC2 only" & KO_direction == "down in KO", gene])
print(sapply(sets, length))

an <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = intersect(uni, keys(org.Hs.eg.db, "SYMBOL")),
                                             keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")))
an <- as.data.table(an)[ONTOLOGYALL == "BP" & !is.na(GOALL), .(gene = SYMBOL, go = GOALL)]
an <- unique(an)
pw <- split(an$gene, an$go)
pw <- pw[lengths(pw) >= 15 & lengths(pw) <= 500]
cat("GO BP terms tested:", length(pw), "; universe genes with a BP term:", uniqueN(an$gene), "of", length(uni), "\n")
res <- rbindlist(lapply(names(sets), function(sn) {
  r <- as.data.table(fora(pathways = pw, genes = intersect(sets[[sn]], uni), universe = uni, minSize = 15,
                          maxSize = 500))
  if (!nrow(r)) return(NULL)
  r[, `:=`(set = sn, set_size = length(intersect(sets[[sn]], uni)),
           genes = sapply(overlapGenes, paste, collapse = ","))][, overlapGenes := NULL]
}))
res[, term := Term(GOTERM[pathway])]
setnames(res, c("pathway", "size"), c("go_id", "term_size"))
setcolorder(res, c("set", "term", "go_id", "overlap", "set_size", "term_size", "pval", "padj", "genes"))
res <- res[order(set, pval)]
fwrite(res, file.path(od, "go_bp.csv"))
for (sn in names(sets)) {
  x <- res[set == sn & padj < 0.05]
  cat(sprintf("\n== %s (%d genes): %d terms padj < 0.05\n", sn, length(sets[[sn]]), nrow(x)))
  if (nrow(x)) print(head(x[, .(term = substr(term, 1, 60), overlap, term_size, padj = signif(padj, 2))], 12))
}
cat("wrote", file.path(od, "go_bp.csv"), "\n")
