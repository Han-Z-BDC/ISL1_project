# Step 7.1o (GO): GO biological-process over-representation of the D22 EC1 beta-axis quadrant categories (after
# 7.1o_ec_beta_model.R; user, 2026-11-23: "GO for each gene category, colour based").
# Sets (quadrants.csv, EC1): beta-high up in KO EC, beta-high down in KO EC, EC-high down in KO EC, EC-high up in KO EC.
# As 7.1n_ec_go.R: GO annotation from org.Hs.eg.db (GOALL, BP), term names from GO.db, fgsea::fora (hypergeometric),
# terms with 15-500 universe genes, BH padj per set. Universe = all genes tested in the 7.1o model (de_wide.csv).
# Writes ${OUT_DIR}/day_analysis/ec_beta_D22/go_bp_quadrants_EC1.csv (set, term, GO id, overlap, set size, term size,
#   pval, padj, genes) and prints the top terms per set.
# Run: sbatch analysis/scripts/7.1o_ec_beta_go.sb   (isl1_downstream env)
suppressPackageStartupMessages({ library(data.table); library(fgsea); library(AnnotationDbi); library(org.Hs.eg.db)
                                 library(GO.db) })
od <- file.path(Sys.getenv("OUT_DIR"), "day_analysis", "ec_beta_D22")
uni <- fread(file.path(od, "de_wide.csv"), select = "gene")$gene
q <- fread(file.path(od, "quadrants.csv"))[ec == "EC1"]
CATS <- c("beta-high, up in KO EC", "beta-high, down in KO EC (check ambient)", "EC-high, down in KO EC",
          "EC-high, up in KO EC")
sets <- setNames(lapply(CATS, function(cn) q[quadrant == cn, gene]), CATS)
print(sapply(sets, length))

an <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = intersect(uni, keys(org.Hs.eg.db, "SYMBOL")),
                                             keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")))
an <- unique(as.data.table(an)[ONTOLOGYALL == "BP" & !is.na(GOALL), .(gene = SYMBOL, go = GOALL)])
pw <- split(an$gene, an$go)
pw <- pw[lengths(pw) >= 15 & lengths(pw) <= 500]
cat("GO BP terms tested:", length(pw), "; universe genes with a BP term:", uniqueN(an$gene), "of", length(uni), "\n")
res <- rbindlist(lapply(names(sets), function(sn) {
  r <- as.data.table(fora(pathways = pw, genes = intersect(sets[[sn]], uni), universe = uni, minSize = 15, maxSize = 500))
  if (!nrow(r)) return(NULL)
  r[, `:=`(set = sn, set_size = length(intersect(sets[[sn]], uni)),
           genes = sapply(overlapGenes, paste, collapse = ","))][, overlapGenes := NULL]
}))
res[, term := Term(GOTERM[pathway])]
setnames(res, c("pathway", "size"), c("go_id", "term_size"))
setcolorder(res, c("set", "term", "go_id", "overlap", "set_size", "term_size", "pval", "padj", "genes"))
res <- res[order(set, pval)]
fwrite(res, file.path(od, "go_bp_quadrants_EC1.csv"))
for (sn in names(sets)) {
  x <- res[set == sn & padj < 0.05]
  cat(sprintf("\n== %s (%d genes): %d terms padj < 0.05\n", sn, length(sets[[sn]]), nrow(x)))
  if (nrow(x)) print(head(x[, .(term = substr(term, 1, 60), overlap, term_size, padj = signif(padj, 2),
                                genes = substr(genes, 1, 70))], 12))
}
cat("wrote", file.path(od, "go_bp_quadrants_EC1.csv"), "\n")
