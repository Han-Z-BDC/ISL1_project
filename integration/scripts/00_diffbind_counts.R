#!/usr/bin/env Rscript
# Step 00: ISL1 site counts and normalization from step 27's DiffBind object, so every integration step uses the
# same counts as the differential binding calls.
#   Input: <isl1_diff>/diffbind.counted.rds (step 27: 9 control ISL1 libraries, paired-end fragments counted on the
#   14,431 sites, duplicates kept, no control subtraction) and <isl1_diff>/sites_used.bed (site names).
#   Normalization as step 27's main analysis: DESeq2 native, background (size factors from 15 kb genome bins).
#   Normalized signal per site and library = reads / size factor (the size factors are DiffBind's norm.factors).
# usage: Rscript 00_diffbind_counts.R <step 27 isl1_diff dir> <out_dir>
# Output in <out_dir>: counts_raw.tsv (site x library reads), norm_factors.tsv (sample, lib_size, norm_factor),
#   counts_norm.tsv (DiffBind's normalized scores, kept to check the reads / size factor values)

suppressPackageStartupMessages({ library(DiffBind); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
d <- args[1]; od <- args[2]
dir.create(od, showWarnings = FALSE, recursive = TRUE)

db <- readRDS(file.path(d, "diffbind.counted.rds"))
db$config$cores <- as.integer(Sys.getenv("DB_CORES", "4"))
sites <- fread(file.path(d, "sites_used.bed"), header = FALSE,
               col.names = c("chr", "start", "end", "name", "score", "strand", "summit_day", "days_called", "n_days"))
name_of <- setNames(sites$name, paste(sites$chr, sites$start + 1, sep = ":"))
smp <- fread(file.path(d, "samples.tsv"))

dbn <- dba.normalize(db, method = DBA_DESEQ2, normalize = DBA_NORM_NATIVE, library = DBA_LIBSIZE_BACKGROUND,
                     background = TRUE)
f <- dba.normalize(dbn, method = DBA_DESEQ2, bRetrieve = TRUE)
nf <- data.table(sample = dbn$samples$SampleID, lib_size = f$lib.sizes, norm_factor = f$norm.factors)
stopifnot(setequal(nf$sample, smp$sample))
nf <- nf[match(smp$sample, sample)]
fwrite(nf, file.path(od, "norm_factors.tsv"), sep = "\t", quote = FALSE)

# raw reads from each library's count table in the object (re-scoring with dba.count would apply DiffBind's
# default filter again and drop the 4 sites with no reads in any library)
raw <- rbindlist(lapply(seq_along(dbn$peaks), function(i) {
  p <- as.data.table(dbn$peaks[[i]])
  data.table(name = name_of[paste(p$Chr, p$Start, sep = ":")], sample = dbn$samples$SampleID[i], reads = p$Reads)
}))
stopifnot(!anyNA(raw$name), raw[, .N, by = sample][, all(N == nrow(sites))])
raw <- dcast(raw, name ~ sample, value.var = "reads")[match(sites$name, name), c("name", smp$sample), with = FALSE]
x <- as.data.table(dba.peakset(dba.count(dbn, peaks = NULL, score = DBA_SCORE_NORMALIZED, filter = 0), bRetrieve = TRUE,
                               DataType = DBA_DATA_FRAME))
cat(sprintf("normalized scores retrieved for %d of %d sites\n", nrow(x), nrow(sites)))
norm <- data.table(name = name_of[paste(x$CHR, x$START, sep = ":")], x[, smp$sample, with = FALSE])
norm <- norm[match(sites$name, name)]   # sites DiffBind does not return stay NA
norm[, name := sites$name]
fwrite(raw, file.path(od, "counts_raw.tsv"), sep = "\t", quote = FALSE)
fwrite(norm, file.path(od, "counts_norm.tsv"), sep = "\t", quote = FALSE)

# check: DiffBind's normalized scores against reads / size factor
r <- as.matrix(raw[, -1]); n <- as.matrix(norm[, -1])
mine <- sweep(r, 2, nf$norm_factor, "/")
ratio <- n / mine
cat(sprintf("sites %d, libraries %d; DiffBind normalized / (reads / size factor): median %.4f, range %.4f-%.4f\n",
            nrow(r), ncol(r), median(ratio, na.rm = TRUE), min(ratio, na.rm = TRUE), max(ratio, na.rm = TRUE)))
cat(sprintf("Pearson r (log, all values): %.4f\n", cor(log1p(as.vector(n)), log1p(as.vector(mine)), use = "complete.obs")))
print(nf)
