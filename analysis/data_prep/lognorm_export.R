#!/usr/bin/env Rscript
# Depth-normalize a raw genes x cells count CSV for the SMART-seq2 FUCCI dataset
# (non-UMI -> SCTransform inappropriate). One normalized matrix feeds BOTH the xGATE
# graph and the competing methods.
#   mode "cp10k"   : library-size normalized COUNTS (count-scale, no log) -- the SMART-seq2
#                    analog of SCT corrected counts; correct input for xGATE's count-scale
#                    SifiNet pipeline. [default]
#   mode "lognorm" : log1p(CP10K) -- Seurat standard log-scale data slot.
# Usage: lognorm_export.R <raw_counts_csv> <out_csv> [mode]
suppressMessages({ library(Seurat); library(Matrix); library(data.table) })
a <- commandArgs(trailingOnly = TRUE); in_csv <- a[1]; out_csv <- a[2]
mode <- ifelse(length(a) >= 3, a[3], "cp10k")
dt <- fread(in_csv); genes <- as.character(dt[[1]])
m <- as.matrix(dt[, -1]); rownames(m) <- genes
keep <- !is.na(genes) & nzchar(trimws(genes)) & !duplicated(genes)
m <- m[keep, , drop = FALSE]
storage.mode(m) <- "double"; m <- round(m); m[m < 0] <- 0   # RSEM expected counts -> integer
cat(sprintf("raw: %d genes x %d cells\n", nrow(m), ncol(m)))
obj <- CreateSeuratObject(counts = as(m, "dgCMatrix"))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
ln <- as.matrix(GetAssayData(obj, layer = "data"))          # log1p(CP10K)
gkeep <- rowSums(ln > 0) >= 0.05 * ncol(ln)
ckeep <- colSums(ln > 0) >= 0.05 * nrow(ln)
ln <- ln[gkeep, ckeep, drop = FALSE]
cat(sprintf("LogNorm (filtered): %d genes x %d cells\n", nrow(ln), ncol(ln)))
out <- data.table(gene = rownames(ln), as.data.table(ln))
fwrite(out, out_csv)
cat("wrote", out_csv, "\n")
