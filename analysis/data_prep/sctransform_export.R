#!/usr/bin/env Rscript
# SCTransform (vst.flavor="v2") a raw genes x cells count CSV -> SCT corrected counts CSV
# (manuscript-standard preprocessing for xGATE + competing methods).
# Usage: sctransform_export.R <raw_counts_csv> <out_sct_csv>
suppressMessages({ library(Seurat); library(sctransform); library(Matrix); library(data.table) })
options(future.globals.maxSize = 8 * 1024^3)   # large donors exceed default 500 MiB future limit
a <- commandArgs(trailingOnly = TRUE); in_csv <- a[1]; out_csv <- a[2]
dt <- fread(in_csv); genes <- as.character(dt[[1]])
m <- as.matrix(dt[, -1]); rownames(m) <- genes
keep <- !is.na(genes) & nzchar(trimws(genes)) & !duplicated(genes)   # drop empty/NA/dup gene ids
m <- m[keep, , drop = FALSE]
storage.mode(m) <- "double"; m <- round(m)            # SCTransform expects integer counts
m[m < 0] <- 0
cat(sprintf("raw: %d genes x %d cells\n", nrow(m), ncol(m)))
obj <- CreateSeuratObject(counts = as(m, "dgCMatrix"))
obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE)
sct <- GetAssayData(obj, assay = "SCT", layer = "counts")   # corrected counts (as USER used)
sct <- as.matrix(sct)
# 5% / 5% filter (after normalization, per manuscript)
gkeep <- rowSums(sct > 0) >= 0.05 * ncol(sct)
ckeep <- colSums(sct > 0) >= 0.05 * nrow(sct)
sct <- sct[gkeep, ckeep, drop = FALSE]
cat(sprintf("SCT (filtered): %d genes x %d cells\n", nrow(sct), ncol(sct)))
out <- data.table(gene = rownames(sct), as.data.table(sct))
fwrite(out, out_csv)
cat("wrote", out_csv, "\n")
