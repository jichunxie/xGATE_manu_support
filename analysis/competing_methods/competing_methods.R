#!/usr/bin/env Rscript
# Competing pathway-activity methods (AUCell, ORA, GESECA) on a benchmark, using the
# SAME entrez gene sets as xGATE. Args: counts_csv genesets_json out_csv
suppressMessages({
  library(AUCell); library(clusterProfiler); library(fgsea)
  library(jsonlite); library(data.table)
})
`%||%` <- function(a, b) if (!is.null(a)) a else b
args <- commandArgs(trailingOnly = TRUE)
counts_csv <- args[1]; gs_json <- args[2]; out_csv <- args[3]

message("[load] counts ", counts_csv)
dt <- fread(counts_csv)
genes <- as.character(dt[[1]]); counts <- as.matrix(dt[,-1]); rownames(counts) <- genes
storage.mode(counts) <- "double"
message("  counts: ", nrow(counts), " genes x ", ncol(counts), " cells")

gs <- fromJSON(gs_json, simplifyVector = FALSE)
labels <- sapply(gs, function(x) x$label)
detected <- rownames(counts)
genesets <- lapply(gs, function(x) intersect(as.character(unlist(x$genes)), detected))
genesets <- genesets[sapply(genesets, length) >= 5]
paths <- names(genesets)
message("  ", length(paths), " gene sets with >=5 detected genes")

## ---- AUCell (raw counts) ----
message("[AUCell] building rankings ...")
ar <- AUCell_buildRankings(counts, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)
auc <- AUCell_calcAUC(genesets, ar, verbose = FALSE)
am <- getAUC(auc)                                  # gene set x cell
frac_enriched <- function(cutoff_q = 0.5) {
  # a cell is 'enriched' for a set if its AUC exceeds the per-set assignment threshold
  apply(am, 1, function(a) mean(a > quantile(a, cutoff_q)))   # placeholder; replaced below
}
# AUCell assignment thresholds (auto) -> fraction of cells assigned per set
assign <- tryCatch({
  th <- AUCell_exploreThresholds(auc, plotHist = FALSE, assignCells = TRUE, verbose = FALSE)
  sapply(paths, function(p) length(th[[p]]$assignment) / ncol(counts))
}, error = function(e) { message("  threshold fallback: ", conditionMessage(e)); rep(NA, length(paths)) })
names(assign) <- paths

## ---- ORA (top-expressed genes -> hypergeometric via enricher) ----
message("[ORA] ...")
pb <- rowMeans(counts)
top <- names(sort(pb, decreasing = TRUE))[1:round(0.10 * length(pb))]
t2g <- rbindlist(lapply(paths, function(p) data.table(term = p, gene = genesets[[p]])))
ora <- tryCatch(as.data.frame(enricher(top, TERM2GENE = t2g, universe = detected,
        pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 3000)),
        error = function(e) data.frame(ID=character(), p.adjust=numeric()))
ora_p <- setNames(rep(1, length(paths)), paths); if (nrow(ora)) ora_p[ora$ID] <- ora$p.adjust

## ---- GESECA (log-normalized) ----
message("[GESECA] ...")
logmat <- log1p(sweep(counts, 2, pmax(colSums(counts), 1), "/") * 1e4)
ges <- tryCatch(geseca(genesets, logmat, minSize = 5, maxSize = 3000),
        error = function(e) { message("  geseca err: ", conditionMessage(e)); data.table(pathway=character(), padj=numeric()) })
ges_p <- setNames(rep(1, length(paths)), paths)
if (nrow(ges)) ges_p[ges$pathway] <- ges$padj

## ---- PAGODA (pagoda2 testPathwayOverdispersion; raw counts) ----
message("[PAGODA] ...")
pagoda_z <- setNames(rep(NA_real_, length(paths)), paths)
try({
  suppressMessages(library(pagoda2))
  cs <- as(counts, "dgCMatrix")
  penv <- list2env(genesets)
  set.seed(1234)
  p2 <- Pagoda2$new(cs, log.scale = TRUE, n.cores = 1)
  p2$adjustVariance(plot = FALSE, gam.k = 10)
  p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3000)
  p2$testPathwayOverdispersion(setenv = penv, verbose = FALSE,
                               recalculate.pca = FALSE, min.pathway.size = 5)
  for (pw in paths) if (!is.null(p2$misc$pwpca[[pw]]))
    pagoda_z[pw] <- if (!is.null(p2$misc$pwpca[[pw]]$z)) p2$misc$pwpca[[pw]]$z else NA_real_
})
pagoda_active <- !is.na(pagoda_z) & (pagoda_z > 1.96)

## ---- scGSEA (gficf::runScGSEA) if available ----
scgsea_active <- setNames(rep(NA, length(paths)), paths)
if (requireNamespace("gficf", quietly = TRUE)) {
  message("[scGSEA] ...")
  try({
    library(gficf)
    data <- gficf::gficf(M = counts, normalize = TRUE, verbose = FALSE)
    sg <- gficf::runScGSEA(data, geneID = "ensgene", pathway.list = genesets, verbose = FALSE)
    sc_scores <- sg$scgsea$x %||% sg$x
    med <- apply(sc_scores, 2, median)        # median per pathway
    cut <- median(med)                         # fixed cutoff (also report F1-opt downstream)
    for (pw in paths) if (pw %in% names(med)) scgsea_active[pw] <- med[pw] > cut
  })
} else message("[scGSEA] gficf not installed -> skipped")

out <- data.frame(
  pathway = paths, label = labels[paths], n_genes = sapply(genesets, length),
  AUCell_frac_enriched = assign[paths],
  AUCell_active_20 = assign[paths] >= 0.20,
  AUCell_active_10 = assign[paths] >= 0.10,
  AUCell_active_30 = assign[paths] >= 0.30,
  AUCell_active_50 = assign[paths] >= 0.50,
  ORA_padj = ora_p[paths], ORA_active = ora_p[paths] < 0.05,
  GESECA_padj = ges_p[paths], GESECA_active = ges_p[paths] < 0.05,
  PAGODA_z = pagoda_z[paths], PAGODA_active = pagoda_active[paths],
  scGSEA_active = scgsea_active[paths]
)
write.csv(out, out_csv, row.names = FALSE)
message("[done] -> ", out_csv)
print(out[, c("pathway","label","AUCell_frac_enriched","AUCell_active_20","ORA_active","GESECA_active")])
