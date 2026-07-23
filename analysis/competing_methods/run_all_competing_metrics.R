# =============================================================================
# Comprehensive competing-method benchmark for the xGATE revision.
# Runs AUCell, ORA, scGSEA (gficf), GESECA (fgsea), PAGODA (pagoda2) on one dataset
# and computes ALL benchmark metrics (precision, recall, F1, specificity, accuracy,
# MCC, AUCPR) per method. ENSEMBL gene IDs throughout (gficf requires ensembl).
#
# ---------- HOW TO RUN IN RSTUDIO ----------
# 1. Open this file in RStudio. Make sure you're in your R 4.4 with your ~/R library
#    (the one that loads AUCell/pagoda2/scde/fgsea/GSVA).
# 2. One-time: install gficf if missing (you have rtracklayer/rhdf5/msigdbr already):
#       remotes::install_github("gambalab/gficf")
# 3. Source the whole file (Ctrl+Shift+S, or click "Source"). This loads the libraries
#    and defines run_benchmark().
# 4. In the Console, run it for each dataset:
#
#    run_benchmark(
#      counts_csv = "/path/to/xGATE/data/fucci_u2os/fucci_counts_ensembl.csv",
#      gs_json    = "/path/to/xGATE/data/fucci_u2os/fucci_genesets_ensembl.json",
#      out_prefix = "/path/to/xGATE/results/fucci")
#
#    run_benchmark(
#      counts_csv = "/path/to/xGATE/data/ts_fibroblast/ts_counts_ensembl.csv",
#      gs_json    = "/path/to/xGATE/data/ts_fibroblast/ts_genesets_ensembl.json",
#      out_prefix = "/path/to/xGATE/results/ts")
#
# Outputs -> <out_prefix>_competing_metrics.csv  and  <out_prefix>_competing_percall.csv
# =============================================================================

suppressMessages({
  library(jsonlite); library(Matrix); library(data.table)
  library(AUCell); library(clusterProfiler); library(fgsea); library(pagoda2)
})

## area under precision-recall curve from a continuous score
auprc <- function(sc, y) {
  sc <- as.numeric(sc); ok <- is.finite(sc); sc <- sc[ok]; y <- y[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  o <- order(sc, decreasing = TRUE); y <- y[o]
  tp <- cumsum(y); fp <- cumsum(!y); rec <- tp/sum(y); prec <- tp/(tp+fp)
  sum(diff(c(0, rec)) * prec)
}

run_benchmark <- function(counts_csv, gs_json, out_prefix, aucell_cutoff = 0.20) {
  message("== load: ", counts_csv)
  dt <- fread(counts_csv); genes <- as.character(dt[[1]])
  counts <- as.matrix(dt[, -1]); rownames(counts) <- genes; storage.mode(counts) <- "double"
  counts <- counts[!duplicated(rownames(counts)), ]
  message(sprintf("   counts: %d genes x %d cells", nrow(counts), ncol(counts)))
  gs <- fromJSON(gs_json, simplifyVector = FALSE)
  labels <- sapply(gs, function(x) x$label); detected <- rownames(counts)
  genesets <- lapply(gs, function(x) intersect(as.character(unlist(x$genes)), detected))
  genesets <- genesets[sapply(genesets, length) >= 5]
  paths <- names(genesets); labels <- labels[paths]; is_pos <- labels == "positive"
  message(sprintf("   %d gene sets (%d pos / %d neg)", length(paths), sum(is_pos), sum(!is_pos)))
  calls <- list(); scores <- list()

  message("== AUCell =="); try({
    ar <- AUCell_buildRankings(counts, plotStats = FALSE, splitByBlocks = TRUE, verbose = FALSE)
    auc <- AUCell_calcAUC(genesets, ar, verbose = FALSE)
    th  <- AUCell_exploreThresholds(auc, plotHist = FALSE, assignCells = TRUE, verbose = FALSE)
    frac <- sapply(paths, function(p) length(th[[p]]$assignment) / ncol(counts))
    for (q in c(0.10,0.20,0.30,0.50)) calls[[sprintf("AUCell_%d", q*100)]] <- frac[paths] >= q
    calls[["AUCell"]] <- frac[paths] >= aucell_cutoff; scores[["AUCell"]] <- frac[paths]
  })

  message("== ORA =="); try({
    pb <- rowMeans(counts); top <- names(sort(pb, decreasing = TRUE))[1:round(0.10*length(pb))]
    t2g <- rbindlist(lapply(paths, function(p) data.table(term = p, gene = genesets[[p]])))
    ora <- as.data.frame(enricher(top, TERM2GENE = t2g, universe = detected,
              pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 5000))
    op <- setNames(rep(1, length(paths)), paths); if (nrow(ora)) op[ora$ID] <- ora$p.adjust
    calls[["ORA"]] <- op[paths] < 0.05; scores[["ORA"]] <- 1 - op[paths]
  })

  message("== scGSEA (gficf) ==")
  if (requireNamespace("gficf", quietly = TRUE)) { try({
    library(gficf)
    data <- gficf::gficf(M = as(counts, "dgCMatrix"), verbose = FALSE)
    # if this errors on geneID, switch "ensembl" -> "ensgene" (USER used "ensamble")
    data <- gficf::runScGSEA(data = data, geneID = "ensamble", species = "human",
              pathway.list = genesets, nmf.k = 100, fdr.th = 0.1, rescale = "none", verbose = FALSE)
    sx <- as.matrix(data[["scgsea"]][["x"]])
    med <- apply(sx, 2, median); med <- med[intersect(names(med), paths)]
    sc_full <- setNames(rep(NA_real_, length(paths)), paths); sc_full[names(med)] <- med
    scores[["scGSEA"]] <- sc_full
    calls[["scGSEA_fixed"]] <- sc_full > median(sc_full, na.rm = TRUE)
    f1_at <- function(thr){ c<-sc_full>thr; tp<-sum(c&is_pos,na.rm=T); fp<-sum(c&!is_pos,na.rm=T); fn<-sum(!c&is_pos,na.rm=T)
      p<-ifelse(tp+fp>0,tp/(tp+fp),0); r<-ifelse(tp+fn>0,tp/(tp+fn),0); ifelse(p+r>0,2*p*r/(p+r),0) }
    grid <- sort(unique(sc_full[is.finite(sc_full)])); best <- grid[which.max(sapply(grid, f1_at))]
    calls[["scGSEA"]] <- sc_full > best
  }) } else message("   gficf NOT installed -> remotes::install_github('gambalab/gficf')")

  message("== GESECA =="); try({
    logmat <- log1p(sweep(counts, 2, pmax(colSums(counts),1), "/") * 1e4)
    ges <- geseca(genesets, logmat, minSize = 5, maxSize = 5000)
    gp <- setNames(rep(1, length(paths)), paths); if (nrow(ges)) gp[ges$pathway] <- ges$padj
    calls[["GESECA"]] <- gp[paths] < 0.05; scores[["GESECA"]] <- 1 - gp[paths]
  })

  message("== PAGODA =="); try({
    set.seed(1234); p2 <- Pagoda2$new(as(counts, "dgCMatrix"), log.scale = TRUE, n.cores = 1)
    p2$adjustVariance(plot = FALSE, gam.k = 10); p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3000)
    p2$testPathwayOverdispersion(setenv = list2env(genesets), verbose = FALSE,
                                 recalculate.pca = FALSE, min.pathway.size = 5)
    z <- setNames(rep(NA_real_, length(paths)), paths)
    for (pw in paths) if (!is.null(p2$misc$pwpca[[pw]])) z[pw] <- p2$misc$pwpca[[pw]]$z
    calls[["PAGODA"]] <- !is.na(z) & z > 1.96; scores[["PAGODA"]] <- z
  })

  metric_row <- function(m) {
    if (is.null(calls[[m]])) return(NULL)
    c <- as.logical(calls[[m]]); c[is.na(c)] <- FALSE
    tp<-sum(c&is_pos); fp<-sum(c&!is_pos); fn<-sum(!c&is_pos); tn<-sum(!c&!is_pos)
    P<-ifelse(tp+fp>0,tp/(tp+fp),NA); R<-ifelse(tp+fn>0,tp/(tp+fn),NA)
    F1<-ifelse(!is.na(P)&!is.na(R)&P+R>0,2*P*R/(P+R),0)
    spec<-ifelse(tn+fp>0,tn/(tn+fp),NA); acc<-(tp+tn)/length(c)
    d<-sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); mcc<-ifelse(d>0,(as.double(tp)*tn-as.double(fp)*fn)/d,NA)
    aup<-if(!is.null(scores[[m]])) auprc(scores[[m]], is_pos) else NA
    data.frame(method=m,precision=P,recall=R,F1=F1,specificity=spec,accuracy=acc,MCC=mcc,AUCPR=aup,TP=tp,FP=fp,FN=fn,TN=tn)
  }
  met <- do.call(rbind, lapply(c("AUCell","ORA","scGSEA","GESECA","PAGODA"), metric_row))
  write.csv(met, paste0(out_prefix,"_competing_metrics.csv"), row.names = FALSE)
  pc <- data.frame(pathway=paths, label=labels, stringsAsFactors=FALSE)
  for(m in names(calls))  pc[[paste0(m,"_active")]] <- as.logical(calls[[m]])
  for(m in names(scores)) pc[[paste0(m,"_score")]]  <- as.numeric(scores[[m]])
  write.csv(pc, paste0(out_prefix,"_competing_percall.csv"), row.names = FALSE)
  message("\n== METRICS =="); print(met)
  message(sprintf("[done] -> %s_competing_metrics.csv + _competing_percall.csv", out_prefix))
  invisible(met)
}

## Rscript fallback (command line): run_all_competing_metrics.R counts gs out
.a <- commandArgs(trailingOnly = TRUE)
if (length(.a) >= 3) run_benchmark(.a[1], .a[2], .a[3])
