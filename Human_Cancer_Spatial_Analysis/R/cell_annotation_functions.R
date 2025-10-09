# R/cell_annotation_functions.R
#
# Contains functions for performing cell type annotation on spatial data
# using a single-cell reference dataset with SingleR.

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleR)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(BiocParallel)
})

# --- Helper Functions ---
# (These are internal to the main function below)

# Heuristic to check if a matrix contains raw counts
.is_count_like <- function(m, sample_n = 100000, tol = 1e-8) {
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  nz <- m@x
  if (length(nz) == 0) return(TRUE)
  if (length(nz) > sample_n) nz <- sample(nz, sample_n)
  mean(abs(nz - round(nz)) < tol) >= 0.95
}

# Robustly build a spatial SingleCellExperiment object
.make_spatial_sce <- function(sp_obj) {
  stopifnot(is.list(sp_obj), !is.null(sp_obj$spatial_countMat_list))
  if (length(sp_obj$spatial_countMat_list) != 1L) {
    stop("This annotator expects a single slice in spatial_countMat_list.")
  }

  mat <- sp_obj$spatial_countMat_list[[1]]
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")

  # Spatial locations to align colnames
  loc <- sp_obj$spatial_location_list[[1]]
  if (is.null(rownames(loc))) stop("spatial_location_list[[1]] must have rownames (cell IDs).")
  loc_ids <- rownames(loc)

  # Ensure colnames exist; if absent, set from loc_ids (must match ncol)
  if (is.null(colnames(mat))) {
    if (ncol(mat) != length(loc_ids)) {
      stop("Normalized/raw matrix has no colnames and ncol != number of spatial locations.")
    }
    colnames(mat) <- loc_ids
  }

  # If there is a mismatch, intersect and reorder to loc order
  common_cells <- intersect(loc_ids, colnames(mat))
  if (length(common_cells) == 0) {
    stop("No overlapping cell IDs between matrix columns and spatial locations.")
  }
  # Restrict both to common set
  loc_ids <- common_cells
  mat <- mat[, loc_ids, drop = FALSE]

  # Make column names unique if duplicated
  if (anyDuplicated(colnames(mat))) {
    warning("Duplicated cell IDs detected in matrix columns; making them unique.")
    colnames(mat) <- make.unique(colnames(mat))
    # Keep locations in same order and names
    rownames(loc) <- make.unique(rownames(loc))
  }

  # Ensure rownames (gene IDs) exist; if missing, create synthetic (discouraged)
  if (is.null(rownames(mat))) {
    warning("Matrix has no rownames (genes). Creating synthetic gene IDs (g1..gN). Consider adding real gene IDs.")
    rownames(mat) <- sprintf("g%d", seq_len(nrow(mat)))
  }

  # Make rownames unique if duplicated (avoid assay dimnames conflicts)
  if (anyDuplicated(rownames(mat))) {
    warning("Duplicated gene IDs detected; making them unique with suffixes.")
    rownames(mat) <- make.unique(rownames(mat))
  }

  if (.is_count_like(mat)) {
    sce <- SingleCellExperiment(list(counts = mat))
    sce <- scuttle::logNormCounts(sce)
    attr(sce, "assay.type.test") <- "logcounts"
  } else {
    # Already log-like; put into logcounts directly via constructor
    sce <- SingleCellExperiment(list(logcounts = mat))
    attr(sce, "assay.type.test") <- "logcounts"
  }

  sce
}

# Convert a Seurat object to a SingleCellExperiment reference
.make_ref_from_seurat <- function(seurat_rds_path,
                                  label_col,
                                  assay = "RNA",
                                  slot = "counts",
                                  min_cells_per_label = 25L,
                                  verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required to use a Seurat reference.")
  }
  seu <- readRDS(seurat_rds_path)
  if (!label_col %in% colnames(seu[[]])) {
    stop("Label column '", label_col, "' not found in Seurat meta.data. Available: ",
         paste(colnames(seu[[]]), collapse = ", "))
  }
  mat <- Seurat::GetAssayData(seu, assay = assay, slot = slot)
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")

  # Ensure dimnames
  if (is.null(rownames(mat))) stop("Reference matrix (Seurat) has no gene rownames.")
  if (is.null(colnames(mat))) {
    colnames(mat) <- colnames(seu)
    if (is.null(colnames(mat))) stop("Reference matrix has no column names and could not infer from Seurat.")
  }
  if (anyDuplicated(rownames(mat))) {
    warning("Duplicated gene IDs in reference; making them unique.")
    rownames(mat) <- make.unique(rownames(mat))
  }
  if (anyDuplicated(colnames(mat))) {
    warning("Duplicated cell IDs in reference; making them unique.")
    colnames(mat) <- make.unique(colnames(mat))
  }

  ref_sce <- SingleCellExperiment(list(counts = mat))
  SummarizedExperiment::colData(ref_sce)$label <- as.character(seu[[label_col]][, 1])

  # Drop empty and rare labels
  keep <- nzchar(SummarizedExperiment::colData(ref_sce)$label)
  ref_sce <- ref_sce[, keep, drop = FALSE]
  if (min_cells_per_label > 0) {
    tab <- table(SummarizedExperiment::colData(ref_sce)$label)
    keep_labs <- names(tab)[tab >= min_cells_per_label]
    dropped <- setdiff(names(tab), keep_labs)
    if (isTRUE(verbose) && length(dropped)) {
      message("Dropping rare labels (<", min_cells_per_label, " cells): ", paste(dropped, collapse = ", "))
    }
    ref_sce <- ref_sce[, SummarizedExperiment::colData(ref_sce)$label %in% keep_labs, drop = FALSE]
  }

  ref_sce <- scuttle::logNormCounts(ref_sce)
  if (isTRUE(verbose)) {
    message("Reference (Seurat) SCE: genes=", nrow(ref_sce), " cells=", ncol(ref_sce),
            " labels=", length(unique(SummarizedExperiment::colData(ref_sce)$label)))
  }
  ref_sce
}

# Prepare a generic SingleCellExperiment reference
.prepare_ref_sce <- function(ref_sce, label_col = "label", min_cells_per_label = 25L, verbose = TRUE) {
  stopifnot(inherits(ref_sce, "SingleCellExperiment"))

  if (!label_col %in% colnames(SummarizedExperiment::colData(ref_sce))) {
    stop("label_col '", label_col, "' not found in colData(ref_sce). Available: ",
         paste(colnames(SummarizedExperiment::colData(ref_sce)), collapse = ", "))
  }
  SummarizedExperiment::colData(ref_sce)$label <- as.character(SummarizedExperiment::colData(ref_sce)[[label_col]])
  keep <- nzchar(SummarizedExperiment::colData(ref_sce)$label)
  ref_sce <- ref_sce[, keep, drop = FALSE]

  if (min_cells_per_label > 0) {
    tab <- table(SummarizedExperiment::colData(ref_sce)$label)
    keep_labs <- names(tab)[tab >= min_cells_per_label]
    dropped <- setdiff(names(tab), keep_labs)
    if (isTRUE(verbose) && length(dropped)) {
      message("Dropping rare labels (<", min_cells_per_label, " cells): ", paste(dropped, collapse = ", "))
    }
    ref_sce <- ref_sce[, SummarizedExperiment::colData(ref_sce)$label %in% keep_labs, drop = FALSE]
  }

  # Ensure logcounts present
  if (!"logcounts" %in% assayNames(ref_sce)) {
    if ("counts" %in% assayNames(ref_sce)) {
      ref_sce <- scuttle::logNormCounts(ref_sce)
    } else {
      stop("Reference SCE has neither 'logcounts' nor 'counts' assays.")
    }
  }

  # Ensure names are unique
  if (anyDuplicated(rownames(ref_sce))) {
    warning("Duplicated gene IDs in reference SCE; making them unique.")
    rownames(ref_sce) <- make.unique(rownames(ref_sce))
  }

  if (isTRUE(verbose)) {
    message("Reference (SCE) prepared: genes=", nrow(ref_sce), " cells=", ncol(ref_sce),
            " labels=", length(unique(SummarizedExperiment::colData(ref_sce)$label)))
  }
  ref_sce
}

# Aggregate a reference to pseudo-bulk profiles
.aggregate_ref <- function(ref_sce, verbose = TRUE) {
  aggr <- SingleR::aggregateReference(ref_sce, labels = SummarizedExperiment::colData(ref_sce)$label)
  if (isTRUE(verbose)) {
    labs <- as.character(SummarizedExperiment::colData(aggr)$label)
    message("Aggregated reference: ", length(labs), " pseudo-bulk profiles for ",
            length(unique(labs)), " labels.")
  }
  aggr
}

# Align genes between the test and reference datasets
.align_features <- function(test_sce, ref_se) {
  common <- intersect(rownames(test_sce), rownames(ref_se))
  if (length(common) < 200) {
    warning("Low gene overlap (", length(common), "). Consider harmonizing gene IDs (symbols vs Ensembl).")
  }
  test_sce <- test_sce[common, , drop = FALSE]
  ref_se   <- ref_se[common, , drop = FALSE]
  list(test = test_sce, ref = ref_se)
}


# --- Main Annotation Function ---

#' Annotate Spatial Data with a Single-Cell Reference
#'
#' Uses SingleR to perform cell type annotation on a spatial dataset (`sp_input` format)
#' using either a Seurat or SingleCellExperiment object as a reference.
#'
#' @param spatial_rds Path to the spatial data object to be annotated.
#' @param out_rds Path where the final annotated spatial object will be saved.
#' @param ref_seurat_rds Path to a Seurat object (`.rds`) to use as a reference.
#' @param ref_seurat_label_col The column in the Seurat metadata with cell type labels.
#' @param ref_sce_rds Path to a SingleCellExperiment object (`.rds`) to use as a reference.
#' @param ref_sce_label_col The column in the SCE colData with cell type labels.
#' @param min_cells_per_label Minimum number of cells required for a label to be kept in the reference.
#' @param aggregate_reference Logical, whether to aggregate the reference to pseudo-bulk profiles.
#' @param batch_size The number of cells to process in each chunk.
#' @param workers The number of parallel workers to use.
#' @param verbose Logical, whether to print progress messages.
#' @return Invisibly returns the annotated spatial object.
annotate_spatial_with_singleR <- function(
    spatial_rds,
    out_rds,
    ref_seurat_rds = NULL,
    ref_seurat_label_col = NULL,
    ref_sce_rds = NULL,
    ref_sce_label_col = "label",
    min_cells_per_label = 25L,
    aggregate_reference = TRUE,
    batch_size = 50000,
    workers = NULL,
    verbose = TRUE
) {
  sp_obj <- readRDS(spatial_rds)

  # Build spatial SCE with correct dimnames handling
  test_sce <- .make_spatial_sce(sp_obj)
  assay_test <- attr(test_sce, "assay.type.test")
  if (is.null(assay_test)) assay_test <- "logcounts"

  # Build reference
  if (!is.null(ref_seurat_rds)) {
    if (is.null(ref_seurat_label_col)) stop("Provide ref_seurat_label_col for Seurat reference.")
    ref_sce <- .make_ref_from_seurat(
      seurat_rds_path = ref_seurat_rds,
      label_col = ref_seurat_label_col,
      min_cells_per_label = min_cells_per_label,
      verbose = verbose
    )
  } else if (!is.null(ref_sce_rds)) {
    ref_sce <- readRDS(ref_sce_rds)
    ref_sce <- .prepare_ref_sce(
      ref_sce, label_col = ref_sce_label_col,
      min_cells_per_label = min_cells_per_label, verbose = verbose
    )
  } else {
    stop("Provide either ref_seurat_rds+ref_seurat_label_col or ref_sce_rds+ref_sce_label_col.")
  }

  ref_final <- if (aggregate_reference) .aggregate_ref(ref_sce, verbose = verbose) else ref_sce

  # Align genes
  aligned <- .align_features(test_sce, ref_final)
  test_sce <- aligned$test
  ref_final <- aligned$ref

  # Labels vector
  cd <- SummarizedExperiment::colData(ref_final)
  labels_vec <- if ("label.main" %in% colnames(cd)) cd$label.main else {
    if ("label" %in% colnames(cd)) cd$label else stop("Reference colData lacks 'label' or 'label.main'.")
  }

  # Parallel backend
  BPPARAM <- if (is.null(workers)) BiocParallel::bpparam() else {
    if (.Platform$OS.type == "unix") BiocParallel::MulticoreParam(workers = as.integer(workers))
    else BiocParallel::SnowParam(workers = as.integer(workers))
  }

  # Chunk and run SingleR
  n <- ncol(test_sce)
  idx_chunks <- split(seq_len(n), ceiling(seq_len(n) / batch_size))
  if (isTRUE(verbose)) message("Annotating ", n, " cells in ", length(idx_chunks), " chunk(s).")

  run_chunk <- function(cols) {
    pred <- SingleR::SingleR(
      test = test_sce[, cols, drop = FALSE],
      ref = ref_final,
      labels = labels_vec,
      assay.type.test = assay_test,
      assay.type.ref = if ("logcounts" %in% assayNames(ref_final)) "logcounts" else "counts",
      BPPARAM = BPPARAM
    )
    labs <- if ("pruned.labels" %in% colnames(pred)) pred$pruned.labels else pred$labels
    data.frame(cell_id = colnames(test_sce)[cols],
               general_cell_type = as.character(labs),
               stringsAsFactors = FALSE)
  }
  ann_list <- lapply(idx_chunks, run_chunk)
  anno <- do.call(rbind, ann_list)

  # Attach labels to spatial locations
  loc <- sp_obj$spatial_location_list[[1]]
  loc_df <- as.data.frame(loc)
  if (!"cell_id" %in% colnames(loc_df)) loc_df$cell_id <- rownames(loc_df)
  loc_df <- merge(loc_df, anno, by = "cell_id", all.x = TRUE, sort = FALSE)
  rownames(loc_df) <- loc_df$cell_id
  sp_obj$spatial_location_list[[1]] <- loc_df

  # Summary
  if (isTRUE(verbose)) {
    tab <- sort(table(loc_df$general_cell_type), decreasing = TRUE)
    message("Top assigned types:")
    print(utils::head(tab, 20))
    n_ok <- sum(!is.na(loc_df$general_cell_type))
    message(sprintf("Annotated %d / %d cells (%.1f%%).", n_ok, nrow(loc_df), 100 * n_ok / max(1, nrow(loc_df))))
  }

  saveRDS(sp_obj, out_rds)
  message("Saved annotated object: ", out_rds)
  invisible(sp_obj)
}