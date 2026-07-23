# R/data_ingestion_functions.R
#
# Contains functions for reading and processing raw spatial transcriptomics
# data from various formats (e.g., 10x Genomics output) into the list
# structure required by downstream tools like IRIS.

suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
  library(dplyr)
  library(stringr)
})

#' Build a Single-Slice `sp_input` Object from 10x-like Files
#'
#' Reads raw matrix, barcode, feature, and spatial location files to construct
#' the standard `sp_input` list format used for IRIS analysis.
#'
#' @param matrix_fp Path to the matrix.mtx[.gz] file.
#' @param barcodes_fp Path to the barcodes.tsv[.gz] file.
#' @param features_fp Path to the features.tsv[.gz] file.
#' @param spatial_loc_fp Path to the spatial coordinates CSV file.
#' @param save_rds Optional. If a file path is provided, the resulting
#'   `sp_input` object will be saved to this location.
#' @return A list containing `spatial_countMat_list` and `spatial_location_list`.
build_sp_input_slice <- function(matrix_fp, barcodes_fp, features_fp, spatial_loc_fp,
                                 save_rds = NULL) {
  # --- 1. Read count matrix and feature/barcode IDs ---
  mat <- Matrix::readMM(matrix_fp)
  barcodes <- readr::read_tsv(barcodes_fp, col_names = FALSE, show_col_types = FALSE)[[1]]
  feats <- readr::read_tsv(features_fp, col_names = FALSE, show_col_types = FALSE)
  
  # Assume standard 3-column features file (ID, Name, Type)
  gene_names <- feats[[2]]
  
  # Ensure matrix is genes x barcodes
  if (nrow(mat) == length(gene_names) && ncol(mat) == length(barcodes)) {
    # Dimensions are correct
  } else if (ncol(mat) == length(gene_names) && nrow(mat) == length(barcodes)) {
    mat <- Matrix::t(mat) # Transpose if needed
  } else {
    stop("Matrix dimensions are incompatible with feature/barcode files.")
  }
  rownames(mat) <- make.unique(gene_names)
  colnames(mat) <- barcodes

  # --- 2. Read and standardize spatial locations ---
  loc <- readr::read_csv(spatial_loc_fp, show_col_types = FALSE, progress = FALSE)
  
  # Check for standard 10x Visium/Xenium headers first
  if (all(c("barcode", "pxl_row_in_fullres", "pxl_col_in_fullres") %in% names(loc))) {
    loc_df <- loc %>%
      dplyr::transmute(
        cell_id = as.character(.data$barcode),
        x = .data$pxl_col_in_fullres,
        y = .data$pxl_row_in_fullres
      )
  } else if ("x" %in% names(loc) && "y" %in% names(loc) && "cell_id" %in% names(loc)) {
    loc_df <- loc %>% dplyr::select(cell_id, x, y)
  } else {
    stop("Could not find required columns in spatial location file. Need 'barcode', 'pxl_row_in_fullres', and 'pxl_col_in_fullres', OR 'cell_id', 'x', and 'y'.")
  }
  
  # --- 3. Align matrix and locations by common barcodes ---
  common_barcodes <- intersect(colnames(mat), loc_df$cell_id)
  if (length(common_barcodes) == 0) {
    stop("No overlapping barcodes found between the count matrix and spatial locations.")
  }
  
  mat_aligned <- mat[, common_barcodes]
  loc_aligned <- loc_df[match(common_barcodes, loc_df$cell_id), ]
  rownames(loc_aligned) <- loc_aligned$cell_id

  # --- 4. Construct the final `sp_input` object ---
  sp_input <- list(
    spatial_countMat_list = list(mat_aligned),
    spatial_location_list = list(as.data.frame(loc_aligned[, c("x", "y")]))
  )

  # --- 5. Save the output if requested ---
  if (!is.null(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(sp_input, save_rds)
    message("Saved sp_input object to: ", save_rds)
  }

  return(sp_input)
}