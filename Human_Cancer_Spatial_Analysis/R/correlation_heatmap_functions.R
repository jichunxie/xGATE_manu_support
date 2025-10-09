# R/correlation_heatmap_functions.R
#
# Contains functions for building a data matrix from pathway results,
# shortening pathway labels for display, and plotting a correlation heatmap
# using the ComplexHeatmap package.

suppressPackageStartupMessages({
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


#' Shorten Pathway Labels for Display
#'
#' Cleans up long pathway names by removing common suffixes, applying standard
#' abbreviations, and optionally truncating to a maximum character length.
#'
#' @param labels A character vector of pathway labels.
#' @param max_chars An integer specifying the maximum number of characters.
#'   Longer labels will be truncated with an ellipsis. Set to Inf to disable.
#' @return A character vector of shortened labels.
shorten_pathway_labels <- function(labels, max_chars = 40) {
  x <- as.character(labels)

  # Normalize separators and whitespace
  x <- gsub("[\\p{Pd}_/\\.]+", " ", x, perl = TRUE)
  x <- stringr::str_squish(x)

  # Remove trailing 'pathway' or 'signature'
  x <- sub("\\s*\\b(signaling pathway|signature)\\b\\s*$", "", x, ignore.case = TRUE)
  x <- stringr::str_squish(x)

  # Specific phrase replacements (case-insensitive)
  x <- gsub("\\bEGFR tyrosine kinase inhibitor resistance\\b", "EGFR TKI Resistance", x, ignore.case = TRUE)
  x <- gsub("\\bJAK-STAT signaling pathway\\b", "JAK-STAT", x, ignore.case = TRUE)

  # Common global abbreviations inside strings (case-insensitive)
  x <- gsub("nf\\s*kappa\\s*b", "NF-κB", x, ignore.case = TRUE)
  x <- gsub("\\btyrosine kinase inhibitor\\b", "TKI", x, ignore.case = TRUE)
  x <- stringr::str_squish(x)

  # Optional max char truncation with smart word break + ellipsis
  if (is.finite(max_chars) && max_chars > 3) {
    x <- vapply(x, function(s) {
      if (nchar(s) <= max_chars) return(s)
      # Trim to just under max_chars and find the last space
      trimmed <- substr(s, 1, max_chars - 1)
      last_space <- rev(gregexpr(" ", trimmed)[[1]])[1]
      if (!is.na(last_space) && last_space > 1) {
        trimmed <- substr(trimmed, 1, last_space - 1)
      }
      paste0(trimmed, "…")
    }, character(1))
  }
  return(x)
}


#' Build a Pathways x Subdomains Matrix
#'
#' Pivots a long-format data frame of pathway results into a wide matrix where
#' rows are pathways and columns are subdomains.
#'
#' @param df The input data frame.
#' @param pathway_col The name of the column containing pathway names.
#' @param subdomain_key A vector of unique identifiers for each subdomain.
#' @param value_vec The numeric vector of values (e.g., z-scores) to fill the matrix.
#' @param min_non_na_cols Minimum number of non-NA values for a pathway to be kept.
#' @return A numeric matrix with pathways as rows and subdomains as columns.
build_matrix <- function(df, pathway_col, subdomain_key, value_vec,
                         what = "value", min_non_na_cols = 2) {

  df_use <- tibble::tibble(
    pathway   = as.character(df[[pathway_col]]),
    subdomain = subdomain_key,
    value     = as.numeric(value_vec)
  ) %>%
    dplyr::filter(!is.na(pathway) & nzchar(pathway),
                  !is.na(subdomain) & nzchar(subdomain))

  # Average any duplicates (though ideally there are none)
  mat_df <- df_use %>%
    dplyr::group_by(pathway, subdomain) %>%
    dplyr::summarize(score = mean(value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = subdomain, values_from = score)

  pathways <- mat_df$pathway
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- pathways

  # Filter by minimal non-NA values
  if (is.finite(min_non_na_cols) && min_non_na_cols > 0) {
    keep_rows <- rowSums(is.finite(mat)) >= min_non_na_cols
    if (any(!keep_rows)) {
      message(
        "Removing ", sum(!keep_rows), " pathways with fewer than ",
        min_non_na_cols, " non-NA values for ", what, "."
      )
    }
    mat <- mat[keep_rows, , drop = FALSE]
  }
  return(mat)
}


#' Plot and Save a Correlation Heatmap
#'
#' Generates and saves a correlation heatmap with specific ordering and formatted labels.
#' This function is self-contained and takes all necessary configuration as arguments.
#'
#' @param cor_mat The square correlation matrix to plot.
#' @param output_dir Directory to save the output JPG file.
#' @param tag A string tag (e.g., "zscore") used in the output filename.
#' @param corr_method The correlation method name (e.g., "pearson") for the legend title.
#' @param selected_pathways A character vector of pathways in the desired order for the heatmap.
#' @param case_insensitive_match Logical, whether to match `selected_pathways` case-insensitively.
#' @param label_max_chars Max character length for labels before truncation.
#' @param output_width_in Width of the output JPG in inches.
#' @param output_height_in Height of the output JPG in inches.
#' @return Invisibly returns NULL. The plot is saved to a file.
plot_and_save_heatmap <- function(cor_mat,
                                  output_dir,
                                  tag,
                                  corr_method,
                                  selected_pathways,
                                  case_insensitive_match = TRUE,
                                  label_max_chars = 50,
                                  output_width_in = 12,
                                  output_height_in = 12) {
  # Guard against invalid correlation values
  cor_mat[cor_mat > 1] <- 1
  cor_mat[cor_mat < -1] <- -1

  # ---- SUBSET AND REORDER MATRIX ----
  # Reorder the matrix to match the exact order specified in `selected_pathways`
  if (isTRUE(case_insensitive_match)) {
    mat_rownames_lower <- tolower(rownames(cor_mat))
    sel_pathways_lower <- tolower(selected_pathways)
    
    # Find the indices in the matrix that correspond to the selected pathways
    row_indices <- match(sel_pathways_lower, mat_rownames_lower)
    row_indices <- row_indices[!is.na(row_indices)] # Remove not found
  } else {
    row_indices <- match(selected_pathways, rownames(cor_mat))
    row_indices <- row_indices[!is.na(row_indices)] # Remove not found
  }
  
  if (length(row_indices) < 2) {
    warning("Fewer than 2 of the selected pathways were found in the correlation matrix. Skipping heatmap for tag: ", tag)
    return(invisible(NULL))
  }
  
  # Get the exact names from the matrix in the desired order
  ordered_pathway_names <- rownames(cor_mat)[row_indices]
  cor_mat_ordered <- cor_mat[ordered_pathway_names, ordered_pathway_names, drop = FALSE]

  # ---- PREPARE LABELS AND STYLING ----
  row_labs <- shorten_pathway_labels(rownames(cor_mat_ordered), max_chars = label_max_chars)
  col_labs <- shorten_pathway_labels(colnames(cor_mat_ordered), max_chars = label_max_chars)

  row_gp <- grid::gpar(fontsize = 10)
  col_gp <- grid::gpar(fontsize = 10)

  # ---- DEFINE HEATMAP ----
  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#3b4cc0", "#f7f7f7", "#b40426"))

  ht <- ComplexHeatmap::Heatmap(
    cor_mat_ordered,
    name = paste0(corr_method, "\ncorr."),
    col = col_fun,
    na_col = "grey90",
    cluster_rows = FALSE,      # Use the pre-defined order
    cluster_columns = FALSE,   # Use the pre-defined order
    row_labels = row_labs,
    column_labels = col_labs,
    row_names_gp = row_gp,
    column_names_gp = col_gp,
    row_names_side = "left",
    column_names_side = "bottom",
    column_names_rot = 45, # Angle column labels for readability
    use_raster = TRUE,
    raster_device = "png"
  )

  # ---- SAVE PLOT TO FILE ----
  out_jpg <- file.path(output_dir, paste0("pathway_correlation_heatmap_", corr_method, "_from_", tag, ".jpg"))

  jpeg(
    filename = out_jpg,
    width = output_width_in,
    height = output_height_in,
    units = "in",
    res = 300,
    quality = 100
  )
  
  # Draw the heatmap with some padding around the edges
  ComplexHeatmap::draw(ht, padding = unit(c(5, 5, 5, 15), "mm")) # T, R, B, L padding
  
  dev.off()
  message("Saved heatmap: ", out_jpg)
  
  return(invisible(NULL))
}