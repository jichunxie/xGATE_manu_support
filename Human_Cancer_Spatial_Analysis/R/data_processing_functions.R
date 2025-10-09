# R/data_processing_functions.R
#
# Contains functions for common data wrangling and preprocessing tasks.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

#' Preprocess and Combine Spatial and Pathway Data
#'
#' Loads spatial data and pathway analysis results, cleans column names,
#' parses domain/cluster information, and returns a list of processed data frames.
#'
#' @param sp_rds_path Path to the spatial data RDS file (with subdomains).
#' @param pathway_csv_path Path to the pathway results CSV file.
#' @return A list containing two data frames: `spatial_df` and `pathway_res`.
preprocess_spatial_pathway_data <- function(sp_rds_path, pathway_csv_path) {
  
  # --- 1. Load and process pathway results ---
  pathway_res <- read.csv(pathway_csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  names(pathway_res) <- make.names(names(pathway_res)) # "p-value" -> "p.value"
  
  # Ensure required columns exist
  required_cols <- c("pathway", "domain", "cluster", "p.value")
  if (!all(required_cols %in% names(pathway_res))) {
    stop("Pathway file is missing one or more required columns: ", paste(required_cols, collapse=", "))
  }
  
  pathway_res <- pathway_res %>%
    mutate(
      domain_num  = as.integer(str_extract(as.character(domain), "\\d+")),
      cluster_num = as.integer(str_extract(as.character(cluster), "\\d+")),
      p.value     = as.numeric(p.value)
    ) %>%
    filter(!is.na(domain_num), !is.na(cluster_num))
  
  # --- 2. Load and process spatial data ---
  sp_obj <- readRDS(sp_rds_path)
  spatial_df <- as.data.frame(sp_obj$spatial_location_list[[1]])
  
  # The `.extract_domain_cluster` function should be in R/utils.R
  # For now, we define it locally if not found.
  if (!exists(".extract_domain_cluster")) {
    .extract_domain_cluster <- function(x) {
      nums <- stringr::str_extract_all(as.character(x), "\\d+")
      d <- suppressWarnings(as.integer(sapply(nums, `[`, 1)))
      c <- suppressWarnings(as.integer(sapply(nums, `[`, 2)))
      data.frame(domain_num = d, cluster_num = c)
    }
  }
  
  if ("subdomain" %in% names(spatial_df)) {
    parsed <- .extract_domain_cluster(spatial_df$subdomain)
    spatial_df$domain_num  <- parsed$domain_num
    spatial_df$cluster_num <- parsed$cluster_num
  } else {
    stop("Spatial data must contain a 'subdomain' column for this analysis.")
  }

  # Drop spots with no domain info
  spatial_df <- spatial_df %>% filter(!is.na(domain_num))

  return(list(spatial_df = spatial_df, pathway_res = pathway_res))
}