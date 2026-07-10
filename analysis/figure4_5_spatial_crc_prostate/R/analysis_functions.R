# R/analysis_functions.R
#
# Contains functions for performing specific statistical analyses, such as
# calculating log odds ratios for cell type enrichment.

#' Compute Log2 Odds Ratio for Cell Type Enrichment in a Domain
#'
#' Calculates the log2 odds ratio of a specific cell type being found inside
#' a given domain versus outside of it.
#'
#' @param domain_id The identifier for the domain of interest.
#' @param cell_type_name The name of the cell type of interest.
#' @param metadata_df A data frame with at least 'domain' and 'general_cell_type' columns.
#' @return The calculated log2 odds ratio as a single numeric value.
compute_log_odds <- function(domain_id, cell_type_name, metadata_df) {
  # Counts for the 2x2 contingency table
  a <- sum(metadata_df$domain == domain_id & metadata_df$general_cell_type == cell_type_name) # In domain, is cell type
  b <- sum(metadata_df$domain == domain_id & metadata_df$general_cell_type != cell_type_name) # In domain, not cell type
  c <- sum(metadata_df$domain != domain_id & metadata_df$general_cell_type == cell_type_name) # Out of domain, is cell type
  d <- sum(metadata_df$domain != domain_id & metadata_df$general_cell_type != cell_type_name) # Out of domain, not cell type
  
  # Add 0.5 pseudocount to prevent division by zero
  a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  
  # Log2 Odds Ratio: log2( (a/b) / (c/d) )
  log2_odds <- log2((a * d) / (b * c))
  
  return(log2_odds)
}