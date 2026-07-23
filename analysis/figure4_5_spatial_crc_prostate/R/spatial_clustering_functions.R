# R/spatial_clustering_functions.R
#
# Contains functions for spatial clustering, such as partitioning large
# domains into smaller, balanced subdomains.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Cluster Large Domains into Balanced Subdomains
#'
#' For each spatial domain in a data frame, this function checks if its size
#' exceeds a minimum threshold. If it does, it partitions the domain into
#' roughly equal-sized subdomains using k-means clustering on spot coordinates.
#'
#' @param loc_df A spatial location data frame containing 'x', 'y', and 'domain' columns.
#' @param min_size The minimum number of spots a domain can have before it is
#'   partitioned into smaller subdomains.
#' @return The input data frame with a new 'subdomain' column added.
cluster_domains_balanced <- function(loc_df, min_size = 5000) {
  loc_df <- loc_df %>% mutate(subdomain = NA_character_)

  # Get a table of domain sizes
  domain_sizes <- table(loc_df$domain)

  for (dom in names(domain_sizes)) {
    # Find all spots belonging to the current domain
    domain_indices <- which(loc_df$domain == dom)
    coords <- loc_df[domain_indices, c("x", "y")]
    size <- domain_sizes[[dom]]

    if (size > min_size) {
      # If the domain is large, cluster it
      n_clusters <- ceiling(size / min_size)
      km_result <- kmeans(coords, centers = n_clusters, nstart = 20)
      loc_df$subdomain[domain_indices] <- paste0(dom, "_", km_result$cluster)
    } else {
      # If the domain is small, it becomes a single subdomain
      loc_df$subdomain[domain_indices] <- paste0(dom, "_1")
    }
  }

  return(loc_df)
}