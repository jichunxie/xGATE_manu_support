# --- 1. SETUP ---
# This script defines spatial subdomains by partitioning large domains
# into smaller, balanced clusters using k-means.
source("scripts/00_config.R")
source("R/spatial_clustering_functions.R")

# --- 2. LOAD PROCESSED DATA ---
message("Loading data with domain and cell type info: ", SUBDOMAIN_INPUT_FILE)
if (!file.exists(SUBDOMAIN_INPUT_FILE)) {
  stop("Input file not found. Please run previous processing scripts first (e.g., 01, 05).")
}
sp_input <- readRDS(SUBDOMAIN_INPUT_FILE)

# --- 3. CLUSTER DOMAINS TO CREATE SUBDOMAINS ---
message("Defining subdomains with a target minimum size of ", SUBDOMAIN_MIN_SIZE, "...")
# Apply the clustering function to each spatial slice in the list
sp_input$spatial_location_list <- lapply(
  sp_input$spatial_location_list,
  function(df) cluster_domains_balanced(df, min_size = SUBDOMAIN_MIN_SIZE)
)

# --- 4. SAVE DATA WITH SUBDOMAINS ---
saveRDS(sp_input, file = PROCESSED_SUBDOMAIN_FILE)
message("Processed data with subdomains saved to: ", PROCESSED_SUBDOMAIN_FILE)
message("Script 02 complete.")