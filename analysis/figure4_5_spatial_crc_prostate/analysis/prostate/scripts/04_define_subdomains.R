# --- 1. SETUP ---
# This script defines spatial subdomains for the prostate dataset by
# partitioning large domains into smaller, balanced clusters.

# Source the shared spatial clustering function
source("../../R/spatial_clustering_functions.R")
# Source the dataset-specific configuration
source("00_config_prostate.R")

# --- 2. LOAD DATA WITH DOMAINS ---
message("Loading prostate data with domains: ", PROSTATE_WITH_DOMAINS_FILE)
if (!file.exists(PROSTATE_WITH_DOMAINS_FILE)) {
  stop("Input file with domains not found. Please run the domain identification script for the prostate dataset first.")
}
sp_input <- readRDS(PROSTATE_WITH_DOMAINS_FILE)

# --- 3. CLUSTER DOMAINS TO CREATE SUBDOMAINS ---
message("Defining subdomains with a target minimum size of ", SUBDOMAIN_MIN_SIZE, "...")
# Apply the shared clustering function to the prostate data
sp_input$spatial_location_list <- lapply(
  sp_input$spatial_location_list,
  function(df) cluster_domains_balanced(df, min_size = SUBDOMAIN_MIN_SIZE)
)

# --- 4. SAVE DATA WITH SUBDOMAINS ---
saveRDS(sp_input, file = PROSTATE_WITH_SUBDOMAINS_FILE)
message("Processed prostate data with subdomains saved to: ", PROSTATE_WITH_SUBDOMAINS_FILE)
message("Prostate subdomain definition complete.")