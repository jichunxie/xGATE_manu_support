# --- 1. SETUP ---
# This script loads and preprocesses the spatial and pathway data for the
# prostate sample, creating a unified data object for visualization.

source("../../R/data_processing_functions.R")
source("00_config_prostate.R")

# --- 2. RUN PREPROCESSING ---
message("Preprocessing prostate spatial and pathway data...")
processed_data <- preprocess_spatial_pathway_data(
  sp_rds_path = PROSTATE_WITH_SUBDOMAINS_FILE,
  pathway_csv_path = PROSTATE_PATHWAY_RESULTS_FILE
)

# --- 3. SAVE PROCESSED DATA ---
saveRDS(processed_data, file = PROSTATE_SPATIAL_PATHWAY_FILE)
message("Saved combined spatial-pathway data to: ", PROSTATE_SPATIAL_PATHWAY_FILE)