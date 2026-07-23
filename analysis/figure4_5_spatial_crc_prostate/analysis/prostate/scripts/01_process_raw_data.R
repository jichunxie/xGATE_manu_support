# --- 1. SETUP ---
# This script processes the raw 10x-like spatial data for the prostate
# sample and converts it into the `sp_input` format required for IRIS.

# Source the shared data ingestion function
# Note the path: we go up one level to the project root to find the R/ folder
source("../../R/data_ingestion_functions.R")

# Source the dataset-specific configuration
source("00_config_prostate.R")

# --- 2. PROCESS RAW DATA ---
message("Building sp_input object for the prostate dataset...")
sp_input <- build_sp_input_slice(
  matrix_fp      = MATRIX_FILE,
  barcodes_fp    = BARCODES_FILE,
  features_fp    = FEATURES_FILE,
  spatial_loc_fp = SPATIAL_LOC_FILE,
  save_rds       = SP_INPUT_FILE
)

message("Prostate raw data processing complete.")