# --- 1. SETUP ---
# This script visualizes the newly defined spatial subdomains for the
# prostate dataset (full view).

# Source the shared visualization functions
source("../../R/visualization_functions.R")
# Source the dataset-specific configuration
source("00_config_prostate.R")

# --- 2. LOAD DATA WITH SUBDOMAINS ---
message("Loading prostate data with defined subdomains: ", PROSTATE_WITH_SUBDOMAINS_FILE)
if (!file.exists(PROSTATE_WITH_SUBDOMAINS_FILE)) {
  stop("Subdomain file not found. Please run '02_define_subdomains.R' for the prostate dataset first.")
}
sp_input <- readRDS(PROSTATE_WITH_SUBDOMAINS_FILE)
plot_df <- sp_input$spatial_location_list[[1]] # Visualizing the first slice

# --- 3. GENERATE AND SAVE PLOT ---
message("Generating full-view subdomain map for the prostate dataset...")
plot_subdomain_map(
  spatial_df = plot_df,
  output_file = SUBDOMAIN_PLOT_FILE,
  zoom_coords = NULL, # Full view
  # Adjust dimensions for the prostate sample's aspect ratio
  plot_width = 12,
  plot_height = 16,
  dpi = 600
)

message("Prostate subdomain visualization complete.")