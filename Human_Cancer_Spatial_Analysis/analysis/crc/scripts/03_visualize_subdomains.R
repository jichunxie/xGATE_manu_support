# --- 1. SETUP ---
# This script visualizes the spatial subdomains (full view).
source("scripts/00_config.R")
source("R/visualization_functions.R") # Source our new function

# --- 2. LOAD DATA WITH SUBDOMAINS ---
message("Loading data with defined subdomains: ", PROCESSED_SUBDOMAIN_FILE)
if (!file.exists(PROCESSED_SUBDOMAIN_FILE)) {
  stop("Subdomain file not found. Please run '02_define_subdomains.R' first.")
}
sp_input <- readRDS(PROCESSED_SUBDOMAIN_FILE)
plot_df <- sp_input$spatial_location_list[[1]]

# --- 3. GENERATE AND SAVE PLOT ---
message("Generating full-view subdomain map...")
plot_subdomain_map(
  spatial_df = plot_df,
  output_file = SUBDOMAIN_PLOT_FILE,
  zoom_coords = NULL # Pass NULL to generate the full view
)

message("Script 03 complete.")