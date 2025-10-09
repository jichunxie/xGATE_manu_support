# --- 1. SETUP ---
# This script generates a zoomed-in visualization of the spatial subdomains
# for a pre-defined region of interest (ROI).
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
# Select the active zoom coordinates from the config file
active_zoom_coords <- ZOOM_COORDS[[ACTIVE_ZOOM_REGION]]
message("Generating zoomed subdomain map for ROI: ", ACTIVE_ZOOM_REGION)

plot_subdomain_map(
  spatial_df = plot_df,
  output_file = SUBDOMAIN_ZOOM_PLOT_FILE,
  zoom_coords = active_zoom_coords, # Pass the zoom coordinates
  plot_width = 10, # Adjust dimensions for a zoomed plot
  plot_height = 10
)

message("Script 04 complete.")