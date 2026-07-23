# --- 1. SETUP ---
# This script generates individual zoomed-in plots for each region of
# interest (ROI) defined in the configuration file.

suppressPackageStartupMessages({
  library(glue)
})

# Source shared functions and dataset-specific config
source("../../R/visualization_functions.R")
source("00_config_prostate.R")

# --- 2. LOAD DATA & PREPARE FOR PLOTTING ---
message("Loading prostate data with subdomains: ", PROSTATE_WITH_SUBDOMAINS_FILE)
sp_input <- readRDS(PROSTATE_WITH_SUBDOMAINS_FILE)
plot_df <- sp_input$spatial_location_list[[1]]

# Create shared color map and legend plot (code is identical to script 03)
plot_df <- plot_df %>% mutate(domain_subdomain_id = paste(domain, subdomain, sep = "_"))
color_map <- # ... (create color_map as in script 03)
legend_plot <- # ... (create legend_plot as in script 03)

# --- 3. GENERATE ZOOM PLOTS FOR EACH ROI ---
message("Generating zoomed plots for ", length(PROSTATE_ZOOM_REGIONS), " ROI(s)...")

for (region_name in names(PROSTATE_ZOOM_REGIONS)) {
  region_coords <- PROSTATE_ZOOM_REGIONS[[region_name]]
  xlim <- region_coords$xlim
  ylim <- region_coords$ylim
  
  # Filter data to the specific zoom rectangle
  df_zoom <- plot_df %>%
    filter(x >= min(xlim), x <= max(xlim), y >= min(ylim), y <= max(ylim))
  
  # Create the zoomed plot using our shared helper function
  zoom_plot <- make_zoom_plot(df_zoom, color_map, xlim, ylim, label = region_name)
  
  # Combine with the universal legend
  combined_zoom_plot <- zoom_plot + legend_plot + plot_layout(widths = c(3, 1))
  
  # Create a descriptive output filename
  output_file <- file.path(RESULTS_DIR, glue("prostate_subdomains_zoom_{region_name}.jpg"))
  
  ggsave(output_file, combined_zoom_plot, width = 12, height = 12, dpi = 300, bg = "white")
  message("Saved zoom plot: ", output_file)
}

message("Zoom plot generation complete.")