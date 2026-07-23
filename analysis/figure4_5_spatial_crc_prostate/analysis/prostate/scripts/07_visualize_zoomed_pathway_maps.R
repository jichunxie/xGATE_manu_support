# --- 1. SETUP ---
# This script generates zoomed-in spatial maps for all pathway activities
# in a specified region of interest (ROI).

source("../../R/visualization_functions.R")
source("00_config_prostate.R")

# --- 2. LOAD PREPROCESSED DATA ---
message("Loading preprocessed spatial-pathway data: ", PROSTATE_SPATIAL_PATHWAY_FILE)
processed_data <- readRDS(PROSTATE_SPATIAL_PATHWAY_FILE)
spatial_df <- processed_data$spatial_df
pathway_res <- processed_data$pathway_res

# --- 3. GENERATE ZOOMED PLOTS FOR EACH PATHWAY ---
# Select the active zoom coordinates from the config file
active_zoom_coords <- PROSTATE_ZOOM_REGIONS[[ACTIVE_PATHWAY_ZOOM_REGION]]
message("Generating zoomed pathway maps for ROI: ", ACTIVE_PATHWAY_ZOOM_REGION)

unique_pathways <- unique(pathway_res$pathway)
for (pw in unique_pathways) {
  # Aggregate pathway results to get one p-value per subdomain
  pathway_subset <- pathway_res %>%
    filter(pathway == pw) %>%
    group_by(domain_num, cluster_num) %>%
    summarize(p.value = min(p.value, na.rm = TRUE), .groups = "drop")
  
  # Join p-values to the spatial data frame
  spatial_df_with_pval <- left_join(spatial_df, pathway_subset, by = c("domain_num", "cluster_num"))
  
  # Define a clean output filename
  safe_pw_name <- gsub("[^A-Za-z0-9]+", "_", pw)
  output_file <- file.path(
    PROSTATE_ZOOMED_PATHWAY_DIR,
    paste0("pathway_", safe_pw_name, "_zoom_", ACTIVE_PATHWAY_ZOOM_REGION, ".jpg")
  )
  
  # Call our shared plotting function
  plot_pathway_map(
    spatial_df_with_pval = spatial_df_with_pval,
    pathway_name = pw,
    output_file = output_file,
    color_values = PVAL_COLOR_VALS,
    color_palette = PVAL_COLORS,
    zoom_coords = active_zoom_coords,
    point_size = 0.2,
    width = 3, height = 3, dpi = 800
  )
}

message("Finished generating zoomed pathway plots.")