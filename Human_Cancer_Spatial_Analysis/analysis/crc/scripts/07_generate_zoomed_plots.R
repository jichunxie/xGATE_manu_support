# --- 1. SETUP ---
source("scripts/00_config.R")
source("R/visualization_functions.R")

# --- 2. LOAD PROCESSED DATA ---
message("Loading pre-processed data...")
processed_data <- readRDS(file.path(PROCESSED_DATA_DIR, "processed_colon_data.rds"))
spatial_df <- processed_data$spatial_df
pathway_agg <- processed_data$pathway_agg

# --- 3. GENERATE AND SAVE PLOTS ---
unique_pathways <- unique(pathway_agg$pathway)
message(sprintf("Generating %d zoomed pathway maps for section '%s'...", length(unique_pathways), ACTIVE_ZOOM_SECTION))

legend_saved <- FALSE
for (pw in unique_pathways) {
  pathway_subset <- filter(pathway_agg, pathway == pw)
  spatial_df_with_pval <- left_join(spatial_df, pathway_subset, by = c("domain_num", "cluster_num"))
  
  base_name <- paste0("pathway_", gsub("[^A-Za-z0-9]+", "_", pw), "_zoomed_", ACTIVE_ZOOM_SECTION)
  out_file <- file.path(ZOOM_PLOT_DIR, paste0(base_name, ".jpg"))
  
  # Generate plot and get the version with a legend back
  p_legend <- plot_pathway_map(
    spatial_df_with_pval = spatial_df_with_pval,
    pathway_name = pw,
    output_file = out_file,
    color_values = PVAL_BREAKS,
    color_palette = PVAL_COLORS,
    zoom_coords = ACTIVE_ZOOM, # Use the zoom coordinates from config
    point_size = 0.3,
    plot_width = 2.5,
    plot_height = 2.5,
    dpi = 800
  )
  cat("Saved:", out_file, "\n")
  
  # Save the legend just once
  if (!legend_saved && !is.null(p_legend)) {
    legend_file <- file.path(ZOOM_PLOT_DIR, "pvalue_legend.jpg")
    save_colorbar_legend(p_legend, legend_file, width = 2.8, height = 0.55, dpi = 800)
    message("Saved shared legend: ", legend_file)
    legend_saved <- TRUE
  }
}
message("Done generating zoomed plots.")