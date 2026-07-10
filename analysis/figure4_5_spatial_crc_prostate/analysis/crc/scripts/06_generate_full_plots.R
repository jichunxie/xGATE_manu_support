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
message(sprintf("Generating %d full-view pathway maps...", length(unique_pathways)))

# Use a loop to generate a plot for each pathway
for (pw in unique_pathways) {
  # Filter data for the current pathway
  pathway_subset <- filter(pathway_agg, pathway == pw)
  
  # Join p-values to spatial data
  spatial_df_with_pval <- left_join(spatial_df, pathway_subset, by = c("domain_num", "cluster_num"))
  
  # Define output file path
  out_file <- file.path(FULL_PLOT_DIR, paste0("pathway_", gsub("[^A-Za-z0-9]+", "_", pw), "_full.jpg"))
  
  # Generate and save plot
  plot_pathway_map(
    spatial_df_with_pval = spatial_df_with_pval,
    pathway_name = pw,
    output_file = out_file,
    color_values = PVAL_BREAKS,
    color_palette = PVAL_COLORS,
    point_size = 0.22,
    plot_width = 6,
    plot_height = 5,
    dpi = 800
  )
  cat("Saved:", out_file, "\n")
}
message("Done generating full plots.")