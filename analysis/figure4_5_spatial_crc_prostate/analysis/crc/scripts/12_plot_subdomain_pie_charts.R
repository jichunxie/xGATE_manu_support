# --- 1. SETUP ---
# This script generates pie charts of cell type composition for specific subdomains.
source("scripts/00_config.R")
source("R/utils.R")

# --- 2. LOAD PROCESSED DATA ---
message("Loading processed data with general cell types: ", PROCESSED_CELLTYPE_FILE)
# Check if the processed file exists from the previous step
if (!file.exists(PROCESSED_CELLTYPE_FILE)) {
  stop("Processed cell type file not found. Please run '05_prepare_celltype_data.R' first.")
}
sp_input <- readRDS(PROCESSED_CELLTYPE_FILE)
spatial_df <- sp_input$spatial_location_list[[1]]

# --- 3. GENERATE PIE CHARTS FOR SELECTED SUBDOMAINS ---
# Find which of the desired subdomains are actually present in the data
present_subdomains <- unique(spatial_df$subdomain)
subdomains_to_plot <- intersect(PIE_CHART_SUBDOMAINS, present_subdomains)

if (length(subdomains_to_plot) == 0) {
  warning("None of the specified subdomains for pie charts were found in the data.")
} else {
  message("Generating pie charts for subdomains: ", paste(subdomains_to_plot, collapse = ", "))
}

for (sd in subdomains_to_plot) {
  # Filter data for the current subdomain and calculate cell type counts
  df_sub <- spatial_df %>%
    filter(subdomain == sd, !is.na(general_cell_type)) %>%
    group_by(general_cell_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n)) # Order by count descending for consistent coloring and labeling

  if (nrow(df_sub) == 0) {
    message("Skipping subdomain '", sd, "' as it contains no cells after filtering.")
    next
  }

  # Prepare elements for the pie chart
  counts <- df_sub$n
  percentages <- round(100 * counts / sum(counts), 1)

  # Create labels for only the top 3 slices
  pie_labels <- rep("", length(counts))
  top_n_labels <- min(3, length(counts))
  pie_labels[seq_len(top_n_labels)] <- paste0(percentages[seq_len(top_n_labels)], "%")

  # Get colors from our central palette
  pie_colors <- CELLTYPE_COLOR_PALETTE[as.character(df_sub$general_cell_type)]
  pie_colors[is.na(pie_colors)] <- "#808080" # Fallback color for any unmapped type

  # Define output file path
  output_file <- file.path(
    PIE_CHART_DIR,
    paste0("pie_chart_", sanitize(sd), ".jpg")
  )

  # Save the plot using base R plotting functions
  jpeg(output_file, width = 600, height = 600, quality = 100)
  par(mar = c(0, 0, 0, 0)) # Set margins to zero to maximize plot area
  pie(
    counts,
    labels = pie_labels,
    col = pie_colors,
    border = "white",
    clockwise = TRUE,
    init.angle = 90,
    cex = 2.0 # Label text size
  )
  dev.off()
}

if (length(subdomains_to_plot) > 0) {
  message("Finished. Pie charts saved to: ", PIE_CHART_DIR)
}