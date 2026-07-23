# --- 1. SETUP ---
source("scripts/00_config.R")
# No custom functions needed, just ggplot2 loaded from config

# --- 2. LOAD PROCESSED DATA ---
message("Loading processed data with general cell types: ", PROCESSED_CELLTYPE_FILE)
sp_input <- readRDS(PROCESSED_CELLTYPE_FILE)
plot_df <- sp_input$spatial_location_list[[1]] %>%
  filter(!is.na(general_cell_type))

# --- 3. GENERATE SPATIAL PLOT ---
message("Generating spatial map of general cell types...")
p <- ggplot(plot_df, aes(x = x, y = y, color = general_cell_type)) +
  geom_point(size = 0.4, alpha = 0.8, shape = 16) +
  scale_color_manual(values = CELLTYPE_COLOR_PALETTE, name = "Cell type") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Spatial Distribution of General Cell Types") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) # Enlarge legend dots

# --- 4. SAVE PLOT ---
ggsave(CELLTYPE_SPATIAL_PLOT_FILE, p, width = 12, height = 9, dpi = 400, bg = "white")
message("Spatial cell type map saved to: ", CELLTYPE_SPATIAL_PLOT_FILE)