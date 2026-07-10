# --- 1. SETUP ---
# This script visualizes the pre-calculated cell type proportions (stacked bar)
# and enrichment (heatmap) for the prostate dataset.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
})

# Source the dataset-specific configuration
source("00_config_prostate.R")

# --- 2. VISUALIZE PROPORTIONS (STACKED BAR) ---
message("Generating cell type proportions bar plot...")
prop_df <- readr::read_csv(PROSTATE_CELLTYPE_PROPS_TABLE, show_col_types = FALSE)

# Define a consistent color palette for cell types
ct_levels <- sort(unique(prop_df$general_cell_type))
ct_palette <- setNames(viridis::turbo(length(ct_levels)), ct_levels)

# Order domains for the plot
domain_order <- sort(unique(prop_df$domain_num))
prop_df$domain_num <- factor(prop_df$domain_num, levels = domain_order)

p_bar <- ggplot(prop_df, aes(x = domain_num, y = prop, fill = general_cell_type)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = ct_palette) +
  labs(x = "Domain", y = "Proportion", fill = "Cell Type", title = "Cell Type Proportions per Prostate Domain") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank(), legend.position = "right")

ggsave(PROSTATE_CELLTYPE_BARPLOT_FILE, p_bar, width = 10, height = 7, dpi = 300, bg = "white")
message("Bar plot saved to: ", PROSTATE_CELLTYPE_BARPLOT_FILE)

# --- 3. VISUALIZE ENRICHMENT (HEATMAP) ---
message("Generating log odds enrichment heatmap...")
log_odds_df <- readr::read_csv(PROSTATE_CELLTYPE_LOGODDS_TABLE, show_col_types = FALSE)

# Pivot to matrix format for the heatmap
log_odds_mat <- log_odds_df %>%
  tidyr::pivot_wider(names_from = domain_num, values_from = log2_odds) %>%
  as.data.frame()
rownames(log_odds_mat) <- log_odds_mat$general_cell_type
log_odds_mat <- as.matrix(log_odds_mat[, -1])

# Filter for positive enrichment and order columns
log_odds_mat[log_odds_mat <= 0] <- NA
log_odds_mat <- log_odds_mat[, as.character(domain_order), drop = FALSE]

# Create domain annotation bar
ha <- HeatmapAnnotation(
  Domain = anno_simple(colnames(log_odds_mat), col = PROSTATE_DOMAIN_COLORS[colnames(log_odds_mat)]),
  annotation_name_side = "left"
)

# Color function for heatmap body
col_fun <- colorRamp2(seq(0, max(log_odds_mat, na.rm = TRUE), 0.01), viridis::magma(256))

ht <- Heatmap(
  log_odds_mat,
  name = "Log2 Odds Ratio",
  col = col_fun,
  na_col = "grey95",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  bottom_annotation = ha,
  row_names_gp = gpar(fontsize = 12)
)

jpeg(PROSTATE_CELLTYPE_HEATMAP_FILE, width = 8, height = 8, units = "in", res = 400)
draw(ht)
dev.off()
message("Heatmap saved to: ", PROSTATE_CELLTYPE_HEATMAP_FILE)
message("Cell type enrichment visualization complete.")