# --- 1. SETUP ---
source("scripts/00_config.R")
source("R/analysis_functions.R")

# --- 2. LOAD PROCESSED DATA ---
message("Loading processed data with general cell types: ", PROCESSED_CELLTYPE_FILE)
sp_input <- readRDS(PROCESSED_CELLTYPE_FILE)
meta_df <- sp_input$spatial_location_list[[1]] %>%
  select(domain, general_cell_type) %>%
  filter(!is.na(domain), !is.na(general_cell_type))

# --- 3. CALCULATE LOG ODDS RATIOS ---
message("Calculating log odds ratios for cell type enrichment in domains...")
# Create a grid of all domain/cell type combinations to test
all_combos <- expand.grid(
  domain = unique(meta_df$domain),
  general_cell_type = unique(meta_df$general_cell_type),
  stringsAsFactors = FALSE
)

# Calculate log odds for each combination
log_odds_df <- all_combos %>%
  rowwise() %>%
  mutate(log2_odds = compute_log_odds(domain, general_cell_type, meta_df)) %>%
  ungroup()

# --- 4. PREPARE MATRIX FOR HEATMAP ---
message("Preparing matrix for heatmap...")
# Pivot to wide format
log_odds_matrix <- log_odds_df %>%
  pivot_wider(names_from = domain, values_from = log2_odds)

# Convert to a numeric matrix with cell types as rownames
log_odds_mat <- as.matrix(log_odds_matrix[,-1])
rownames(log_odds_mat) <- log_odds_matrix$general_cell_type

# Filter out non-positive values
log_odds_mat[log_odds_mat <= 0] <- NA

# Reorder rows and columns based on config
ordered_rows <- intersect(HEATMAP_CELL_TYPE_ORDER, rownames(log_odds_mat))
ordered_cols <- paste0("domain_", intersect(HEATMAP_DOMAIN_ORDER, gsub("domain_", "", colnames(log_odds_mat))))
log_odds_mat <- log_odds_mat[ordered_rows, ordered_cols]

# Remove rows/columns that are now all NA
log_odds_mat <- log_odds_mat[rowSums(is.finite(log_odds_mat)) > 0, ]
log_odds_mat <- log_odds_mat[, colSums(is.finite(log_odds_mat)) > 0]

# --- 5. GENERATE HEATMAP ---
message("Generating log odds heatmap...")
col_fun <- colorRamp2(seq(0, max(log_odds_mat, na.rm = TRUE), length.out = 256), viridis::magma(256))

domain_ids <- gsub("domain_", "", colnames(log_odds_mat))
ha <- HeatmapAnnotation(
  Domain = anno_simple(domain_ids, col = HEATMAP_DOMAIN_COLORS[domain_ids]),
  annotation_name_side = "left"
)

ht <- Heatmap(
  log_odds_mat,
  name = "Log2 Odds Ratio",
  col = col_fun,
  na_col = "grey95",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  bottom_annotation = ha,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10)
)

# --- 6. SAVE HEATMAP ---
jpeg(filename = CELLTYPE_LOGODDS_HEATMAP_FILE, width = 8, height = 7, units = "in", res = 400)
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
message("Log odds heatmap saved to: ", CELLTYPE_LOGODDS_HEATMAP_FILE)