# --- 1. SETUP ---
source("scripts/00_config.R")
source("R/utils.R")

# --- 2. LOAD RAW DATA ---
message("Loading raw spatial data with cell types: ", CELLTYPE_INPUT_FILE)
sp_input <- readRDS(CELLTYPE_INPUT_FILE)
spatial_df <- sp_input$spatial_location_list[[1]]

# --- 3. MAP TO GENERAL CELL TYPES ---
message("Mapping detailed cell types to general categories...")
spatial_df$general_cell_type <- map_to_general(
  detailed_cell_types = spatial_df$cell_type,
  cell_type_map = CELL_TYPE_MAP
)

# --- 4. SAVE PROCESSED DATA ---
# Replace the original data frame in the list structure before saving
sp_input$spatial_location_list[[1]] <- spatial_df
saveRDS(sp_input, file = PROCESSED_CELLTYPE_FILE)
message("Processed data with general cell types saved to: ", PROCESSED_CELLTYPE_FILE)