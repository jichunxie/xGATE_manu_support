# --- 1. SETUP ---
# This script performs spatial domain identification using the IRIS package.
# It loads raw single-cell and spatial data, runs the IRIS algorithm,
# saves the resulting object, and generates a visualization of the domains.
source("scripts/00_config.R")

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(IRIS)
  library(doMC)
  library(ggplot2)
})

# Register parallel backend
registerDoMC(cores = N_CORES)
message("Registered ", N_CORES, " cores for parallel processing.")

# --- 2. LOAD RAW DATA ---
message("Loading raw single-cell and spatial data...")
sc_input <- readRDS(RAW_SC_INPUT_FILE)
sp_input <- readRDS(RAW_SP_INPUT_FILE)

# --- 3. RUN IRIS DOMAIN IDENTIFICATION ---
# This section can be skipped if the processed file already exists.
if (file.exists(PROCESSED_IRIS_OBJECT_FILE)) {
  message("Loading existing processed IRIS object from: ", PROCESSED_IRIS_OBJECT_FILE)
  IRIS_object <- readRDS(PROCESSED_IRIS_OBJECT_FILE)
} else {
  message("Creating new IRIS object...")
  IRIS_object <- createIRISObject(
    spatial_countMat_list = sp_input$spatial_countMat_list,
    spatial_location_list = sp_input$spatial_location_list,
    sc_count = sc_input$sc_count,
    sc_meta = sc_input$sc_meta,
    ct.varname = sc_input$ct.varname,
    sample.varname = sc_input$sample.varname,
    minCountGene = 100,
    minCountSpot = 5
  )

  message("Performing integrative reference-informed domain detection for ", IRIS_N_DOMAINS, " domains...")
  IRIS_object <- IRIS_spatial(IRIS_object, numCluster = IRIS_N_DOMAINS)

  message("Saving processed IRIS object to: ", PROCESSED_IRIS_OBJECT_FILE)
  saveRDS(IRIS_object, file = PROCESSED_IRIS_OBJECT_FILE)
}

# --- 4. POST-PROCESS AND VISUALIZE DOMAINS ---
message("Preparing data for visualization...")
# Extract domain labels detected by IRIS
iris_domain_df <- IRIS_object@spatialDomain[, c("Slice", "spotName", "IRIS_domain")]

# Re-level the spatial domains using the map from the config file
iris_domain_df$IRIS_domain <- plyr::mapvalues(
  iris_domain_df$IRIS_domain,
  from = names(IRIS_DOMAIN_MAP),
  to = as.character(IRIS_DOMAIN_MAP)
)

# Extract and transform spatial locations for visualization consistency
# Note: This transformation is for visualization purposes only.
spatial_location <- data.frame(
  x = IRIS_object@spatialDomain$y,
  y = -IRIS_object@spatialDomain$x
)

message("Generating domain visualization plot...")
# Use the consistent domain color palette from the config
domain_plot <- IRIS.visualize.domain(
  iris_domain_df,
  spatial_location,
  colors = HEATMAP_DOMAIN_COLORS,
  numCols = 4
)

# --- 5. SAVE PLOT ---
ggsave(
  filename = IRIS_DOMAIN_PLOT_FILE,
  plot = domain_plot,
  width = 8,
  height = 8,
  dpi = 600,
  bg = "white"
)
message("Domain visualization saved to: ", IRIS_DOMAIN_PLOT_FILE)
message("Script 01 complete.")