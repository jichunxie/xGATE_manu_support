# --- 1. SETUP ---
# This script performs spatial domain identification on the prostate dataset using IRIS.
# It uses a pre-annotated spatial object and a single-cell reference to define domains.

suppressPackageStartupMessages({
  library(dplyr)
  library(IRIS)
  library(ggplot2)
})

# Source the dataset-specific configuration
source("00_config_prostate.R")

# --- 2. LOAD DATA ---
message("Loading annotated spatial data and single-cell reference...")
sp_input_prostate <- readRDS(SP_INPUT_ANNOTATED_FILE)
sc_input_prostate <- readRDS(REF_SC_INPUT_FILE)

# --- 3. PREPARE SINGLE-CELL REFERENCE ---
# This mapping function is specific to the prostate cancer reference data.
map_to_general_cancer <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("T cell", "B cell", "NK cell", "MNP", "Mast cell") ~ "Immune",
    x %in% c("LE-KLK3", "LE-KLK4") ~ "Luminal Epithelial (Tumor)",
    x %in% c("BE", "HE", "CE") ~ "Basal/Secretory Epithelial",
    x %in% c("Fibroblast", "Endothelial") ~ "Stroma",
    TRUE ~ "Other"
  )
}

# Apply the mapping to create a new metadata column for IRIS
sc_input_prostate$sc_meta <- sc_input_prostate$sc_meta %>%
  mutate(general_cell_type_cancer = map_to_general_cancer(cell_type))

# --- 4. RUN IRIS DOMAIN IDENTIFICATION ---
# This section can be skipped if the processed IRIS object already exists.
if (file.exists(PROSTATE_IRIS_OBJECT_FILE)) {
  message("Loading existing processed IRIS object from: ", PROSTATE_IRIS_OBJECT_FILE)
  IRIS_object <- readRDS(PROSTATE_IRIS_OBJECT_FILE)
} else {
  message("Creating new IRIS object...")
  IRIS_object <- createIRISObject(
    spatial_countMat_list = sp_input_prostate$spatial_countMat_list,
    spatial_location_list = sp_input_prostate$spatial_location_list,
    sc_count = sc_input_prostate$sc_count,
    sc_meta = sc_input_prostate$sc_meta,
    ct.varname = 'general_cell_type_cancer',
    sample.varname = sc_input_prostate$sample.varname,
    minCountGene = 100,
    minCountSpot = 5
  )

  message("Performing domain detection for ", IRIS_N_DOMAINS, " domains...")
  IRIS_object <- IRIS_spatial(IRIS_object, numCluster = IRIS_N_DOMAINS)

  message("Saving processed IRIS object to: ", PROSTATE_IRIS_OBJECT_FILE)
  saveRDS(IRIS_object, file = PROSTATE_IRIS_OBJECT_FILE)
}

# --- 5. VISUALIZE DOMAINS ---
message("Preparing data for visualization...")
# Extract domain labels and re-level them using the map from the config file
iris_domain_df <- IRIS_object@spatialDomain[, c("Slice", "spotName", "IRIS_domain")]
iris_domain_df$IRIS_domain <- plyr::mapvalues(
  iris_domain_df$IRIS_domain,
  from = names(PROSTATE_DOMAIN_MAP),
  to = as.character(PROSTATE_DOMAIN_MAP)
)

# Extract and transform spatial coordinates for visualization
spatial_location <- data.frame(
  x = IRIS_object@spatialDomain$y,
  y = -IRIS_object@spatialDomain$x
)

message("Generating domain visualization plot...")
domain_plot <- IRIS.visualize.domain(
  iris_domain_df,
  spatial_location,
  colors = PROSTATE_DOMAIN_COLORS,
  numCols = 4
)

# --- 6. SAVE PLOT ---
ggsave(
  filename = PROSTATE_DOMAIN_PLOT_FILE,
  plot = domain_plot,
  width = 6,
  height = 12,
  dpi = 300,
  bg = "white"
)

message("Prostate domain visualization saved to: ", PROSTATE_DOMAIN_PLOT_FILE)
message("Prostate domain identification complete.")