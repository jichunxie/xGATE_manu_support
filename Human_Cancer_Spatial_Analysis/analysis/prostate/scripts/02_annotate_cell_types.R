# --- 1. SETUP ---
# This script performs cell type annotation on the prostate spatial data
# using a single-cell reference dataset and the SingleR package.

# Source the shared cell annotation functions from the main R/ directory
source("../../R/cell_annotation_functions.R")
# Source the dataset-specific configuration for prostate
source("00_config_prostate.R")

# --- 2. RUN ANNOTATION ---
message("Starting cell type annotation for the prostate dataset...")

# Check if the input file from the previous step exists
if (!file.exists(SP_INPUT_FILE)) {
  stop("Input file '", SP_INPUT_FILE, "' not found. Please run '01_process_raw_data.R' first.")
}

# Call the main annotation function with parameters from the config file
annotate_spatial_with_singleR(
  spatial_rds = SP_INPUT_FILE,
  out_rds = SP_INPUT_ANNOTATED_FILE,
  
  # Using an SCE reference as defined in the config file
  ref_sce_rds = REF_SCE_FILE,
  ref_sce_label_col = REF_SCE_LABEL_COL,
  
  # You could also specify a Seurat reference here instead, for example:
  # ref_seurat_rds = REF_SEURAT_FILE,
  # ref_seurat_label_col = REF_SEURAT_LABEL_COL,
  
  # Pass other parameters from the config
  workers = ANNOTATION_WORKERS,
  aggregate_reference = ANNOTATION_AGGREGATE_REF,
  min_cells_per_label = 25
)

message("Cell type annotation for the prostate dataset is complete.")