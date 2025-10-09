# Configuration file for the Prostate Cancer spatial analysis

# Define base directory for this specific analysis
PROSTATE_BASE_DIR <- "analysis/prostate"

# --- Input File Paths ---
# Raw data files are located in the 'data/raw' subdirectory
RAW_DATA_DIR <- file.path(PROSTATE_BASE_DIR, "data/raw")
MATRIX_FILE <- file.path(RAW_DATA_DIR, "matrix.mtx")
BARCODES_FILE <- file.path(RAW_DATA_DIR, "barcodes.tsv")
FEATURES_FILE <- file.path(RAW_DATA_DIR, "features.tsv")
SPATIAL_LOC_FILE <- file.path(RAW_DATA_DIR, "spatial_location.csv")

# --- Processed Data Paths ---
PROCESSED_DATA_DIR <- file.path(PROSTATE_BASE_DIR, "data/processed")
SP_INPUT_FILE <- file.path(PROCESSED_DATA_DIR, "sp_input_prostate.rds")
# Note: For this script, we need an input that already has domains identified.
# We will assume a file `sp_input_prostate_with_domains.rds` is created by a domain script.
PROSTATE_WITH_DOMAINS_FILE <- file.path(PROCESSED_DATA_DIR, "sp_input_prostate_with_domains.rds") 
PROSTATE_WITH_SUBDOMAINS_FILE <- file.path(PROCESSED_DATA_DIR, "sp_input_prostate_with_subdomains.rds") 
PROSTATE_IRIS_OBJECT_FILE <- file.path(PROCESSED_DATA_DIR, "iris_object_prostate.rds") 
PROSTATE_SPATIAL_PATHWAY_FILE <- file.path(PROCESSED_DATA_DIR, "prostate_spatial_pathway_data.rds") 
SP_INPUT_ANNOTATED_FILE <- file.path(PROCESSED_DATA_DIR, "sp_input_prostate_annotated.rds")

# --- Results Paths ---
RESULTS_DIR <- file.path(PROSTATE_BASE_DIR, "results")
SUBDOMAIN_PLOT_FILE <- file.path(RESULTS_DIR, "prostate_subdomains_map.jpg")
SUBDOMAIN_OVERVIEW_PLOT_FILE <- file.path(RESULTS_DIR, "prostate_subdomains_map_with_rois.jpg")
PROSTATE_DOMAIN_PLOT_FILE <- file.path(RESULTS_DIR, "prostate_iris_domains.jpg")
PROSTATE_CELLTYPE_BARPLOT_FILE <- file.path(RESULTS_DIR, "prostate_celltype_proportions_bar.jpg") 
PROSTATE_CELLTYPE_HEATMAP_FILE <- file.path(RESULTS_DIR, "prostate_celltype_enrichment_heatmap.jpg")
PROSTATE_ZOOMED_PATHWAY_DIR <- file.path(RESULTS_DIR, "zoomed_pathway_plots") 
dir.create(PROSTATE_ZOOMED_PATHWAY_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Pathway Analysis Inputs ---
PROSTATE_PATHWAY_RESULTS_FILE <- file.path(RAW_DATA_DIR, "prostate_pathway_analysis_results.csv") 

# --- Analysis Parameters ---
# Subdomain definition
SUBDOMAIN_MIN_SIZE <- 3000
# Parameters for cell annotation
ANNOTATION_WORKERS <- 8
ANNOTATION_AGGREGATE_REF <- TRUE

# Parameters for IRIS domain identification
IRIS_N_DOMAINS <- 10
# Mapping to re-level the domain IDs for visualization
PROSTATE_DOMAIN_MAP <- c(
  "0" = "5", "1" = "2", "2" = "3", "3" = "1", "4" = "0", "5" = "4", "6" = "6"
)
# Color palette for visualizing the domains
PROSTATE_DOMAIN_COLORS <- c("#ebe5c2", "#D57358", "#023047", "#F7CA71", "#1697a6", 
                          "#8bc6cc", "#C9DEC3", "#A45EE5", "#FF6F61", "#4CAF50")

# --- Reference Data Paths ---
# Location of the single-cell reference data
REF_DATA_DIR <- file.path(PROSTATE_BASE_DIR, "data", "reference")
# Using an SCE reference for this example
REF_SCE_FILE <- file.path(REF_DATA_DIR, "prostate_portal_300921.RDS")
REF_SCE_LABEL_COL <- "celltype"
REF_SC_INPUT_FILE <- file.path(REF_DATA_DIR, "sc_input_ProstatePortal.rds")

                          
# --- Intermediate Data Tables ---
# It's good practice to save calculated tables for easy access
TABLES_DIR <- file.path(RESULTS_DIR, "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
PROSTATE_CELLTYPE_PROPS_TABLE <- file.path(TABLES_DIR, "prostate_celltype_proportions.csv") 
PROSTATE_CELLTYPE_LOGODDS_TABLE <- file.path(TABLES_DIR, "prostate_celltype_logodds.csv") 

# --- Zoom Visualization Parameters ---
# Set to TRUE to generate the overview plot with zoom boxes
CREATE_OVERVIEW_WITH_BOXES <- FALSE

# --- Visualization Palettes & Parameters ---
# Two-phase color scale for p-values
PVAL_COLOR_VALS <- c(0, 0.02, 0.05, 0.0500001, 0.2, 0.5, 1)
PVAL_COLORS <- c("#7f0000", "darkorange3", "orange", "#cfe8ff", "#76b4e6", "#2f6fb3", "#0b2c55")


# Define multiple zoom regions of interest (ROIs) for the prostate sample.
# The script `04_visualize_subdomains_zoomed.R` will generate a plot for each.
PROSTATE_ZOOM_REGIONS <- list(
  Section_A = list(xlim = c( 8500, 13000), ylim = c(   0,  4000))
)

# Set the active ROI for the zoomed pathway plots script
ACTIVE_PATHWAY_ZOOM_REGION <- "Section_A"