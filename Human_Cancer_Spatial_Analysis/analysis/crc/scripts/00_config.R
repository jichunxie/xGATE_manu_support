# A. Spatial Domain Identification (IRIS)

# -- File Paths --
# Raw input data for IRIS
RAW_SC_INPUT_FILE <- file.path(RAW_DATA_DIR, "sc_input_crc.RDS")
RAW_SP_INPUT_FILE <- file.path(RAW_DATA_DIR, "countList_spatial_crc.RDS")

# Processed output: the full IRIS object after domain identification
# Saving this avoids re-running the time-consuming analysis
PROCESSED_IRIS_OBJECT_FILE <- file.path(PROCESSED_DATA_DIR, "processed_iris_object.rds")

# Final figure output
IRIS_DOMAIN_PLOT_FILE <- file.path(RESULTS_DIR, "figures", "iris_spatial_domains.jpg")

# -- Analysis Parameters --
# Number of parallel cores to use for analysis
N_CORES <- 40

# Number of spatial domains (clusters) to identify
IRIS_N_DOMAINS <- 10

# -- Visualization & Remapping --
# Mapping to re-level the domain IDs for consistent visualization
# This maps the original IRIS domain ID (name) to the desired final ID (value)
IRIS_DOMAIN_MAP <- c(
  "0" = "5", "1" = "2", "2" = "3", "3" = "1", "4" = "0", "5" = "4", "6" = "6"
)
# Note: The 'HEATMAP_DOMAIN_COLORS' palette defined in Section E will be used for consistency.

# B. Subdomain Definition

# -- File Paths --
# Input for this step is the output from the cell type prep script,
# as it contains the most complete spatial dataframe.
SUBDOMAIN_INPUT_FILE <- PROCESSED_CELLTYPE_FILE # From Section E
# Output file with subdomains added
PROCESSED_SUBDOMAIN_FILE <- file.path(PROCESSED_DATA_DIR, "spatial_data_with_subdomains.rds")

# Final figure output
SUBDOMAIN_PLOT_FILE <- file.path(RESULTS_DIR, "figures", "spatial_subdomains_map.jpg")

# -- Analysis Parameters --
# The target minimum size for a subdomain. Domains larger than this
# will be split into balanced clusters via k-means.
SUBDOMAIN_MIN_SIZE <- 5000

# Final figure output
SUBDOMAIN_PLOT_FILE <- file.path(RESULTS_DIR, "figures", "spatial_subdomains_map.jpg")
SUBDOMAIN_ZOOM_PLOT_FILE <- file.path(RESULTS_DIR, "figures", "spatial_subdomains_map_zoomed.jpg") 

# C. File and Directory Paths
# Using relative paths makes the project portable.
# This assumes you run the scripts from the project's root directory.
RAW_DATA_DIR <- "data/raw"
PROCESSED_DATA_DIR <- "data/processed"
RESULTS_DIR <- "results"

# Input files
PATHWAY_FILE <- file.path(RAW_DATA_DIR, "colon_pathway_analysis_results_cpact_subdomain.csv")
SPATIAL_FILE <- file.path(RAW_DATA_DIR, "sp_input_with_subdomains.rds")

# Output directories
FULL_PLOT_DIR <- file.path(RESULTS_DIR, "figures", "full_view")
ZOOM_PLOT_DIR <- file.path(RESULTS_DIR, "figures", "zoomed_view")

# Create output directories if they don't exist
dir.create(FULL_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ZOOM_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# D. Plotting Parameters

# Color scale definition (defined ONCE)
PVAL_COLORS <- c("#7f0000", "darkorange3", "orange", "#cfe8ff", "#76b4e6", "#2f6fb3", "#0b2c55")
PVAL_BREAKS <- c(0, 0.02, 0.05, 0.0500001, 0.2, 0.5, 1)

# Zoom coordinates for different sections
ZOOM_COORDS <- list(
  A = list(xlim = c(53000, 56300) , ylim = c(-283, 4000)),
  B = list(xlim = c(56000, 65000), ylim = c(1300, 9000)),
  C = list(xlim = c(60000, 64500), ylim = c(9000, 14000))
)

# Select which zoom section to use for the analysis
# This makes it easy to switch between sections A, B, C
ACTIVE_ZOOM_SECTION <- "C"
ACTIVE_ZOOM <- ZOOM_COORDS[[ACTIVE_ZOOM_SECTION]]


# E. Load Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(scales)
  library(cowplot)
  library(grid)
})


# F. Correlation Analysis Parameters

# Directory for correlation results
CORR_PLOT_DIR <- file.path(RESULTS_DIR, "figures", "correlation_heatmaps")
dir.create(CORR_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Correlation settings
CORR_METHOD         <- "pearson"  # "pearson" or "spearman"
CORR_MIN_NON_NA     <- 2          # Min non-NA subdomains for a pathway to be included
CORR_PVAL_TRANSFORM <- "minuslog10" # "minuslog10" or "raw"

# A specific list of pathways to highlight and order in the heatmap
SELECTED_PATHWAYS <- c(
  "Immunosuppressive",
  "JAK-STAT signaling pathway",
  "NF-kappa B signaling pathway",
  "TLS",
  "EGFR tyrosine kinase inhibitor resistance"
)
# Match pathways case-insensitively?
CASE_INSENSITIVE_MATCH <- TRUE

# Heatmap label display settings
HEATMAP_LABEL_MAX_CHARS <- 40 # Truncate long pathway names


# G. Cell Type Analysis Parameters

# -- File Paths --
# Input file with detailed cell type annotations
CELLTYPE_INPUT_FILE <- file.path(RAW_DATA_DIR, "sp_input_with_domains_with_celltype_and_cluster.rds")
# Processed output file with general cell types mapped
PROCESSED_CELLTYPE_FILE <- file.path(PROCESSED_DATA_DIR, "spatial_data_with_general_celltypes.rds")

# Output plot paths
CELLTYPE_SPATIAL_PLOT_FILE <- file.path(RESULTS_DIR, "figures", "spatial_celltype_map.jpg")
CELLTYPE_LOGODDS_HEATMAP_FILE <- file.path(RESULTS_DIR, "figures", "celltype_logodds_heatmap.jpg")

# -- Cell Type Mapping --
# Mapping from detailed cell types to general categories
CELL_TYPE_MAP <- list(
  "T/NK cells" = c("CD4+ T cells", "CD8+ T cells", "Regulatory T cells",
                   "T follicular helper cells", "T helper 17 cells", "gamma delta T cells",
                   "NK cells"),
  "B cells" = c("CD19+CD20+ B"),
  "Plasma cells" = c("IgA+ Plasma", "IgG+ Plasma"),
  "Endothelial" = c("Lymphatic ECs", "Proliferative ECs", "Stalk-like ECs", "Tip-like ECs"),
  "Mature Enterocytes" = c("Mature Enterocytes type 1", "Mature Enterocytes type 2"),
  "Goblet cells" = c("Goblet cells"),
  "Epithelial" = c("Stem-like/TA"),
  "Fibroblasts" = c("Myofibroblasts"),
  "Glial cells" = c("Enteric glial cells"),
  "Smooth muscle" = c("Smooth muscle cells"),
  "Pericytes" = c("Pericytes"),
  "Dendritic cells" = c("cDC"),
  "Stromal cells" = c("Stromal 1", "Stromal 2", "Stromal 3"),
  "Tumor cells" = c("CMS1", "CMS2", "CMS3", "SPP1+", "Pro-inflammatory", "Intermediate", "Proliferating", "Unknown"),
  "Mast cells" = c("Mast cells")
)

# -- Plotting Palettes and Orderings --
# Custom color palette for the spatial cell type plot
CELLTYPE_COLOR_PALETTE <- c(
  "T/NK cells" = "#1f77b4", "B cells" = "#2ca02c", "Plasma cells" = "#17becf",
  "Endothelial" = "#9467bd", "Mature Enterocytes" = "#bcbd22", "Goblet cells" = "#ff7f0e",
  "Epithelial" = "#8c564b", "Fibroblasts" = "#e377c2", "Glial cells" = "#7f7f7f",
  "Smooth muscle" = "#000000", "Pericytes" = "#aec7e8", "Dendritic cells" = "#ffbb78",
  "Stromal cells" = "#b5bd61", "Tumor cells" = "#e41a1c", "Mast cells" = "#f781bf",
  "Other" = "#bbbbbb"
)

# Order for heatmap rows
HEATMAP_CELL_TYPE_ORDER <- c(
  "T/NK cells", "B cells", "Plasma cells", # Immune
  "Epithelial", "Mature Enterocytes", "Goblet cells", # Epithelial
  "Fibroblasts", "Pericytes", "Smooth muscle", "Endothelial", "Glial cells", # Stromal
  "Tumor cells", "Dendritic cells", "Mast cells",
  "Stromal cells", "Other"
)

# Order and colors for heatmap columns (domains)
HEATMAP_DOMAIN_ORDER <- c(5, 2, 3, 1, 0, 4, 6, 7, 8, 9)
HEATMAP_DOMAIN_COLORS <- c(
  "0"="#1697a6", "1"="#F7CA71", "2"="#D57358", "3"="#023047", "4"="#8bc6cc",
  "5"="#ebe5c2", "6"="#C9DEC3", "7"="#A45EE5", "8"="#FF6F61", "9"="#4CAF50"
)

# H. Subdomain Pie Chart Parameters

# Output directory for the pie charts
PIE_CHART_DIR <- file.path(RESULTS_DIR, "figures", "subdomain_pie_charts")
dir.create(PIE_CHART_DIR, recursive = TRUE, showWarnings = FALSE)

# A specific list of subdomains to generate pie charts for
PIE_CHART_SUBDOMAINS <- c(
  "domain_9_3",
  "domain_8_19",
  "domain_6_1",
  "domain_7_5",
  "domain_2_5"
)
