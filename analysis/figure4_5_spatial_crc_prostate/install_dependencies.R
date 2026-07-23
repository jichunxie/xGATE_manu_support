# install_dependencies.R 

# --- CRAN Packages ---
cran_packages <- c(
  "ggplot2", "readr", "dplyr", "stringr", "scales", "cowplot", "grid",
  "tidyr", "ComplexHeatmap", "circlize", "RColorBrewer", "Polychrome",
  "viridis", "Seurat", "Matrix", "IRIS", "doMC", "plyr", "dbscan", "patchwork"
)

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

# --- Bioconductor Packages (New Section) ---
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
bioc_packages <- c(
  "SingleR",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scuttle",
  "BiocParallel"
)

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

message("All required packages are installed.")