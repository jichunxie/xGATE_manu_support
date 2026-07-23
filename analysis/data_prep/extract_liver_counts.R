# =============================================================================
# Extract hepatocyte SCTransform corrected counts from liver Seurat object
# Seurat v5-compatible version: uses layer= instead of slot=
# Output CSV format:
#   first column = gene_id
#   remaining columns = cell IDs
# =============================================================================

suppressMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(sctransform)
})

set.seed(1234)

# ----------------------------
# Paths
# ----------------------------

liver_rds <- "/path/to/xGATE/data/liver/human liver hepatocyte cells.rds"
out_csv   <- "/path/to/xGATE/data/liver/liver_counts_ensembl.csv"

# ----------------------------
# Step 1: Load Seurat object
# ----------------------------

liver_data <- readRDS(liver_rds)

metadata <- liver_data@meta.data

cat("Metadata columns:\n")
print(colnames(metadata))

cat("\nCell type distribution:\n")
print(table(metadata$cell_type, useNA = "ifany"))

# ----------------------------
# Step 2: Select hepatocyte cells
# ----------------------------

hepatocyte_cells <- rownames(metadata)[metadata$cell_type == "hepatocyte"]

cat("\nNumber of hepatocyte cells in metadata:", length(hepatocyte_cells), "\n")

# Seurat v5-compatible access to RNA counts
rna_counts <- GetAssayData(
  object = liver_data,
  assay = "RNA",
  layer = "counts"
)

hepatocyte_cells <- intersect(hepatocyte_cells, colnames(rna_counts))

cat("Number of hepatocyte cells in RNA matrix:", length(hepatocyte_cells), "\n")

liver_raw_counts <- rna_counts[, hepatocyte_cells, drop = FALSE]

cat(
  "Raw hepatocyte matrix:",
  nrow(liver_raw_counts), "genes x",
  ncol(liver_raw_counts), "cells\n"
)

# ----------------------------
# Step 3: Create hepatocyte-only Seurat object
# ----------------------------

so <- CreateSeuratObject(counts = liver_raw_counts)

# ----------------------------
# Step 4: Run SCTransform
# ----------------------------

so <- SCTransform(
  object = so,
  vst.flavor = "v2",
  verbose = FALSE
)

# In Seurat v5:
# SCT layer "counts"     = corrected counts
# SCT layer "data"       = log1p corrected counts
# SCT layer "scale.data" = Pearson residuals / scaled values
#
# For the current competing-metrics benchmark, use SCT corrected counts.
liver_sct_counts <- GetAssayData(
  object = so,
  assay = "SCT",
  layer = "counts"
)

cat(
  "SCT corrected count matrix:",
  nrow(liver_sct_counts), "genes x",
  ncol(liver_sct_counts), "cells\n"
)

# ----------------------------
# Step 5: Check gene IDs
# ----------------------------

gene_ids <- rownames(liver_sct_counts)

cat("\nFirst 10 gene IDs:\n")
print(head(gene_ids, 10))

is_ensembl <- mean(grepl("^ENSG", gene_ids)) > 0.5

if (!is_ensembl) {
  warning(
    "Most rownames do not look like ENSEMBL IDs. ",
    "They may be gene symbols. If your gene-set JSON is ENSEMBL-based, ",
    "we need to map symbols to ENSEMBL or regenerate the gene sets as symbols."
  )
} else {
  cat("Gene IDs look like ENSEMBL IDs.\n")
}

# ----------------------------
# Step 6: Remove duplicated gene IDs if present
# ----------------------------

dup_genes <- duplicated(rownames(liver_sct_counts))

if (any(dup_genes)) {
  cat("Removing duplicated gene IDs:", sum(dup_genes), "\n")
  liver_sct_counts <- liver_sct_counts[!dup_genes, , drop = FALSE]
}

# ----------------------------
# Step 7: Write to CSV
# ----------------------------

# data.table::fwrite cannot directly write Matrix objects.
# Convert to data.table and add gene_id as first column.
out_dt <- as.data.table(as.matrix(liver_sct_counts))
out_dt <- cbind(gene_id = rownames(liver_sct_counts), out_dt)

fwrite(out_dt, out_csv)

cat("\nWrote:", out_csv, "\n")
cat(
  "Final exported matrix:",
  nrow(liver_sct_counts), "genes x",
  ncol(liver_sct_counts), "cells\n"
)