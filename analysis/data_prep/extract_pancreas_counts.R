# =============================================================================
# Extract pancreatic beta-cell SCTransform corrected counts
# from pancreas_human.rds -> CSV for run_all_competing_metrics.R
#
# Input:
#   /path/to/xGATE/data/pancreas/pancreas_human.rds
#
# Output:
#   /path/to/xGATE/data/pancreas/pancreas_counts_ensembl.csv
#
# Format:
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

pancreas_rds <- "/path/to/xGATE/data/pancreas/pancreas_human.rds"
out_csv      <- "/path/to/xGATE/data/pancreas/pancreas_counts_ensembl.csv"

# ----------------------------
# Step 1: Load Seurat object
# ----------------------------

pancreas_data <- readRDS(pancreas_rds)

metadata <- pancreas_data@meta.data

cat("Metadata columns:\n")
print(colnames(metadata))

cat("\nMetadata preview:\n")
print(head(metadata))

# ----------------------------
# Step 2: Identify cell-type column
# ----------------------------

candidate_celltype_cols <- c(
  "cell_type",
  "celltype",
  "cell.types",
  "celltype_major",
  "CellType",
  "cell_type1",
  "cell_type2",
  "annotation",
  "Annotation",
  "labels",
  "seurat_annotations"
)

celltype_col <- candidate_celltype_cols[candidate_celltype_cols %in% colnames(metadata)][1]

if (is.na(celltype_col)) {
  stop(
    "Could not find a cell-type annotation column. ",
    "Please inspect colnames(metadata) and set celltype_col manually."
  )
}

cat("\nUsing cell-type column:", celltype_col, "\n")

cat("\nCell-type distribution:\n")
print(table(metadata[[celltype_col]], useNA = "ifany"))

# ----------------------------
# Step 3: Filter pancreatic beta cells
# ----------------------------

celltype_values <- as.character(metadata[[celltype_col]])

beta_pattern <- "beta|β|Beta|BETA"

beta_cells <- rownames(metadata)[grepl(beta_pattern, celltype_values)]

cat("\nNumber of beta cells detected:", length(beta_cells), "\n")

if (length(beta_cells) == 0) {
  stop(
    "No beta cells found using pattern: ", beta_pattern, "\n",
    "Please check table(metadata[[celltype_col]]) and update beta_pattern."
  )
}

# ----------------------------
# Step 4: Extract raw RNA counts
# ----------------------------

# Seurat v5-compatible assay access
if (!"RNA" %in% names(pancreas_data@assays)) {
  stop("RNA assay not found. Available assays: ", paste(names(pancreas_data@assays), collapse = ", "))
}

cat("\nRNA assay layers:\n")
print(Layers(pancreas_data[["RNA"]]))

# If there are multiple count layers, JoinLayers may be needed.
# This is safe for standard Seurat v5 objects with split layers.
if (length(grep("^counts", Layers(pancreas_data[["RNA"]]), value = TRUE)) > 1) {
  cat("\nMultiple RNA count layers detected. Joining RNA layers...\n")
  pancreas_data[["RNA"]] <- JoinLayers(pancreas_data[["RNA"]])
}

rna_counts <- GetAssayData(
  object = pancreas_data,
  assay = "RNA",
  layer = "counts"
)

beta_cells <- intersect(beta_cells, colnames(rna_counts))

cat("Number of beta cells in RNA count matrix:", length(beta_cells), "\n")

pancreas_raw_counts <- rna_counts[, beta_cells, drop = FALSE]

cat(
  "Raw beta-cell count matrix:",
  nrow(pancreas_raw_counts), "genes x",
  ncol(pancreas_raw_counts), "cells\n"
)

# ----------------------------
# Step 5: Create beta-cell-only Seurat object
# ----------------------------

so <- CreateSeuratObject(counts = pancreas_raw_counts)

# ----------------------------
# Step 6: Run SCTransform
# ----------------------------

so <- SCTransform(
  object = so,
  vst.flavor = "v2",
  verbose = FALSE
)

cat("\nSCT assay layers:\n")
print(Layers(so[["SCT"]]))

# SCT layer "counts" = corrected counts
pancreas_sct_counts <- GetAssayData(
  object = so,
  assay = "SCT",
  layer = "counts"
)

cat(
  "SCT corrected count matrix:",
  nrow(pancreas_sct_counts), "genes x",
  ncol(pancreas_sct_counts), "cells\n"
)

# ----------------------------
# Step 7: Check gene IDs
# ----------------------------

gene_ids <- rownames(pancreas_sct_counts)

cat("\nFirst 10 gene IDs:\n")
print(head(gene_ids, 10))

is_ensembl <- mean(grepl("^ENSG", gene_ids)) > 0.5

if (!is_ensembl) {
  warning(
    "Most rownames do not look like ENSEMBL IDs. ",
    "They may be gene symbols. If your pancreas gene-set JSON is ENSEMBL-based, ",
    "we need to map symbols to ENSEMBL or regenerate the gene sets as symbols."
  )
} else {
  cat("Gene IDs look like ENSEMBL IDs.\n")
}

# ----------------------------
# Step 8: Remove duplicated gene IDs if present
# ----------------------------

dup_genes <- duplicated(rownames(pancreas_sct_counts))

if (any(dup_genes)) {
  cat("Removing duplicated gene IDs:", sum(dup_genes), "\n")
  pancreas_sct_counts <- pancreas_sct_counts[!dup_genes, , drop = FALSE]
}

# ----------------------------
# Step 9: Write CSV
# ----------------------------

# fwrite cannot directly write Matrix objects.
# Convert to data.table and add gene_id as first column.
out_dt <- as.data.table(as.matrix(pancreas_sct_counts))
out_dt <- cbind(gene_id = rownames(pancreas_sct_counts), out_dt)

fwrite(out_dt, out_csv)

cat("\nWrote:", out_csv, "\n")
cat(
  "Final exported matrix:",
  nrow(pancreas_sct_counts), "genes x",
  ncol(pancreas_sct_counts), "cells\n"
)