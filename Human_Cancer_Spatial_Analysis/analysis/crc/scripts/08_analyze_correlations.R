# --- 1. SETUP ---
# Load configuration, utilities, and plotting functions
source("scripts/00_config.R")
source("R/utils.R")
source("R/correlation_heatmap_functions.R")

# --- 2. LOAD AND FILTER DATA ---
message("Reading pathway data: ", PATHWAY_FILE)
df_raw <- readr::read_csv(PATHWAY_FILE, show_col_types = FALSE, guess_max = 1e6)
df <- to_lower_names(df_raw)

# Dynamically find the correct column names
pathway_col <- first_present(names(df), c("pathway", "term", "gs_name"))
pval_col    <- first_present(names(df), c("p-value", "p.value", "pvalue", "padj"))
z_col       <- first_present(names(df), c("z-score", "zscore", "z_score"))
domain_col  <- first_present(names(df), c("domain"))
cluster_col <- first_present(names(df), c("cluster"))
subdomain_col <- if ("subdomain" %in% names(df)) "subdomain" else NA_character_

if (is.na(pathway_col)) stop("Could not find a pathway column.")

# Filter for selected pathways from the config file
message("Filtering for ", length(SELECTED_PATHWAYS), " selected pathways.")
df <- df %>%
  filter(if (CASE_INSENSITIVE_MATCH) tolower(.data[[pathway_col]]) %in% tolower(SELECTED_PATHWAYS) else .data[[pathway_col]] %in% SELECTED_PATHWAYS)

if (dplyr::n_distinct(df[[pathway_col]]) < 2) {
  stop("Fewer than 2 selected pathways were found in the data. Cannot compute correlations.")
}

# Create the subdomain key for pivoting
subdomain_key <- make_subdomain_key(df, domain_col, cluster_col, subdomain_col)


# --- 3. Z-SCORE BASED CORRELATION ---
if (!is.na(z_col)) {
  message("Calculating correlations based on z-scores...")
  zvec <- as.numeric(df[[z_col]])
  mat_z <- build_matrix(df, pathway_col, subdomain_key, zvec, what = "zscore", min_non_na_cols = CORR_MIN_NON_NA)

  if (nrow(mat_z) >= 2) {
    cor_z <- cor(t(mat_z), use = "pairwise.complete.obs", method = CORR_METHOD)
    cor_z <- sanitize_cormat(cor_z)
    plot_and_save_heatmap(
      cor_mat = cor_z,
      output_dir = CORR_PLOT_DIR,
      tag = "zscore",
      corr_method = CORR_METHOD,
      selected_pathways = SELECTED_PATHWAYS,
      case_insensitive_match = CASE_INSENSITIVE_MATCH,
      label_max_chars = HEATMAP_LABEL_MAX_CHARS
    )
  }
} else {
  message("Z-score column not found. Skipping z-score based correlations.")
}

message("Correlation analysis complete.")