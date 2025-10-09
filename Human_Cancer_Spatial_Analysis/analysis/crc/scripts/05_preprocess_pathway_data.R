# --- 1. SETUP ---
# Load configuration and utility functions
source("scripts/00_config.R")
source("R/utils.R")

# --- 2. LOAD RAW DATA ---
message("Loading raw data...")
pathway_res <- read.csv(PATHWAY_FILE, check.names = FALSE, stringsAsFactors = FALSE)
sp_obj <- readRDS(SPATIAL_FILE)
spatial_df <- as.data.frame(sp_obj$spatial_location_list[[1]])

# --- 3. PREPROCESS PATHWAY DATA ---
message("Processing pathway results...")
# Ensure column names are clean (e.g., "p-value" -> "p.value")
names(pathway_res) <- make.names(names(pathway_res))
pathway_res <- pathway_res %>%
  mutate(
    domain_num  = as.integer(str_extract(as.character(domain), "\\d+")),
    cluster_num = as.integer(str_extract(as.character(cluster), "\\d+")),
    p.value     = as.numeric(p.value)
  ) %>%
  filter(!is.na(domain_num), !is.na(cluster_num))

# Aggregate pathway results: take the minimum p-value for any duplicate domain/cluster pairs
pathway_agg <- pathway_res %>%
  group_by(pathway, domain_num, cluster_num) %>%
  summarize(p.value = min(p.value, na.rm = TRUE), .groups = "drop")

# --- 4. PREPROCESS SPATIAL DATA ---
message("Processing spatial data...")
# Parse domain/cluster from the 'subdomain' column
parsed_coords <- extract_domain_cluster(spatial_df$subdomain)
spatial_df$domain_num <- parsed_coords$domain_num
spatial_df$cluster_num <- parsed_coords$cluster_num

# Filter out spots with no domain information
n_before <- nrow(spatial_df)
spatial_df <- spatial_df %>% filter(!is.na(domain_num))
n_after <- nrow(spatial_df)
message(sprintf("Filtered out unknown domains: kept %d of %d spots (%.1f%%).",
                n_after, n_before, 100 * n_after / max(1, n_before)))

# --- 5. COMBINE AND SAVE ---
# We will combine pathway p-values with spatial data in the plotting scripts
# For now, let's save the key processed components
processed_data <- list(
  spatial_df = spatial_df,
  pathway_agg = pathway_agg
)

saveRDS(processed_data, file = file.path(PROCESSED_DATA_DIR, "processed_colon_data.rds"))
message("Processed data saved to: ", file.path(PROCESSED_DATA_DIR, "processed_colon_data.rds"))