# --- 1. SETUP ---
# This script calculates cell type proportions and log2 odds ratio enrichment
# for each spatial domain in the prostate dataset.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Source shared utility functions and dataset-specific config
source("../../R/utils.R")
source("00_config_prostate.R")

# --- 2. LOAD & PREPARE DATA ---
message("Loading annotated prostate data: ", SP_INPUT_ANNOTATED_FILE)
sp_obj <- readRDS(SP_INPUT_ANNOTATED_FILE)
loc_df <- sp_obj$spatial_location_list[[1]]

# Use helper functions to robustly extract domain and cell type info
meta <- loc_df %>%
  mutate(
    domain_num = resolve_domain_num(loc_df),
    general_cell_type = resolve_general_cell_type(loc_df)
  ) %>%
  filter(!is.na(domain_num), !is.na(general_cell_type), nzchar(general_cell_type)) %>%
  # Filter out specific cell types if needed
  filter(!grepl("sperm", general_cell_type, ignore.case = TRUE))

if (nrow(meta) == 0) stop("No data remaining after preparing metadata.")

# --- 3. CALCULATE CELL TYPE PROPORTIONS ---
message("Calculating cell type proportions per domain...")
prop_df <- meta %>%
  dplyr::count(domain_num, general_cell_type, name = "n") %>%
  group_by(domain_num) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Save the proportions table
readr::write_csv(prop_df, PROSTATE_CELLTYPE_PROPS_TABLE)
message("Proportions table saved to: ", PROSTATE_CELLTYPE_PROPS_TABLE)

# --- 4. CALCULATE LOG2 ODDS RATIO ENRICHMENT ---
message("Calculating log2 odds ratio for enrichment...")
# Pre-calculate totals for efficiency
in_counts  <- dplyr::count(meta, domain_num, general_cell_type, name = "a")
tot_domain <- dplyr::count(meta, domain_num, name = "n_domain")
tot_type   <- dplyr::count(meta, general_cell_type, name = "n_type")
n_total    <- nrow(meta)

# Create a complete grid of all domain/cell type combinations
all_combos <- tidyr::expand_grid(
  domain_num = sort(unique(meta$domain_num)),
  general_cell_type = sort(unique(meta$general_cell_type))
)

log_odds_df <- all_combos %>%
  left_join(in_counts, by = c("domain_num", "general_cell_type")) %>%
  left_join(tot_domain, by = "domain_num") %>%
  left_join(tot_type, by = "general_cell_type") %>%
  mutate(
    a = ifelse(is.na(a), 0, a), # if a group has 0 cells, 'a' will be NA
    b = n_domain - a,
    c = n_type - a,
    d = n_total - n_domain - c,
    # Add pseudocount and calculate log2 odds
    log2_odds = log2(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
  ) %>%
  select(general_cell_type, domain_num, log2_odds)

# Save the enrichment table
readr::write_csv(log_odds_df, PROSTATE_CELLTYPE_LOGODDS_TABLE)
message("Log odds table saved to: ", PROSTATE_CELLTYPE_LOGODDS_TABLE)
message("Cell type enrichment analysis complete.")