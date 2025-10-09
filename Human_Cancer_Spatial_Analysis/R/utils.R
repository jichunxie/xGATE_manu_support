# Helper function to parse domain and cluster numbers from a string
# (This was the internal .extract_domain_cluster function)
extract_domain_cluster <- function(x) {
  x_chr <- as.character(x)
  has_unknown <- grepl("unknown", x_chr, ignore.case = TRUE)
  nums <- stringr::str_extract_all(x_chr, "\\d+")
  d <- suppressWarnings(as.integer(sapply(nums, function(v) if (length(v) >= 1) v[1] else NA_character_)))
  c <- suppressWarnings(as.integer(sapply(nums, function(v) if (length(v) >= 2) v[2] else NA_character_)))
  d[has_unknown] <- NA_integer_
  data.frame(domain_num = d, cluster_num = c)
}


# Convert data frame names to lower case
to_lower_names <- function(df) {
  names(df) <- tolower(names(df))
  df
}

# Find the first column name present in a data frame from a list of choices
first_present <- function(cols, choices) {
  hit <- choices[choices %in% cols]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# Create a unique key for subdomains
make_subdomain_key <- function(df, domain_col, cluster_col, subdomain_col) {
  if (!is.na(subdomain_col) && subdomain_col %in% names(df)) {
    return(as.character(df[[subdomain_col]]))
  }
  if (is.na(domain_col) || is.na(cluster_col) || !domain_col %in% names(df) || !cluster_col %in% names(df)) {
    stop("Cannot construct subdomain key: missing both 'subdomain' and ('domain','cluster') columns.")
  }
  d <- as.character(df[[domain_col]])
  c <- as.character(df[[cluster_col]])
  paste0("D", d, "_", c)
}

# Sanitize a correlation matrix by replacing non-finite values with 0
sanitize_cormat <- function(cormat) {
  cormat[is.nan(cormat) | is.infinite(cormat) | is.na(cormat)] <- 0
  cormat
}


#' Map Detailed Cell Types to General Categories
#'
#' Takes a vector of detailed cell types and maps them to broader categories
#' based on a provided mapping list.
#'
#' @param detailed_cell_types A character vector of cell types to map.
#' @param cell_type_map A named list where names are general categories and
#'   values are vectors of detailed cell types belonging to that category.
#' @return A character vector of the same length with general cell type names.
map_to_general <- function(detailed_cell_types, cell_type_map) {
  # Create a reverse map for faster lookup
  reverse_map <- stack(cell_type_map)
  names(reverse_map) <- c("general", "detailed")
  
  # Match detailed types to the general category
  general_types <- reverse_map$general[match(detailed_cell_types, reverse_map$detailed)]
  
  # Assign "Other" to any types that were not found in the map
  general_types[is.na(general_types)] <- "Other"
  
  return(general_types)
}

#' Sanitize a String for Use in a Filename
#'
#' Replaces any characters that are not alphanumeric, underscore, or hyphen
#' with an underscore.
#'
#' @param x The input string.
#' @return A sanitized string suitable for a filename.
sanitize <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}



# ... (at the end of the file) ...

#' Flatten a List-Column to a Vector
#'
#' Safely converts a list-column within a data frame to a simple vector.
#'
#' @param x The list-column to flatten.
#' @param to The target type ("character" or "integer").
#' @return A flattened vector of the specified type.
flatten_one <- function(x, to = c("character","integer")) {
  to <- match.arg(to)
  if (!is.list(x)) return(x) # Return as-is if not a list
  
  .safe_convert <- function(z, type) {
    if (length(z) == 0 || is.null(z) || all(is.na(z))) {
      return(if (type == "character") NA_character_ else NA_integer_)
    }
    if (type == "character") as.character(z[[1]]) else as.integer(z[[1]])
  }
  
  vapply(x, .safe_convert, FUN.VALUE = if (to == "character") character(1) else integer(1), type = to)
}


#' Resolve Domain Number from a Data Frame
#'
#' Finds and extracts a numeric domain identifier from common column names
#' like 'IRIS_domain', 'domain', or 'subdomain'.
#'
#' @param df The input data frame.
#' @return An integer vector of domain numbers.
resolve_domain_num <- function(df) {
  # Prefer 'IRIS_domain' if present and numeric-like; else parse from 'domain' or 'subdomain'
  if ("IRIS_domain" %in% names(df)) {
    dn <- df$IRIS_domain
    dn <- flatten_one(dn, "integer")
    # If it's not integer, parse digits
    if (is.character(dn)) dn <- suppressWarnings(as.integer(stringr::str_extract(dn, "\\d+")))
    return(dn)
  }
  col <- NULL
  if ("domain" %in% names(df)) col <- "domain"
  if (is.null(col) && "subdomain" %in% names(df)) col <- "subdomain"
  if (is.null(col)) stop("No 'domain', 'IRIS_domain', or 'subdomain' column found.")
  dom_chr <- flatten_one(df[[col]], "character")
  unknown <- grepl("unknown", dom_chr, ignore.case = TRUE)
  dn <- suppressWarnings(as.integer(stringr::str_extract(dom_chr, "\\d+")))
  dn[unknown] <- NA_integer_
  dn
}


#' Resolve General Cell Type from a Data Frame
#'
#' Finds and extracts cell type labels from common column names like
#' 'general_cell_type', 'celltype', 'annotation', etc.
#'
#' @param df The input data frame.
#' @return A character vector of cell type labels.
resolve_general_cell_type <- function(df) {
  # Use general_cell_type if present; else try alternates
  if ("general_cell_type" %in% names(df)) {
    g <- df$general_cell_type
    return(flatten_one(g, "character"))
  }
  alt_cols <- c("cell_type", "celltype", "predicted_celltype", "predicted_cell_type", "annotation", "broad_type", "label", "labels")
  hit <- alt_cols[alt_cols %in% names(df)]
  if (length(hit) == 0) stop("No cell type column found (looked for general_cell_type, ", paste(alt_cols, collapse=", "), ").")
  flatten_one(df[[hit[1]]], "character")
}
