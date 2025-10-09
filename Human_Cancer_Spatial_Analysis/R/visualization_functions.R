# R/visualization_functions.R
#
# Contains reusable, high-level functions for creating the main visualizations
# for the project, including pathway maps and subdomain maps.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(cowplot)
  library(grid)
  library(patchwork)
  library(viridis)
})


# =========================================================================
# ==                          PATHWAY VISUALIZATION                      ==
# =========================================================================

#' Save a Standalone ggplot Colorbar Legend
#'
#' Extracts a legend from a ggplot object and saves it as a separate file.
#' Particularly useful for creating shared legends for multi-panel figures.
#'
#' @param p_with_legend A ggplot object that includes the legend to be saved.
#' @param file The full path for the output image file.
#' @param width The width of the output image in inches.
#' @param height The height of the output image in inches.
#' @param dpi The resolution of the output image.
#' @return The plot is saved to a file.
save_colorbar_legend <- function(p_with_legend, file, width = 3, height = 0.6, dpi = 600) {
  leg <- cowplot::get_legend(
    p_with_legend +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = ggplot2::margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
      ) +
      ggplot2::guides(
        color = ggplot2::guide_colorbar(
          direction = "horizontal", title.position = "top", title.hjust = 0.5,
          label.position = "bottom", barheight = grid::unit(3, "mm"),
          barwidth  = grid::unit(35, "mm"), ticks.colour = "black"
        )
      )
  )
  cowplot::ggsave2(filename = file, plot = cowplot::ggdraw(leg),
                   width = width, height = height, dpi = dpi, bg = "white")
}


#' Plot a Spatial Pathway p-value Map
#'
#' Generates a spatial scatter plot colored by pathway p-values.
#' Handles both full and zoomed views.
#'
#' @param spatial_df_with_pval A data frame with 'x', 'y', and 'p.value' columns.
#' @param pathway_name The name of the pathway for the plot title.
#' @param output_file The full path where the output image will be saved.
#' @param color_values The numeric vector of breaks for the color scale.
#' @param color_palette The character vector of colors for the color scale.
#' @param zoom_coords An optional list with 'xlim' and 'ylim' for a zoomed view.
#' @param ... Additional arguments passed to ggsave (e.g., plot_width, plot_height).
#' @return A ggplot object with a legend, which can be used by `save_colorbar_legend`.
plot_pathway_map <- function(spatial_df_with_pval,
                             pathway_name,
                             output_file,
                             color_values,
                             color_palette,
                             zoom_coords = NULL,
                             point_size = 0.3,
                             ...) {

  p <- ggplot(spatial_df_with_pval, aes(x = x, y = y, color = p.value)) +
    geom_point(size = point_size, shape = 15, stroke = 0, alpha = 0.95) +
    scale_color_gradientn(
      colors = color_palette,
      values = scales::rescale(color_values),
      limits = c(0, 1),
      na.value = "gray90",
      name = "p-value",
      breaks = c(0, 0.05, 1),
      labels = c("0", "0.05", "1")
    ) +
    theme_minimal(base_size = 8) +
    ggtitle(pathway_name) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none" # Legend is removed from the saved plot
    )

  coord_system <- if (!is.null(zoom_coords)) {
    coord_fixed(xlim = zoom_coords$xlim, ylim = zoom_coords$ylim, expand = FALSE)
  } else {
    coord_fixed()
  }
  p <- p + coord_system

  ggsave(filename = output_file, plot = p, bg = "white", ...)

  # Return the plot object with a legend attached, for potential saving
  return(p + theme(legend.position = "right"))
}


# =========================================================================
# ==                         SUBDOMAIN VISUALIZATION                     ==
# =========================================================================

#' Plot a Spatial Subdomain Map
#'
#' Generates a spatial map of subdomains with a custom tiled legend. Can produce
#' either a full view or a zoomed-in view based on provided coordinates.
#'
#' @param spatial_df A data frame with 'x', 'y', 'domain', and 'subdomain' columns.
#' @param output_file The full path where the final JPG image will be saved.
#' @param zoom_coords An optional list with 'xlim' and 'ylim' vectors to specify the zoom region.
#' @param ... Additional arguments passed to ggsave (e.g., width, height, dpi).
#' @return Invisibly returns the combined ggplot object.
plot_subdomain_map <- function(spatial_df,
                               output_file,
                               zoom_coords = NULL,
                               ...) {

  # --- 1. Prepare data and color map ---
  plot_df <- spatial_df %>%
    mutate(domain_subdomain_id = paste(domain, subdomain, sep = "_"))

  n_colors <- length(unique(plot_df$domain_subdomain_id))
  if (n_colors == 0) {
    warning("No data to plot. Skipping figure generation.")
    return(invisible(NULL))
  }
  color_palette <- viridis(n_colors, option = "turbo")
  color_map <- setNames(color_palette, sort(unique(plot_df$domain_subdomain_id)))

  # --- 2. Create main spatial plot ---
  coord_system <- if (!is.null(zoom_coords)) {
    coord_fixed(xlim = zoom_coords$xlim, ylim = zoom_coords$ylim, expand = FALSE)
  } else {
    coord_fixed()
  }

  spatial_plot <- ggplot(plot_df, aes(x = x, y = y, color = domain_subdomain_id)) +
    geom_point(size = 0.4, shape = 16) +
    scale_color_manual(values = color_map) +
    coord_system +
    theme_void() +
    theme(legend.position = "none")

  # --- 3. Create custom legend ---
  legend_data <- plot_df %>%
    distinct(domain, subdomain, domain_subdomain_id) %>%
    arrange(domain, subdomain) %>%
    mutate(color = color_map[domain_subdomain_id]) %>%
    group_by(domain) %>%
    mutate(subdomain_numeric = as.numeric(factor(subdomain))) %>%
    ungroup()

  legend_plot <- ggplot(legend_data, aes(x = subdomain_numeric, y = domain, fill = I(color))) +
    geom_tile(color = "white", linewidth = 0.5) +
    labs(x = "Subdomain Index", y = "Domain") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 9)
    )

  # --- 4. Combine plots and save ---
  combined_plot <- spatial_plot + legend_plot + plot_layout(widths = c(3, 1))

  ggsave(filename = output_file, plot = combined_plot, bg = "white", ...)
  
  message("Subdomain map saved to: ", output_file)
  return(invisible(combined_plot))
}


#' Create a Single Zoomed-in Spatial Plot
#'
#' A helper function to generate one ggplot object for a specific zoomed-in region.
#'
#' @param df_zoom A data frame already filtered to the zoom region.
#' @param color_map A named vector mapping domain_subdomain_id to colors.
#' @param xlim The x-axis limits for the plot's coordinate system.
#' @param ylim The y-axis limits for the plot's coordinate system.
#' @param label A character string to use as the plot title.
#' @return A ggplot object representing the zoomed-in plot.
make_zoom_plot <- function(df_zoom, color_map, xlim, ylim, label) {
  if (nrow(df_zoom) == 0) {
    # Return an empty plot if no points are in the zoom box
    return(ggplot() + theme_void() + ggtitle(glue::glue("{label}: (no data points)")))
  }
  
  ggplot(df_zoom, aes(x = x, y = y, color = domain_subdomain_id)) +
    geom_point(size = 0.55, shape = 15, stroke = 0) +
    scale_color_manual(values = color_map, drop = FALSE) +
    coord_fixed(xlim = range(xlim), ylim = range(ylim), expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    ggtitle(glue::glue("{label}"))
}