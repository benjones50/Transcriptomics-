plot_metadata_highlights_2 <- function(
    object,
    group.by,
    save_dir,
    reduction,
    pt.size = 0.02,
    selected_alpha = 1,
    overview_alpha = 0.4,
    shuffle_points = TRUE,
    shuffle_seed = 42,
    background_alpha = 0.02,
    background_color = "white",
    combine_plots = TRUE,
    save_individuals = TRUE, #TODO this is non functional, keep off
    width = 16,
    height = 10,
    dpi = 600,
    color.by = NULL,
    colors = NULL,
    show_cluster_labels = FALSE
) {
  
  
  # ==========================================================
  # Create output directory structure
  # ==========================================================
  
  if (!dir.exists(save_dir)) {
    dir.create(
      save_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  group_dir <- file.path(save_dir, paste0("plots_grouped_by_", group.by))
  
  dir.create(
    group_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (!dir.exists(group_dir)) {
    stop(
      "Could not create output directory:\n",
      group_dir,
      call. = FALSE
    )
  }
  
  # ==========================================================
  # Internal helper: sanitize filenames
  # ==========================================================
  
  sanitize_filename <- function(x) {
    x <- as.character(x)
    
    # Replace whitespace with underscores.
    x <- gsub("\\s+", "_", x)
    
    # Remove filesystem-unfriendly characters.
    x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)
    
    # Collapse repeated underscores.
    x <- gsub("_+", "_", x)
    
    # Trim leading/trailing underscores.
    x <- gsub("^_", "", x)
    x <- gsub("_$", "", x)
    
    if (nchar(x) == 0) {
      x <- "unnamed"
    }
    
    x
  }
  
  # ==========================================================
  # Internal helper: make a more readable panel title
  # ==========================================================
  
  pretty_title <- function(x) {
    x <- as.character(x)
    
    # Keep double underscores as a field separator, then make
    # the title more readable by turning them into line breaks.
    x <- gsub("__", "\n", x, fixed = TRUE)
    
    # Convert remaining underscores to spaces.
    x <- gsub("_", " ", x, fixed = TRUE)
    
    x
  }
  
  # ==========================================================
  # Internal helper: extract metadata values and levels
  # ==========================================================
  
  extract_metadata_column <- function(metadata_name) {
    metadata_values <- object[[metadata_name]][rownames(coords), 1, drop = TRUE]
    
    # Preserve the ordering of factor levels.
    # Otherwise, fall back to alphabetical ordering for character metadata.
    if (is.factor(metadata_values)) {
      metadata_levels <- levels(metadata_values)
    } else {
      #orders numerically groupings
      unique_values <- unique(metadata_values)
      if (all(grepl("^[0-9]+$", unique_values))) {
        metadata_levels <- as.character(
          sort(as.integer(unique_values))
        )
      } else {
        metadata_levels <- sort(unique_values)
      }
    }
    
    # Convert metadata values to character so downstream
    # comparisons and plotting behave consistently regardless of
    # whether the original metadata was stored as a factor or
    # character vector.
    metadata_values <- as.character(metadata_values)
    
    # Give missing metadata values an explicit label so they can
    # be plotted rather than silently dropped.
    metadata_values[is.na(metadata_values)] <- "NA"
    
    list(
      values = metadata_values,
      levels = as.character(metadata_levels)
    )
  }
  

  
  # ==========================================================
  # Extract reduction coordinates and metadata
  # ==========================================================
  # Prepare plotting data
  #
  # Goal
  # ----
  # Build a standalone data frame containing everything needed
  # for plotting.
  #
  # Steps
  # -----
  # * Extract the first two dimensions from the requested
  #   dimensional reduction.
  # * Verify that at least two dimensions are available for
  #   plotting.
  # * Identify the x and y coordinate column names.
  # * Retrieve the requested metadata column for every cell/bin.
  # * Convert metadata values to character so downstream code
  #   handles factors and characters consistently.
  # * Replace missing metadata values with an explicit "NA"
  #   label so they can be displayed as their own group.
  # * Store the metadata alongside the coordinates so every
  #   plotting function only needs to work with a single data
  #   frame.
  # ==========================================================
  
  #    UMAP_1 | UMAP_2 | tissue_name / group.by
  #    -----: | -----: | -----------
  #      -5.1 |    2.8 | WT
  #      -4.9 |    3.1 | WT
  #       1.3 |   -6.4 | cKO
  #
  
  coords <- as.data.frame(Embeddings(object, reduction = reduction))
  
  if (ncol(coords) < 2) {
    stop(
      "Reduction '",
      reduction,
      "' must have at least two dimensions.",
      call. = FALSE
    )
  }
  
  x_col <- colnames(coords)[1]
  y_col <- colnames(coords)[2]
  
  group_info <- extract_metadata_column(group.by)
  group_values <- group_info$values
  group_levels <- group_info$levels
  
  # Store the metadata alongside the reduction coordinates so
  # all subsequent plotting operates on a single data frame.
  coords[[group.by]] <- group_values
  
  if (!is.null(color.by)) {
    color_info <- extract_metadata_column(color.by)
    color_values <- color_info$values
    color_levels <- color_info$levels
    
    # Store the second metadata column if color.by is enabled.
    coords[[color.by]] <- color_values
  } else {
    color_values <- group_values
    color_levels <- group_levels
  }
  
  
  
  # ----------------------------------------------------------
  # Build the plotting color palette
  #
  # If a named color vector is supplied, use it directly.
  # Otherwise, generate a Seurat discrete palette and assign
  # colors to metadata values in the order of color_levels.
  # ----------------------------------------------------------
  if (is.null(colors)) {
    
    plot_colors <- Seurat::DiscretePalette(
      length(color_levels)
    )
    
    names(plot_colors) <- color_levels
    
  } else {
    
    plot_colors <- colors
  }
  
  
  # Legend title should reflect whichever metadata column is
  # currently driving point colors.
  color_legend_title <- if (is.null(color.by)) group.by else color.by
  
  # ==========================================================
  # Randomize plotting order of individual points
  #
  # Randomizing the plotting order prevents one metadata group
  # from consistently being drawn on top of another while
  # preserving every cell's coordinates and metadata.
  # ==========================================================
  if (shuffle_points) {
    
    if (!is.null(shuffle_seed)) {
      set.seed(shuffle_seed)
    }
    
    coords <- coords[
      sample(nrow(coords)),
      ,
      drop = FALSE
    ]
    
  }
  
  if (length(group_levels) == 0) {
    stop(
      "No values were found in metadata column '",
      group.by,
      "'.",
      call. = FALSE
    )
  }
  
  if (!is.null(color.by) && length(color_levels) == 0) {
    stop(
      "No values were found in metadata column '",
      color.by,
      "'.",
      call. = FALSE
    )
  }
  
  # ==========================================================
  # Internal helper: build one plot for one selected group
  # ==========================================================
  
  make_highlight_plot <- function(selected_group) {
    
    selected_idx <- coords[[group.by]] == selected_group
    
    background_df <- coords[!selected_idx, , drop = FALSE]
    selected_df <- coords[selected_idx, , drop = FALSE]
    
    # ----------------------------------------------------------
    # Build the layers that differ between the two cases
    # ----------------------------------------------------------
    if (is.null(color.by)) {
      
      selected_layer <- ggplot2::geom_point(
        data = selected_df,
        ggplot2::aes(
          x = .data[[x_col]],
          y = .data[[y_col]]
        ),
        color = plot_colors[[selected_group]],
        alpha = selected_alpha,
        size = pt.size
      )
      
      color_scale <- NULL
      
      extra_theme <- ggplot2::theme()
      
    } else {
      
      selected_layer <- ggplot2::geom_point(
        data = selected_df,
        ggplot2::aes(
          x = .data[[x_col]],
          y = .data[[y_col]],
          color = .data[[color.by]]
        ),
        alpha = selected_alpha,
        size = pt.size,
        show.legend = FALSE
      )
      
      color_scale <- scale_color_manual(
        values = plot_colors,
        drop = FALSE
      )
      
      extra_theme <- ggplot2::theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    }
    
    # ----------------------------------------------------------
    # Build the common plot once
    # ----------------------------------------------------------
    highlight_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = background_df,
        ggplot2::aes(
          x = .data[[x_col]],
          y = .data[[y_col]]
        ),
        color = background_color,
        alpha = background_alpha,
        size = pt.size
      ) +
      selected_layer
    
    # Add the color scale only when color.by is used
    if (!is.null(color_scale)) {
      highlight_plot <- highlight_plot + color_scale
    }
    
    # Add the shared finishing layers
    highlight_plot <- highlight_plot +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = paste0(pretty_title(selected_group))
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(
          t = 5,
          r = 5,
          b = 5,
          l = 5
        )
      ) +
      extra_theme
    
    highlight_plot
  }
  
  
  # ==========================================================
  # Overview plot
  # ==========================================================
  message("making overview plot")
  
  # ----------------------------------------------------------
  # Choose the metadata column and legend settings once
  # ----------------------------------------------------------
  if (is.null(color.by)) {
    color_column <- group.by
    legend_title <- group.by
  } else {
    color_column <- color.by
    legend_title <- color_legend_title
  }
  
  # ----------------------------------------------------------
  # Build overview plot
  # ----------------------------------------------------------
  overview_plot <- ggplot2::ggplot(coords) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[x_col]],
        y = .data[[y_col]],
        color = .data[[color_column]]
      ),
      size = pt.size,
      alpha = overview_alpha
    ) +
    ggplot2::scale_color_manual(
      values = plot_colors,
      breaks = color_levels,
      drop = FALSE,
      guide = ggplot2::guide_legend(
        override.aes = list(
          size = 5,
          alpha = 1
        )
      )
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0(group.by, ": overall"),
      color = legend_title
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.size = grid::unit(0.5, "cm")
    )
  
  # ----------------------------------------------------------
  # Add cluster labels to the overview plot
  #
  # Labels are placed at the median coordinate of each group
  # represented by color_column. Using the median makes label
  # positions less sensitive to outlying cells/bins.
  # ----------------------------------------------------------
  if (show_cluster_labels) {
    
    label_positions <- coords %>%
      dplyr::filter(
        !is.na(.data[[color_column]])
      ) %>%
      dplyr::group_by(
        cluster = .data[[color_column]]
      ) %>%
      dplyr::summarise(
        x = stats::median(.data[[x_col]]),
        y = stats::median(.data[[y_col]]),
        .groups = "drop"
      )
    
    overview_plot <- overview_plot +
      ggplot2::geom_label(
        data = label_positions,
        mapping = ggplot2::aes(
          x = x,
          y = y,
          label = cluster
        ),
        inherit.aes = FALSE,
        size = 4.5,
        fill = "white",
        color = "black",
        fontface = "bold",
        alpha = 0.5,
        linewidth = 0,
        show.legend = FALSE
      )
  }
  
  if (save_individuals) {
  # Save the overview plot.
  message("saving overview plot")
  save_plot(
    plot_object = overview_plot,
    file_stub = "overall",
    save_dir = group_dir,
    width = width,
    height = height,
    dpi = dpi
  )
  }
  # ==========================================================
  # Build and save one highlighted panel per unique group
  # ==========================================================
  
  plots <- list()
  plots$overall <- overview_plot
  
  # Total number of highlighted panels that will be created.
  n_groups <- length(group_levels)
  
  for (i in seq_along(group_levels)) {
    
    selected_group <- group_levels[i]
    
    message(
      "Creating highlight plot ",
      i, " of ", n_groups, ": ", selected_group,
      "\n"
    )
    
    panel_plot <- make_highlight_plot(selected_group)
    
    if (save_individuals) {
      message(
        "Saving highlight plot ",
        i, " of ",
        n_groups
      )
      
      save_plot(
        plot_object = panel_plot,
        file_stub = sanitize_filename(selected_group),
        save_dir = group_dir,
        width = width,
        height = height,
        dpi = dpi
      )
    }
    
    plots[[selected_group]] <- panel_plot
  }
  
  # ==========================================================
  # Combine plots if requested
  # ==========================================================
  
  if (combine_plots) {
    
    message("Combining plots")
    
    # ========================================================
    # Build the panel of highlighted groups.
    # ========================================================
    
    #makes # of plots per row of right panel
    n_groups <- length(group_levels)
    if (n_groups <= 3) {
      plots_per_row <- n_groups
    } else {
      plots_per_row <- ceiling(sqrt(n_groups))
    }
    
    print(names(plots)) #debug
    print(group_levels)

    
    right_panel <- patchwork::wrap_plots(
      plots[group_levels],
      ncol = plots_per_row
    )
    
    # ========================================================
    # Create the final layout with a large overview plot on
    # the left and the highlighted panels on the right.
    # ========================================================
    
    n_rows <- ceiling(n_groups / plots_per_row)
    
    overview_width <- 2.0

    
    combined_plot <-
      overview_plot +
      right_panel +
      patchwork::plot_layout(
        widths = c(overview_width, plots_per_row)
      )
    # if (save_individuals) {   
    #   save_plot(
    #     plot_object = right_panel,
    #     file_stub = "right_panel",
    #     save_dir = group_dir,
    #     width = width,
    #     height = height,
    #     dpi = dpi
    #   )
    # }
    
    rm(right_panel)
    
    
    # if (save_individuals) {
    #   save_plot(
    #     plot_object = overview_plot,
    #     file_stub = "overview_plot",
    #     save_dir = group_dir,
    #     width = width,
    #     height = height,
    #     dpi = dpi
    #   )
    # }
    
    rm(overview_plot)
    
    
    save_plot(
      plot_object = combined_plot,
      file_stub = if (is.null(color.by)) "combined" else "combined_colored",
      save_dir = group_dir,
      width = 10 + 2.0 * plots_per_row,
      height = 5 + 2.0 * n_rows,
      dpi = dpi
    )
    
    plots$combined <- combined_plot
  }
  
  # ==========================================================
  # Message and return
  # ==========================================================
  
  message("----------------------------------------")
  message("Saved metadata highlight plots to:")
  message("  ", group_dir)
  message("----------------------------------------")
  
 
  return()
  
}