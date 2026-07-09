# ============================================================
# PCA VARIANCE PLOT
#
# Goal
# ----
# Plot the proportion of variance explained by each principal
# component and the cumulative variance explained.
#
# Responsibilities
# ----------------
# * Validate the requested PCA reduction.
# * Compute variance explained.
# * Generate the variance plots.
# * Automatically save the figure as:
#       <save_dir>/variance.pdf
# * Return the combined plot.
#
# This function assumes:
# * `save_dir` is the specific PCA variance figure directory
#   (typically something like .../figures/variance/).
# * The PCA reduction already exists in the object.
# ============================================================

plot_pca_variance <- function(
    object,
    save_dir,
    reduction = "pca.sketch",
    dims_max = NULL,
    plot_title = "Sketch PCA",
    width = 10,
    height = 5
) {
  
  # ----------------------------------------------------------
  # Validate inputs
  # ----------------------------------------------------------
  
  if (missing(save_dir) || is.null(save_dir)) {
    stop(
      "save_dir must be supplied.",
      call. = FALSE
    )
  }
  
  if (!dir.exists(save_dir)) {
    stop(
      "save_dir does not exist:\n",
      save_dir,
      call. = FALSE
    )
  }
  
  if (!reduction %in% Reductions(object)) {
    stop(
      "Reduction '",
      reduction,
      "' not found.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------------------------
  # Extract standard deviations from the PCA reduction
  # ----------------------------------------------------------
  
  stdev <- object[[reduction]]@stdev
  var_explained <- stdev^2 / sum(stdev^2)
  
  pca_df <- data.frame(
    PC = seq_along(stdev),
    Variance = var_explained,
    Cumulative = cumsum(var_explained)
  )
  
  # ----------------------------------------------------------
  # Variance explained by each PC
  # ----------------------------------------------------------
  
  p_pca_var <- ggplot(
    pca_df,
    aes(PC, Variance)
  ) +
    geom_line() +
    geom_point() +
    scale_x_continuous(
      breaks = seq(0, max(pca_df$PC), by = 5)
    ) +
    scale_y_continuous(
      breaks = seq(
        0,
        max(pca_df$Variance),
        by = 0.01
      ),
      labels = scales::percent
    ) +
    labs(
      title = paste(
        plot_title,
        "Variance Explained"
      ),
      x = "Principal Component",
      y = "Proportion of Variance"
    )
  
  if (!is.null(dims_max)) {
    p_pca_var <- p_pca_var +
      geom_vline(
        xintercept = dims_max,
        linetype = "dashed",
        color = "red"
      )
  }
  
  # ----------------------------------------------------------
  # Cumulative variance explained
  # ----------------------------------------------------------
  
  p_pca_cumulative <- ggplot(
    pca_df,
    aes(PC, Cumulative)
  ) +
    geom_line() +
    geom_point() +
    scale_x_continuous(
      breaks = seq(0, max(pca_df$PC), by = 5)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent
    ) +
    labs(
      title = paste(
        "Cumulative",
        plot_title,
        "Variance"
      ),
      x = "Principal Component",
      y = "Cumulative Proportion of Variance"
    )
  
  if (!is.null(dims_max)) {
    p_pca_cumulative <- p_pca_cumulative +
      geom_vline(
        xintercept = dims_max,
        linetype = "dashed",
        color = "red"
      )
  }
  
  # ----------------------------------------------------------
  # Combine plots
  # ----------------------------------------------------------
  
  variance_plot <-
    p_pca_var |
    p_pca_cumulative
  
  # ----------------------------------------------------------
  # Save the plot using the canonical filename
  # ----------------------------------------------------------
  
  save_plot(
    plot_object = variance_plot,
    file_stub = "variance",
    save_dir = save_dir,
    width = width,
    height = height
  )
  
  # ----------------------------------------------------------
  # Return the combined plot
  # ----------------------------------------------------------
  
  return(variance_plot)
}



# ============================================================
# PCA DIAGNOSTIC PLOTS
#
# Goal
# ----
# Generate and automatically save the standard diagnostic plots
# associated with a completed sketch PCA analysis.
#
# Responsibilities
# ----------------
# * Create PCA loading plots.
# * Create an elbow plot.
# * Create a variance explained plot.
# * Automatically save every figure into canonical locations.
# * Return the ggplot objects.
#
# Assumptions
# -----------
# * `run_sketch_pca()` has already been executed.
# * The requested PCA reduction already exists.
# * `pca_dir` is the PCA analysis directory, e.g.
#     .../dimension_reduction/pca/dims50/
# ============================================================

plot_pca_diagnostics <- function(
    object,
    pca_dir,
    reduction = "pca.sketch",
    dims = 1:6,
    ndims = 50
) {
  
  # ----------------------------------------------------------
  # Validate inputs
  # ----------------------------------------------------------
  
  if (missing(pca_dir) || is.null(pca_dir)) {
    stop(
      "pca_dir must be supplied.",
      call. = FALSE
    )
  }
  
  if (!dir.exists(pca_dir)) {
    stop(
      "pca_dir does not exist:\n",
      pca_dir,
      call. = FALSE
    )
  }
  
  if (!reduction %in% Reductions(object)) {
    stop(
      "Reduction '",
      reduction,
      "' not found in object.",
      call. = FALSE
    )
  }
  
  # ----------------------------------------------------------
  # Resolve figure output directories
  # ----------------------------------------------------------
  figures_output_dir <- get_figures_dir(pca_dir)
    
  
  # ----------------------------------------------------------
  # PCA loading plot
  # ----------------------------------------------------------
  
  
  
  cat("figures_output_dir:\n")
  print(figures_output_dir)
  
  cat("\nDirectory exists before save:\n")
  print(dir.exists(figures_output_dir))
  
  
  loading_plot <- VizDimLoadings(
    object = object,
    dims = dims,
    reduction = reduction
  )
  
  save_plot(
    plot_object = loading_plot,
    file_stub = "loadings",
    save_dir = figures_output_dir,
    width = 10,
    height = 20
  )
  
  # ----------------------------------------------------------
  # Elbow plot
  # ----------------------------------------------------------
  
  elbow_plot <- ElbowPlot(
    object = object,
    reduction = reduction,
    ndims = ndims
  ) +
    ggtitle("Variance Explained by Sketch Principal Components")
  
  save_plot(
    plot_object = elbow_plot,
    file_stub = "elbow",
    save_dir = figures_output_dir,
    width = 10,
    height = 5
  )
  
  # ----------------------------------------------------------
  # Variance explained plot
  # ----------------------------------------------------------
  
  variance_plot <- plot_pca_variance(
    object = object,
    save_dir = figures_output_dir,
    reduction = reduction,
    dims_max = NULL,
    plot_title = "Sketch PCA",
    width = 10,
    height = 5
  )
  
  # ----------------------------------------------------------
  # Finished
  # ----------------------------------------------------------
  
  message("----------------------------------------")
  message("PCA diagnostic plots saved:")
  message("  To : ", figures_output_dir)
  message("----------------------------------------")
  
  return(
    list(
      loading_plot = loading_plot,
      elbow_plot = elbow_plot,
      variance_plot = variance_plot
    )
  )
}