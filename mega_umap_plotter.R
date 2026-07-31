#TODO refactor into 1 plotting function 

# ------------------------------------------------------------
# Plot UMAPs colored by metadata
#
# Generates one DimPlot for each metadata column supplied,
# combines them into a single figure, saves the combined
# figure, and also saves each individual plot.
#
# The figure size is chosen procedurally based on the number
# of plots and the number of columns, so the output is less
# likely to look stretched or squished.
#
# Examples of metadata:
#   - tissue
#   - slide
#   - orig.ident
#   - seurat_cluster.projected
# ------------------------------------------------------------
plot_umap_metadata <- function(
    object,
    variables,
    save_dir,
    reduction = "full.umap.sketch",
    ncol = NULL,
    panel_size = 6,
    individual_size = 7
) {
  
  # ----------------------------------------------------------
  # Decide how many columns to use.
  # If ncol is not supplied, choose a roughly square layout.
  # ----------------------------------------------------------
  nplots <- length(variables)
  
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(nplots))
  }
  
  nrow <- ceiling(nplots / ncol)
  
  # ----------------------------------------------------------
  # Generate one DimPlot per metadata variable.
  # Each plot is stored in a list so it can be combined and
  # saved individually later.
  # ----------------------------------------------------------
  plots <- lapply(variables, function(var) {
    
    DimPlot(
      object,
      reduction = reduction,
      group.by = var
    ) +
      ggtitle(var)
    
  })
  
  names(plots) <- variables
  
  # ----------------------------------------------------------
  # Combine all plots into a single patchwork figure.
  # ----------------------------------------------------------
  combined_plot <- wrap_plots(
    plots,
    ncol = ncol
  )
  
  # ----------------------------------------------------------
  # Save the combined figure.
  # ----------------------------------------------------------
  save_plot(
    plot_object = combined_plot,
    file_stub = "metadata_umaps",
    save_dir = save_dir,
    width = ncol * panel_size,
    height = nrow * panel_size
  )
  
  # ----------------------------------------------------------
  # Save each individual plot.
  # ----------------------------------------------------------
  for (var in variables) {
    
    save_plot(
      plot_object = plots[[var]],
      file_stub = paste0("metadata_", var),
      save_dir = save_dir,
      width = individual_size,
      height = individual_size
    )
    
  }
  
  return(combined_plot)
  
}



# ------------------------------------------------------------
# Plot UMAPs colored by continuous features
#
# Uses FeaturePlot to visualize:
#   - QC metrics
#   - gene expression
#   - continuous metadata stored as features
#
# Saves both a combined figure and each individual plot.
# Figure size is chosen procedurally from the number of plots.
# ------------------------------------------------------------
plot_umap_features <- function(
    object,
    features,
    save_dir,
    reduction = "full.umap.sketch",
    ncol = NULL,
    panel_size = 6,
    individual_size = 7
) {
  
  # ----------------------------------------------------------
  # Decide how many columns to use.
  # If ncol is not supplied, choose a roughly square layout.
  # ----------------------------------------------------------
  nplots <- length(features)
  
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(nplots))
  }
  
  nrow <- ceiling(nplots / ncol)
  
  # ----------------------------------------------------------
  # Generate one FeaturePlot per feature.
  # Each plot is stored in a list so it can be combined and
  # saved individually later.
  # ----------------------------------------------------------
  plots <- lapply(features, function(feature) {
    
    FeaturePlot(
      object,
      reduction = reduction,
      features = feature
    ) +
      ggtitle(feature)
    
  })
  
  names(plots) <- features
  
  # ----------------------------------------------------------
  # Combine all plots into a single patchwork figure.
  # ----------------------------------------------------------
  combined_plot <- wrap_plots(
    plots,
    ncol = ncol
  )
  
  # ----------------------------------------------------------
  # Save the combined figure.
  # ----------------------------------------------------------
  save_plot(
    plot_object = combined_plot,
    file_stub = "feature_umaps",
    save_dir = save_dir,
    width = ncol * panel_size,
    height = nrow * panel_size
  )
  
  # ----------------------------------------------------------
  # Save each individual plot.
  # ----------------------------------------------------------
  # for (feature in features) {
  #   
  #   save_plot(
  #     plot_object = plots[[feature]],
  #     file_stub = paste0("feature_", feature),
  #     save_dir = save_dir,
  #     width = individual_size,
  #     height = individual_size
  #   )
  #   
  # }
  
  return(combined_plot)
  
}



# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
# 
# plot_umap_metadata(
#   mega_obj,
#   variables = c(
#     "slide",
#     "tissue",
#     "orig.ident",
#     "seurat_cluster.projected"
#   ),
#   save_dir = file.path(sample_output_dir, "UMAPs", "Metadata")
# )
# 
# plot_umap_features(
#   mega_obj,
#   features = c(
#     "nCount_Spatial.008um",
#     "nFeature_Spatial.008um",
#     "percent.mt"
#   ),
#   save_dir = file.path(sample_output_dir, "UMAPs", "QC")
# )
# 
# plot_umap_features(
#   mega_obj,
#   features = c(
#     "EPCAM",
#     "KRT8",
#     "KRT19",
#     "MUC1"
#   ),
#   save_dir = file.path(sample_output_dir, "UMAPs", "Markers")
# )