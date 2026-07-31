# ------------------------------------------------------------
# Extract DimPlot color palette
#
# Reconstruct the discrete color palette used by Seurat's
# DimPlot() for a given metadata column. This is useful for
# ensuring custom ggplot figures use the exact same colors as
# Seurat visualizations.
#
# The function builds a temporary DimPlot, extracts the colors
# assigned to each plotted point, and recovers the one-to-one
# mapping between metadata values and their assigned colors.
#
# Returns
# -------
# A named character vector where:
#
#   names(cluster_colors) = metadata values
#   values(cluster_colors) = hexadecimal color codes
#
# Example:
#
#   "0" = "#F8766D"
#   "1" = "#EC823C"
#   "2" = "#DD8D00"
# ------------------------------------------------------------
extract_dimplot_colors <- function(
    object,
    reduction = "umap.sketch",
    group.by = "seurat_cluster.sketched"
) {
  
  # Build a temporary DimPlot using Seurat's default
  # discrete color assignment.
  p <- DimPlot(
    object = object,
    reduction = reduction,
    group.by = group.by
  )
  
  # Build the ggplot object so the discrete color scale
  # is fully trained.
  gb <- ggplot_build(p)
  
  # Extract the trained color scale.
  color_scale <- gb$plot$scales$get_scales("colour")
  
  # Get the metadata levels in the order used by the
  # trained discrete scale.
  cluster_levels <- color_scale$get_limits()
  
  # Generate the colors assigned to those levels.
  cluster_colors <- color_scale$map(cluster_levels)
  
  # Create a named vector:
  #
  #   names  = cluster / metadata values
  #   values = corresponding DimPlot colors
  cluster_colors <- setNames(
    cluster_colors,
    cluster_levels
  )
  
  return(cluster_colors)
}