#full credit to https://romanhaa.github.io/projects/scrnaseq_workflow/

# ------------------------------------------------------------
# Build and plot cluster tree
#
# Construct a hierarchical tree of clusters using a specified
# dimensional reduction and generate a dendrogram showing the
# relationships between clusters.
#
# Returns:
#   - plot: ggtree dendrogram
# ------------------------------------------------------------

plot_cluster_tree <- function(
    object,
    reduction,
    ndims,
    cluster_colors = NULL
) {
  
  # Build hierarchical cluster tree
  object <- BuildClusterTree(
    object = object,
    reduction = reduction,
    dims = 1:ndims,
    reorder = FALSE,
    reorder.numeric = FALSE
  )
  
  # Extract phylogenetic tree
  tree <- object@tools$BuildClusterTree
  
  # Save original cluster IDs
  tree$tip.cluster <- tree$tip.label
  
  # Pretty labels
  tree$tip.label <- paste0("Cluster ", tree$tip.label)
  
  p <- ggtree::ggtree(tree) +
    scale_y_reverse() +
    ggtree::geom_tree()
  
  if (!is.null(cluster_colors)) {
    
    # Match colors to the tree ordering
    tip_df <- p$data[p$data$isTip, ]
    tip_df$cluster <- tree$tip.cluster
    
    p <- p +
      ggtree::geom_tippoint(
        data = tip_df,
        aes(color = cluster),
        size = 4
      ) +
      ggtree::geom_tiplab(
        data = tip_df,
        aes(color = cluster),
        offset = 1
      ) +
      scale_color_manual(values = cluster_colors)
    
  } else {
    
    p <- p +
      ggtree::geom_tippoint(size = 4) +
      ggtree::geom_tiplab(offset = 1)
    
  }
  
  p +
    ggtree::theme_tree() +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = unit(c(0, 2.5, 0, 0), "cm"),
      legend.position = "none"
    )
}