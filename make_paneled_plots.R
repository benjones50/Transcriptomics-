source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")

#makes panelled umap plots


# Creates standardized experimental metadata columns from the
# original tissue names. This includes:
#   - experimental_group
#   - tissue_id
#   - experimental_group_tissue
#   - genotype
#   - treatment
#   - timepoint
#   - phenotype
mega_obj <- add_experimental_metadata(mega_obj)



# ------------------------------------------------------------
# Plot dir
# ------------------------------------------------------------

#saving loading plot
plot_dir <- paste0(get_figures_dir(projection_dir),"/paneled_plots")

# ==========================================================
# Plot metadata highlights
#
# Responsibilities
# ----------------
# * Inspect available metadata columns and grouping variables
# * Verify metadata levels, classes, and sample counts
# * Generate UMAP highlight plots for key biological metadata
# * Save overview and individual highlight plots (if enabled)
#
# Notes
# -----
# * Each unique metadata value is highlighted against all other cells.
# * Metadata columns should be factors when a specific plotting order is
#   desired. Character columns will be plotted alphabetically.
# ==========================================================

levels(mega_obj$experimental_group)

class(mega_obj$experimental_group)

table(mega_obj$experimental_group)

names(mega_obj@meta.data)


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
cluster_colors <- extract_dimplot_colors( #pulls colors from umap
  object = mega_obj,
  reduction = "umap.sketch",
  group.by = "seurat_cluster.sketched"
)
# p_sketch <- DimPlot(
#   mega_obj,
#   reduction = "umap.sketch",
#   group.by = "seurat_cluster.sketched",
#   label = FALSE
# )
# p_sketch | p_sketch + scale_color_manual(
#   values = cluster_colors
# )
# 








#decides if saving each individual plot made, FALSE for testing
save_individuals = FALSE

exper_group_tissue_colored <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "experimental_group_tissue",
  color.by = "seurat_cluster.projected",
  colors = cluster_colors,
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = FALSE,
  show_cluster_labels = TRUE
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "slide_id",
  color.by = "seurat_cluster.projected",
  colors = cluster_colors,
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = FALSE,
  show_cluster_labels = TRUE
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "experimental_group",
  color.by = "seurat_cluster.projected",
  colors = cluster_colors,
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = FALSE,
  show_cluster_labels = TRUE
)



#not colored by the clusters


plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "experimental_group_tissue",
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = save_individuals
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "slide_id",
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = save_individuals,
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "experimental_group",
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = save_individuals
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "genotype",
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = save_individuals
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "treatment",
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = save_individuals
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "timepoint",
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = save_individuals
)

plot_out <- plot_metadata_highlights_2(
  object = mega_obj,
  group.by = "side",
  reduction = "full.umap.sketch",
  save_dir = plot_dir,
  save_individuals = save_individuals
)









