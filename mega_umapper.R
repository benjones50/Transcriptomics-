source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")



mega_obj <- load_mega_object(
  directory = mega_dir
)




cat("number of bins/cells in mega object:", ncol(mega_obj))



# pca_results <- run_pca_analysis(
#   object = mega_obj,
#   root_dir = mega_dir
# )
# 
# 
# mega_obj_pca <- pca_results$object
# 
# pca_results$loading_plot
# pca_results$elbow_plot
# pca_results$variance_plot



pca_dir <- get_pca_dir(
  root_dir = mega_dir,
  dims = 50
)


mega_obj <- load_or_build_object(
  directory = pca_dir,
  build_function = function() {
    run_sketch_pca(
      object = mega_obj,
      ncells = 50000
    )
  }
)



pca_plots <- plot_pca_diagnostics(
  object = mega_obj,
  pca_dir = pca_dir,
  reduction = "pca.sketch",
  dims = 1:6,
  ndims = 50
)

pca_plots$variance_plot
pca_plots$elbow_plot
pca_plots$loading_plot


#decide on max dims
max_dims <- 30

clustering_dir <- get_clustering_dir(
  root_dir = mega_dir,
  dims = 1:max_dims,
  resolution = 3
)


mega_obj <- load_or_build_object(
  directory = clustering_dir,
  build_function = function() {
    run_sketch_clustering(
      object = mega_obj,
      dims = 1:max_dims,
      resolution = 3
    )
  }
)


umap_dir <- get_umap_dir(
  root_dir = mega_dir,
  dims = 1:max_dims
)


mega_obj <- load_or_build_object(
  directory = umap_dir,
  build_function = function() {
    run_sketch_umap(
      object = mega_obj,
      dims = 1:max_dims
    )
  }
)


projection_dir <- get_projection_dir(
  root_dir = mega_dir,
  dims = 1:max_dims
)

mega_obj <- load_or_build_object(
  directory = projection_dir,
  build_function = function() {
    run_sketch_projection(
      object = mega_obj
    )
  }
)



# ------------------------------------------------------------
# 1) Basic object checks
# ------------------------------------------------------------
cat("===== Assays =====\n")
print(Assays(mega_obj_umapped))

cat("\n===== Reductions =====\n")
print(Reductions(mega_obj_umapped))

cat("\n===== Metadata cluster columns =====\n")
print(grep("cluster", colnames(mega_obj_umapped[[]]), value = TRUE))


# ------------------------------------------------------------
# 2) Sketch size check
# ------------------------------------------------------------
cat("Original cells:", ncol(mega_obj_umapped), "\n")
cat("Sketch cells:", ncol(mega_obj_umapped[["sketch"]]), "\n")
cat(
  sprintf(
    "Sketch represents %.2f%% of all cells\n",
    100 * ncol(mega_obj_umapped[["sketch"]]) / ncol(mega_obj_umapped)
  )
)


# ------------------------------------------------------------
# 3) Projection checks, makes sure every cell projected to
# ------------------------------------------------------------
# proj <- Embeddings(mega_obj_umapped, "projected.pca")
# cat("Projected cells:", nrow(proj), "\n")
# cat("Total cells:", ncol(mega_obj_umapped), "\n")
# stopifnot(nrow(proj) == ncol(mega_obj_umapped))
# 
# umap <- Embeddings(mega_obj_umapped, "full.umap.sketch")
# stopifnot(nrow(umap) == ncol(mega_obj_umapped))



# ------------------------------------------------------------
# 4) Cluster assignment checks
# ------------------------------------------------------------

cat("Projected cluster assignments:\n")

cat(
  "Cells with projected clusters:",
  sum(!is.na(mega_obj_umapped$seurat_cluster.projected)),
  "\n"
)

cat(
  "Cells missing projected clusters:",
  sum(is.na(mega_obj_umapped$seurat_cluster.projected)),
  "\n"
)

cat(
  "\nNumber of projected clusters:",
  length(unique(mega_obj_umapped$seurat_cluster.projected)),
  "\n"
)

# ------------------------------------------------------------
# Prefix for all diagnostic plots #TODO rethink this, seems bad
# ------------------------------------------------------------
plot_prefix <- paste0(
  "mega_",
  analysis_mode,
  if (!is.null(bin_size)) paste0("_", bin_size, "um") else "",
  "_",
  normalization
)


# ------------------------------------------------------------
# Plot dir
# ------------------------------------------------------------

#saving loading plot
plot_dir <- get_figures_dir(projection_dir)


# ------------------------------------------------------------
# 4b) Projection confidence
#
# Visualize the confidence score assigned by ProjectData() for
# each projected cluster assignment. Higher scores indicate
# greater confidence that a cell belongs to its assigned cluster.
# ------------------------------------------------------------

summary(
  mega_obj_umapped$seurat_cluster.projected.score
)

# ------------------------------------------------------------
# Order clusters by abundance (largest -> smallest)
# ------------------------------------------------------------

cluster_counts <- table(mega_obj_umapped$seurat_cluster.projected)

cluster_percent <- 100 * cluster_counts / sum(cluster_counts)

cluster_order <- names(
  sort(
    cluster_counts,
    decreasing = TRUE
  )
)

# reorders based on cell abundance order
mega_obj_umapped$seurat_cluster.projected <- factor(
  mega_obj_umapped$seurat_cluster.projected,
  levels = cluster_order
)
# Tell Seurat to group violins by the projected clusters
Idents(mega_obj_umapped) <- "seurat_cluster.projected"

# ------------------------------------------------------------
# Create informative x-axis labels
#
# Each label shows:
#   Cluster ID
#   Percent of total cells
#   Number of cells
# ------------------------------------------------------------

cluster_label_map <- setNames(
  paste0(
    names(cluster_counts),
    "\n",
    sprintf("%.1f%%", cluster_percent),
    "\n",
    format(cluster_counts, big.mark = ",")
  ),
  names(cluster_counts)
)

# ------------------------------------------------------------
# Projection confidence violin plot
# ------------------------------------------------------------

cluster_proj_score_plot <- VlnPlot(
  mega_obj_umapped,
  features = "seurat_cluster.projected.score",
  pt.size = 0
) +
  scale_x_discrete(
    labels = cluster_label_map
  ) +
  xlab("Projected Cluster\n(% of Total Cells / Cell Count)") +
  ggtitle("Projected Cluster Assignment Confidence") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    )
  )

save_plot(
  plot_object = cluster_proj_score_plot,
  file_stub = paste0(plot_prefix, "_seurat_cluster.projected.score"),
  save_dir = plot_dir,
  width = 48,
  height = 12
)

# allows us to see:
# clusters with consistently high-confidence assignments,
# clusters with low-confidence assignments,
# clusters that may represent transitions or ambiguous cell states.

#TODO: update numbers, removed #5

#undos ordering clusters based on size
mega_obj_umapped$seurat_cluster.projected <- factor(
  mega_obj_umapped$seurat_cluster.projected,
  levels = sort(
    unique(as.numeric(as.character(
      mega_obj_umapped$seurat_cluster.projected
    )))
  )
)

Idents(mega_obj_umapped) <- "seurat_cluster.projected"
# ------------------------------------------------------------
# 6) Plot: slide mixing on projected UMAP
# ------------------------------------------------------------
p_slide_mix <- DimPlot(
  mega_obj_umapped,
  reduction = "full.umap.sketch",
  group.by = "sample_tissue"
) +
  ggtitle("Projected UMAP Colored by Original Visium HD Slide")

save_plot(
  plot_object = p_slide_mix,
  file_stub = paste0(plot_prefix, "_umap_slide_mixing"),
  save_dir = plot_dir,
  width = 9,
  height = 7
)


p_tissue <- DimPlot(
  mega_obj_umapped,
  reduction = "full.umap.sketch",
  group.by = "tissue"
) +
  ggtitle("Projected UMAP Colored by Tissue")

save_plot(
  plot_object = p_tissue,
  file_stub = paste0(plot_prefix, "_umap_per_tissue_mixing"),
  save_dir = plot_dir,
  width = 9,           
  height = 7,
)


# ------------------------------------------------------------
# 7) Plot: sketch cell coverage
# ------------------------------------------------------------
p_sketch_cells <- DimPlot(
  mega_obj_umapped,
  reduction = "full.umap.sketch",
  cells.highlight = colnames(mega_obj_umapped[["sketch"]]),
  cols.highlight = "red",
  cols = "blue"
) +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c(
      "Bins not selected for sketch",
      "Bins selected by SketchData()"
    ),
    name = NULL
  ) +
  ggtitle("Projected UMAP Highlighting Sketch Bins")

save_plot(
  plot_object = p_sketch_cells,
  file_stub = paste0(plot_prefix, "_umap_sketch_cells_highlighted"),
  save_dir = plot_dir,
  width = 9,           
  height = 7,
)



# Save current assay so we can restore it afterward
original_assay <- DefaultAssay(mega_obj_umapped)


# ------------------------------------------------------------
# Information for plot titles
# ------------------------------------------------------------

n_sketch_cells <- ncol(mega_obj_umapped[["sketch"]])
n_total_cells  <- ncol(mega_obj_umapped)

sketch_title <- sprintf(
  "Sketch UMAP\nRepresentative sketch of %s cells (%.1f%% of dataset)",
  format(n_sketch_cells, big.mark = ","),
  100 * n_sketch_cells / n_total_cells
)

projected_title <- sprintf(
  "Projected UMAP\nProjection of all %s cells into sketch UMAP",
  format(n_total_cells, big.mark = ",")
)

# ------------------------------------------------------------
# Sketch UMAP
# ------------------------------------------------------------
DefaultAssay(mega_obj_umapped) <- "sketch"
Idents(mega_obj_umapped) <- "seurat_cluster.sketched"

p_sketch <- DimPlot(
  mega_obj_umapped,
  reduction = "umap.sketch",
  label = FALSE
) +
  ggtitle(sketch_title) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# Projected UMAP
# ------------------------------------------------------------
DefaultAssay(mega_obj_umapped) <- original_assay
Idents(mega_obj_umapped) <- "seurat_cluster.projected"

p_projected <- DimPlot(
  mega_obj_umapped,
  reduction = "full.umap.sketch",
  label = FALSE
) +
  ggtitle(projected_title) +
  theme(legend.position = "bottom")

# Restore original assay
DefaultAssay(mega_obj_umapped) <- original_assay

# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------
p_projection_compare <- p_sketch | p_projected

save_plot(
  plot_object = p_projection_compare,
  file_stub = paste0(plot_prefix, "_umap_projection_comparison"),
  save_dir = plot_dir,
  width = 14,
  height = 10
)




# quantitative view of where the sketch differs from simple proportional sampling.
#Was this neighborhood sampled more or less than expected by SketchData()?
# ------------------------------------------------------------
# Local Sketch Sampling Deviation
#
# Colors each region of the projected UMAP by how much the local
# sketch sampling differs from the expected global sampling rate.
#
# White  = sampled as expected
# Red    = over-sampled
# Blue   = under-sampled
# ------------------------------------------------------------

# ------------------------------------------------------------
# Extract projected UMAP coordinates
# ------------------------------------------------------------

umap_df <- as.data.frame(
  Embeddings(
    mega_obj_umapped,
    "full.umap.sketch"
  )
)

colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)

# Identify sketch cells
umap_df$sketch <- umap_df$cell %in%
  colnames(mega_obj_umapped[["sketch"]])

# ------------------------------------------------------------
# Expected sampling ratio
# ------------------------------------------------------------

expected_ratio <-
  ncol(mega_obj_umapped[["sketch"]]) /
  ncol(mega_obj_umapped)

cat(
  sprintf(
    "Expected sketch sampling ratio: %.2f%%\n",
    100 * expected_ratio
  )
)

# ------------------------------------------------------------
# Bin UMAP into a regular grid
# ------------------------------------------------------------

grid_size <- 100

umap_df <- umap_df %>%
  mutate(
    x_bin = cut(
      UMAP_1,
      breaks = grid_size,
      labels = FALSE
    ),
    y_bin = cut(
      UMAP_2,
      breaks = grid_size,
      labels = FALSE
    )
  )

# ------------------------------------------------------------
# Compute local sampling statistics
# ------------------------------------------------------------

ratio_df <- umap_df %>%
  group_by(
    x_bin,
    y_bin
  ) %>%
  summarise(
    total_cells = n(),
    sketch_cells = sum(sketch),
    sampling_ratio = sketch_cells / total_cells,
    sampling_difference = sampling_ratio - expected_ratio,
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Plot local sampling deviation
# ------------------------------------------------------------

sampling_ratio_plot <- ggplot(
  ratio_df,
  aes(
    UMAP_1,
    UMAP_2,
    color = sampling_difference
  )
) +
  geom_point(
    size = 1.5
  ) +
  scale_color_gradient2(
    low = "royalblue3",
    mid = "white",
    high = "firebrick3",
    midpoint = 0,
    labels = scales::percent,
    name = "Deviation\nfrom Expected"
  ) +
  coord_equal() +
  labs(
    title = "Local Sketch Sampling Deviation",
    subtitle = sprintf(
      "Expected sampling = %.1f%% of all cells",
      100 * expected_ratio
    ),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic()


# ------------------------------------------------------------
# Projected clusters
# ------------------------------------------------------------

p_projected <- DimPlot(
  mega_obj_umapped,
  reduction = "full.umap.sketch",
  group.by = "seurat_cluster.projected",
  label = FALSE
) +
  ggtitle("Projected UMAP Colored by Projected Cluster")

# ------------------------------------------------------------
# Side-by-side comparison
# ------------------------------------------------------------

sampling_qc_plot <- p_projected | sampling_ratio_plot

save_plot(
  plot_object = sampling_qc_plot,
  file_stub = paste0(
    plot_prefix,
    "_sampling_ratio_qc"
  ),
  save_dir = plot_dir,
  width = 16,
  height = 7
)








#additional plots via mega_umapper_plotter
plot_umap_features(
  object = mega_obj_umapped,
  features = c(
    paste0("nCount_", DefaultAssay(mega_obj_umapped)),
    paste0("nFeature_", DefaultAssay(mega_obj_umapped)),
    "percent.mt"
  ),
  save_dir = file.path(
    plot_dir,
    "QC"
  )
)

plot_umap_metadata(
  object = mega_obj_umapped,
  variables = c(
    "sample_tissue",
    "tissue"
  ),
  save_dir = file.path(
    plot_dir,
    "Metadata"
  )
)

# plot_umap_features(
#   object = mega_obj_umapped,
#   features = c(
#     "EPCAM",
#     "KRT8",
#     "KRT19",
#     "MUC1"
#   ),
#   save_dir = file.path(
#     plot_dir,
#     "Markers"
#   )
# )
# 
# options(future.globals.maxSize = 16 * 1024^3)

