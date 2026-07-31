source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")
#this is where ive been running my pca's, harmony, clustering, umapping, and projecting. then lots of umap plots

mega_obj <- load_mega_object(
  directory = mega_dir
)
cat("number of bins/cells in mega object:", ncol(mega_obj))


options(future.globals.maxSize = 1000 * 1024^2)



sketch_label <- "sketch_33pct" #name of sketch
ncells_total <- ncol(mega_obj)
ncells_sketch <- ncells_total%/%3 #defining how sketch actually made

pca_dir <- get_pca_dir(
  root_dir = mega_dir,
  ndims = 50,
  sketch_label = sketch_label
)


mega_obj <- load_or_build_object(
  directory = pca_dir,
  build_function = function() {
    run_sketch_pca(
      object = mega_obj,
      ncells = ncells_sketch
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


# Decide on embedding
max_dims <- 11



if (embedding == "harmony") {
  group.by.vars <- "tissue"
  
  
  harmony_dir <- get_harmony_dir(
    root_dir = mega_dir,
    ndims = max_dims,
    sketch_label = "sketch_33pct",
    group.by.vars = group.by.vars
  )
  
  mega_obj <- load_or_build_object(
    directory = harmony_dir,
    build_function = function() {
      run_harmony_embedding(
        object = mega_obj,
        ndims = max_dims,
        sketch_reduction = "pca.sketch",
        group.by.vars = "tissue"
      )
    },
    force_rebuild = FALSE
  )
  
  embedding_run_dir <- harmony_dir
  
} else {
  
  embedding_run_dir <- pca_dir
  
}



resolution <- 0.3


#alternative directory used for testing
# clustering_dir_normal <- get_clustering_dir(
#   embedding_run_dir = embedding_run_dir,
#   ndims = max_dims,
#   resolution = resolution
# )
# 
# clustering_dir <- get_nested_subdir(
#   clustering_dir_normal,
#   "multi_res_0.5"
# )

#cluster / load clustered object
clustering_dir <- get_clustering_dir(
  embedding_run_dir = embedding_run_dir,
  ndims = max_dims,
  resolution = resolution
)

mega_obj <- load_or_build_object(
  directory = clustering_dir,
  build_function = function() {
    
    object <- run_sketch_clustering(
      object = mega_obj,
      sketch_reduction = embedding_reduction,
      dims = 1:max_dims,
      resolution = resolution
    )
    
    object <- add_silhouette( #approx silhouette using centroids of clusters
      object = object,
      reduction = embedding_reduction, 
      dims = 1:max_dims,
      cluster_column = "seurat_cluster.sketched"
    )
    object
  },
  force_rebuild = TRUE
)
 


umap_dir <- get_umap_dir(
  clustering_dir = clustering_dir
)



# not used bc don't need to save umap
# mega_obj <- load_or_build_object(
#   directory = umap_dir,
#   build_function = function() {
#     run_sketch_umap(
#       object = mega_obj,
#       sketch_reduction = embedding_reduction,
#       dims = 1:max_dims
#     )
#   },
#   force_rebuild = TRUE
# )
mega_obj <- run_sketch_umap(
      object = mega_obj,
      sketch_reduction = embedding_reduction,
      dims = 1:max_dims
    )

projection_dir <- get_projection_dir(
  umap_dir = umap_dir
)


mega_obj <- load_or_build_object(
  directory = projection_dir,
  build_function = function() {
    run_sketch_projection(
      object = mega_obj,
      normalization_method = normalization_method,
      sketch_reduction = embedding_reduction,
      ndims = max_dims
    )
  },
  force_rebuild = FALSE
)






mega_obj <- add_silhouette( #approx silhouette using centroids of clusters
      object = mega_obj,
      reduction = embedding_reduction, 
      dims = 1:max_dims,
      cluster_column = "seurat_cluster.sketched" # "seurat_cluster.harmony" "seurat_cluster.sketched"?
    )




# ------------------------------------------------------------
# Basic object checks
# ------------------------------------------------------------
cat("===== Assays =====\n")
print(Assays(mega_obj))

cat("\n===== Reductions =====\n")
print(Reductions(mega_obj))

cat("\n===== Metadata cluster columns =====\n")
print(grep("cluster", colnames(mega_obj[[]]), value = TRUE))


# ------------------------------------------------------------
# Sketch size check
# ------------------------------------------------------------
cat("Original cells:", ncol(mega_obj), "\n")
cat("Sketch cells:", ncol(mega_obj[["sketch"]]), "\n")
cat(
  sprintf(
    "Sketch represents %.2f%% of all cells\n",
    100 * ncol(mega_obj[["sketch"]]) / ncol(mega_obj)
  )
)


# ------------------------------------------------------------
# Cluster assignment checks
# ------------------------------------------------------------

cat("Projected cluster assignments:\n")

cat(
  "Cells with projected clusters:",
  sum(!is.na(mega_obj$seurat_cluster.projected)),
  "\n"
)

cat(
  "Cells missing projected clusters:",
  sum(is.na(mega_obj$seurat_cluster.projected)),
  "\n"
)

sum(is.na(mega_obj$seurat_cluster.projected.score))
# 4092

cat(
  "\nNumber of projected clusters:",
  length(unique(mega_obj$seurat_cluster.projected)),
  "\n"
)



# ------------------------------------------------------------
# Plot dir
# ------------------------------------------------------------

plot_dir <- paste0(get_figures_dir(projection_dir),"/mega_umapper")

# ------------------------------------------------------------
# Information for plot titles
# ------------------------------------------------------------

n_sketch_cells <- ncol(mega_obj[["sketch"]])
n_total_cells  <- ncol(mega_obj)

sketch_title <- sprintf(
  "Sketch UMAP\nRepresentative sketch of %s cells (%.1f%% of dataset)",
  format(n_sketch_cells, big.mark = ","),
  100 * n_sketch_cells / n_total_cells
)

projected_title <- sprintf(
  "Projected UMAP\nProjection of all %s cells into sketch UMAP",
  format(n_total_cells, big.mark = ",")
)


#based on code by https://romanhaa.github.io/projects/scrnaseq_workflow/

# ------------------------------------------------------------
# Projected UMAP
# ------------------------------------------------------------

# Convert projected clusters to the same factor ordering as the
# sketch clusters so cluster colors remain consistent across plots.
mega_obj$seurat_cluster.projected <- factor(
  mega_obj$seurat_cluster.projected,
  levels = levels(mega_obj$seurat_cluster.sketched)
)

# ------------------------------------------------------------
# Compute cluster label positions
#
# Labels are placed at the median UMAP coordinate of each cluster,
# which is more robust to outliers than using the mean.
# ------------------------------------------------------------

# Extract projected UMAP coordinates
projected_umap_df <- as.data.frame(
  Embeddings(
    mega_obj,
    "full.umap.sketch"
  )
)

# Give the coordinate columns consistent names
colnames(projected_umap_df)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

# Add projected cluster assignments
projected_umap_df$cluster <- mega_obj$seurat_cluster.projected

# Compute median UMAP position for each cluster
projected_umap_centers <- projected_umap_df %>%
  dplyr::filter(
    !is.na(cluster)
  ) %>%
  dplyr::group_by(
    cluster
  ) %>%
  dplyr::summarise(
    x = median(UMAP_1),
    y = median(UMAP_2),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Build projected UMAP
# ------------------------------------------------------------

p_projected <- DimPlot(
  mega_obj,
  reduction = "full.umap.sketch",
  group.by = "seurat_cluster.projected",
  label = FALSE
) +
  ggtitle(projected_title) +
  theme(
    legend.position = "bottom"
  ) +
  
# --------------------------------------------------------
# Add cluster labels
# --------------------------------------------------------
geom_label(
  data = projected_umap_centers,
  mapping = aes(
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
) +
  
# --------------------------------------------------------
# Annotate the total number of cells/bins
# --------------------------------------------------------
annotate(
  geom = "text",
  x = Inf,
  y = -Inf,
  label = paste0(
    "n = ",
    format(
      ncol(mega_obj),
      big.mark = ",",
      trim = TRUE
    )
  ),
  vjust = -1.5,
  hjust = 1.25,
  color = "black",
  size = 2.5
)

save_plot(
  plot_object = p_projected,
  file_stub = paste0("projection_labeled"),
  save_dir = plot_dir,
  width = 7,
  height = 7
)

# ------------------------------------------------------------
# Sketch UMAP
# ------------------------------------------------------------

# ------------------------------------------------------------
# Compute cluster label positions
#
# Labels are placed at the median UMAP coordinate of each
# sketch cluster, which is more robust to outliers than the mean.
# ------------------------------------------------------------

# Extract sketch UMAP coordinates
sketch_umap_df <- as.data.frame(
  Embeddings(
    mega_obj,
    "umap.sketch"
  )
)

# Give the coordinate columns consistent names
colnames(sketch_umap_df)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

# Add sketch cluster assignments for the cells in the embedding
sketch_umap_df$cluster <-
  mega_obj[[]][
    rownames(sketch_umap_df),
    "seurat_cluster.sketched"
  ]

# Compute median UMAP position for each cluster
sketch_umap_centers <- sketch_umap_df %>%
  dplyr::filter(
    !is.na(cluster)
  ) %>%
  dplyr::group_by(
    cluster
  ) %>%
  dplyr::summarise(
    x = median(UMAP_1),
    y = median(UMAP_2),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Build sketch UMAP
# ------------------------------------------------------------

p_sketch <- DimPlot(
  mega_obj,
  reduction = "umap.sketch",
  group.by = "seurat_cluster.sketched",
  label = FALSE
) +
  ggtitle(sketch_title) +
  theme(
    legend.position = "bottom"
  ) +
  
# ----------------------------------------------------------
# Add cluster labels
# ----------------------------------------------------------
geom_label(
  data = sketch_umap_centers,
  mapping = aes(
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
) +
  
# ----------------------------------------------------------
# Annotate the number of sketch cells
# ----------------------------------------------------------
annotate(
  geom = "text",
  x = Inf,
  y = -Inf,
  label = paste0(
    "n = ",
    format(
      nrow(sketch_umap_df),
      big.mark = ",",
      trim = TRUE
    )
  ),
  vjust = -1.5,
  hjust = 1.25,
  color = "black",
  size = 2.5
)

save_plot(
  plot_object = p_sketch,
  file_stub = "sketch_labeled",
  save_dir = plot_dir,
  width = 7,
  height = 7
)

p_projection_compare <- p_sketch | p_projected

save_plot(
  plot_object = p_projection_compare,
  file_stub = paste0("umap_projection_comparison"),
  save_dir = plot_dir,
  width = 14,
  height = 10
)





# ------------------------------------------------------------
# Projection confidence and sketch silhouette comparison
#
# ------------------------------------------------------------

# Summarize projected cluster assignment confidence
cat("\nSummary of Clustering Confidence:\n")
summary(mega_obj$seurat_cluster.projected.score)

# ------------------------------------------------------------
# Build a shared cluster ordering from the sketch clusters
#
# We order clusters by the number of sketch cells in each group
# ------------------------------------------------------------

cluster_counts <- sort(
  table(mega_obj$seurat_cluster.sketched),
  decreasing = TRUE
)

cluster_order <- names(cluster_counts)

cluster_percent <- 100 * cluster_counts / sum(cluster_counts)

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


mega_obj$seurat_cluster.projected.ordered <- factor(
  mega_obj$seurat_cluster.projected,
  levels = cluster_order
)

mega_obj$seurat_cluster.sketched.ordered <- factor(
  mega_obj$seurat_cluster.sketched,
  levels = cluster_order
)

# ------------------------------------------------------------
# Projection confidence violin plot
#
# This shows the confidence score assigned by ProjectData()
# for each projected cluster assignment.
# ------------------------------------------------------------

cluster_proj_score_plot <- VlnPlot(
  mega_obj,
  features = "seurat_cluster.projected.score",
  group.by = "seurat_cluster.projected.ordered",
  pt.size = 0
) +
  scale_x_discrete(labels = cluster_label_map) +
  xlab("Projected Cluster\n(% of Sketch Cells / Cell Count)") +
  ggtitle("Projected Cluster Assignment Confidence") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    )
  )

# ------------------------------------------------------------
# Prepare a sketch-only object for the silhouette plot
#
# Silhouette values exist only for sketch cells, so subset the
# object to cells with non-NA silhouette values.
# ------------------------------------------------------------

sketch_cells <- WhichCells(
  mega_obj,
  expression = !is.na(silhouette)
)

sketch_obj <- subset(
  mega_obj,
  cells = sketch_cells
)

# Add the ordered sketch grouping column to the subset object
sketch_obj$seurat_cluster.sketched.ordered <- factor(
  sketch_obj$seurat_cluster.sketched,
  levels = cluster_order
)

# ------------------------------------------------------------
# Approximate sketch silhouette violin plot
#
# This shows silhouette values for the sketch cells using the
# same cluster order as the projection confidence plot.
# ------------------------------------------------------------

sketch_silhouette_plot <- VlnPlot(
  sketch_obj,
  features = "silhouette",
  group.by = "seurat_cluster.sketched.ordered",
  pt.size = 0
) +
  scale_x_discrete(labels = cluster_label_map) +
  xlab("Sketch Cluster\n(% of Sketch Cells / Cell Count)") +
  ggtitle("Approximate Sketch Silhouette") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    )
  )

# ------------------------------------------------------------
# Combine both plots
# ------------------------------------------------------------
combined_qc_plot <-
  cluster_proj_score_plot /
  sketch_silhouette_plot

combined_qc_plot



save_plot(
  plot_object = cluster_proj_score_plot,
  file_stub = "seurat_cluster.projected.score",
  save_dir = plot_dir,
  width = 12,
  height = 6
)

save_plot(
  plot_object = sketch_silhouette_plot,
  file_stub = "sketch_silhouette",
  save_dir = plot_dir,
  width = 12,
  height = 6
)

save_plot(
  plot_object = combined_qc_plot,
  file_stub = "cluster_quality",
  save_dir = plot_dir,
  width = 12,
  height = 12
)








# ------------------------------------------------------------
#
# Umaps colored by silhouette scores
#
# ------------------------------------------------------------
umap_silhouette_plot_rwb <- FeaturePlot(
  sketch_obj,
  features = "silhouette",
  reduction = "umap.sketch",
  pt.size = 0.25,
  order = TRUE,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_color_gradient2(
    low = "#313695",
    mid = "white",
    high = "#A50026",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Silhouette"
  ) +
  ggtitle("Approximate Sketch Silhouette") +
  theme_classic()

save_plot(
  plot_object = umap_silhouette_plot_rwb,
  file_stub = "umap_silhouette_plot_rwb",
  save_dir = plot_dir,
  width = 12,
  height = 10,
  dpi = 1200
)

umap_silhouette_plot_red_green <- FeaturePlot(
  sketch_obj,
  features = "silhouette",
  reduction = "umap.sketch",
  pt.size = 0.25,
  order = TRUE,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_color_gradientn(
    colors = c(
      "#d73027",  # red
      "#fc8d59",  # orange
      "#fee08b",  # yellow
      "#91cf60",  # light green
      "#1a9850"   # dark green
    ),
    values = scales::rescale(c(-1, 0, 0.25, 0.50, 0.75)),
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Silhouette"
  ) +
  ggtitle("Approximate Sketch Silhouette") +
  theme_classic()

save_plot(
  plot_object = umap_silhouette_plot_red_green,
  file_stub = "umap_silhouette_plot_red_green",
  save_dir = plot_dir,
  width = 12,
  height = 10,
  dpi = 1200
)

#----------------------------------------
# Categorize silhouette scores
#----------------------------------------

sketch_obj$silhouette_category <- cut(
  sketch_obj$silhouette,
  breaks = c(-Inf, 0, 0.25, 0.50, 0.75, Inf),
  labels = c(
    "< 0",
    "0 - 0.25",
    "0.25 - 0.50",
    "0.50 - 0.75",
    "> 0.75"
  ),
  include.lowest = TRUE,
  right = FALSE
)

#----------------------------------------
# UMAP colored by silhouette category
#----------------------------------------

silhouette_umap_categorical <-
  DimPlot(
    sketch_obj,
    reduction = "umap.sketch",
    group.by = "silhouette_category",
    pt.size = 0.01,
    shuffle = TRUE,
    raster = FALSE
  ) +
  scale_color_manual(
    values = c(
      "< 0"         = "#d73027",
      "0 - 0.25"    = "#fc8d59",
      "0.25 - 0.50" = "#fee08b",
      "0.50 - 0.75" = "#91cf60",
      "> 0.75"      = "#1a9850"
    ),
    drop = FALSE,
    name = "Silhouette"
  ) +
  ggtitle("Approximate Sketch Silhouette") +
  theme_classic()


save_plot(
  plot_object = silhouette_umap_categorical,
  file_stub = paste0("silhouette_umap_categorical"),
  save_dir = plot_dir,
  width = 12,
  height = 10,
  dpi = 1200
)

#combine and save
silhouette_umaps <- umap_silhouette_plot_red_green |
  silhouette_umap_categorical

save_plot(
  plot_object = silhouette_umaps,
  file_stub = paste0("combined_silhouette_umaps"),
  save_dir = plot_dir,
  width = 20,
  height = 8,
  dpi = 1200
)



rm(sketch_obj) #cleanup memory, removes obj with only sketch cells









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
# slide mixing on projected UMAP
# ------------------------------------------------------------
p_slide_mix <- DimPlot(
  mega_obj,
  reduction = "full.umap.sketch",
  group.by = "slide_id"
) +
  ggtitle("Projected UMAP Colored by Original Visium HD Slide")

save_plot(
  plot_object = p_slide_mix,
  file_stub = "umap_slide_mixing",
  save_dir = plot_dir,
  width = 9,
  height = 7
)

p_tissue <- DimPlot(
  mega_obj,
  reduction = "full.umap.sketch",
  group.by = "tissue"
) +
  ggtitle("Projected UMAP Colored by Tissue")

save_plot(
  plot_object = p_tissue,
  file_stub = "experimental_group",
  save_dir = plot_dir,
  width = 9,           
  height = 7,
)


# ------------------------------------------------------------
# Projection uncertainty vs tissue identity
#
# Visualize where ProjectData() was most and least confident
# assigning projected cluster labels.
#
# plots help determine whether uncertainty is
# driven by particular tissues or shared across the dataset.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Compute projection uncertainty
# ------------------------------------------------------------

#default version

p_confidence <- FeaturePlot(
  mega_obj,
  reduction = "full.umap.sketch",
  features = "seurat_cluster.projected.score",
  order = TRUE,
  raster = FALSE,
  pt.size = 0.0001,
) +
  ggtitle("Projected Cluster Assignment Confidence")

confidence_vs_tissue_plot <-
  p_confidence |
  p_tissue

save_plot(
  plot_object = confidence_vs_tissue_plot,
  file_stub = "umap_confidence_vs_tissue",
  save_dir = plot_dir,
  width = 15,
  height = 7
)



mega_obj$projection_uncertainty <-
  1 - mega_obj$seurat_cluster.projected.score

# ------------------------------------------------------------
# Projection uncertainty UMAP
# ------------------------------------------------------------

p_projection_uncertainty <- FeaturePlot(
  mega_obj,
  reduction = "full.umap.sketch",
  features = "projection_uncertainty",
  order = TRUE,
  raster = FALSE,
  pt.size = 0.001
) +
  scale_color_gradientn(
    colors = c(
      "#440154",
      "#31688E",
      "#35B779",
      "#FDE725"
    ),
    limits = c(0, 1),
    name = "Projection\nUncertainty"
  ) +
  ggtitle("Projected UMAP Colored by Projection Uncertainty")

# ------------------------------------------------------------
# Tissue identity UMAP
# ------------------------------------------------------------

p_tissue_identity <- DimPlot(
  mega_obj,
  reduction = "full.umap.sketch",
  group.by = "experimental_group"
) +
  ggtitle("Projected UMAP Colored by Tissue")

# ------------------------------------------------------------
# Side-by-side comparison
# ------------------------------------------------------------

projection_uncertainty_plot <-
  p_projection_uncertainty |
  p_tissue_identity

save_plot(
  plot_object = projection_uncertainty_plot,
  file_stub = "umap_projection_uncertainty_vs_tissue",
  save_dir = plot_dir,
  width = 18,
  height = 7
)

# ------------------------------------------------------------
# 7) Plot: sketch cell coverage
# ------------------------------------------------------------
p_sketch_cells <- DimPlot(
  mega_obj,
  reduction = "full.umap.sketch",
  cells.highlight = colnames(mega_obj[["sketch"]]),
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
  file_stub = "umap_sketch_cells_highlighted",
  save_dir = plot_dir,
  width = 9,           
  height = 7,
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
    mega_obj,
    "full.umap.sketch"
  )
)

colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)

# Identify sketch cells
umap_df$sketch <- umap_df$cell %in%
  colnames(mega_obj[["sketch"]])

# ------------------------------------------------------------
# Expected sampling ratio
# ------------------------------------------------------------

expected_ratio <-
  ncol(mega_obj[["sketch"]]) /
  ncol(mega_obj)

cat(
  sprintf(
    "Expected sketch sampling ratio: %.2f%%\n",
    100 * expected_ratio
  )
)

# ------------------------------------------------------------
# Bin UMAP into a grid
# ------------------------------------------------------------

grid_size <- 1000

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
    size = 0.05
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
# Side-by-side comparison
# ------------------------------------------------------------

sampling_qc_plot <- p_projected | sampling_ratio_plot

save_plot(
  plot_object = sampling_qc_plot,
  file_stub ="sampling_ratio_qc",
  save_dir = plot_dir,
  width = 16,
  height = 7
)






# ------------------------------------------------------------
# Local UMAP Density
#
# Colors each region of the UMAP by the number of cells within
# a square neighborhood. Smaller bin_width values produce
# higher-resolution density maps.
# ------------------------------------------------------------
# Extract UMAP coordinates
# ------------------------------------------------------------

# ------------------------------------------------------------
# Extract UMAP coordinates
# ------------------------------------------------------------

umap_df <- as.data.frame(
  Embeddings(
    mega_obj,
    "full.umap.sketch"
  )
)

colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# ------------------------------------------------------------
# Bin UMAP into a regular grid
#
# Smaller bin_width = more bins = higher resolution
# Larger bin_width  = fewer bins = smoother density
# ------------------------------------------------------------

bin_width <- 0.1

umap_df <- umap_df %>%
  mutate(
    x_bin = floor(UMAP_1 / bin_width),
    y_bin = floor(UMAP_2 / bin_width)
  )

# ------------------------------------------------------------
# Count cells within each grid square
# ------------------------------------------------------------

density_df <- umap_df %>%
  group_by(
    x_bin,
    y_bin
  ) %>%
  summarise(
    n_cells = n(),
    UMAP_1 = (first(x_bin) + 0.5) * bin_width,
    UMAP_2 = (first(y_bin) + 0.5) * bin_width,
    .groups = "drop"
  )

# ------------------------------------------------------------
# Clip color scale to improve contrast
#
# Change this percentile if desired:
#   0.99  = strongest contrast
#   0.995 = moderate
#   0.999 = conservative
# ------------------------------------------------------------

clip_percentile <- 0.995

upper_limit <- quantile(
  density_df$n_cells,
  clip_percentile
)

# ------------------------------------------------------------
# Plot local density
# ------------------------------------------------------------

density_plot <- ggplot(
  density_df,
  aes(
    x = UMAP_1,
    y = UMAP_2
  )
) +
  geom_tile(
    aes(
      fill = pmin(n_cells, upper_limit)
    ),
    width = bin_width,
    height = bin_width
  ) +
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(0, upper_limit),
    oob = scales::squish,
    name = "Cells\nper bin"
  ) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Local UMAP Density",
    subtitle = sprintf(
      "Bin width = %.2f UMAP units | Color scale clipped at %.1fth percentile (%.0f cells/bin)",
      bin_width,
      clip_percentile * 100,
      upper_limit
    ),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic()

density_plot_combined <- p_projected | density_plot


save_plot(
  plot_object = density_plot_combined,
  file_stub = "density_plot",
  save_dir = plot_dir,
  width = 14,
  height = 7
)






Idents(mega_obj) <- "seurat_cluster.sketched"

cluster_colors <- extract_dimplot_colors(
  mega_obj,
  reduction = "umap.sketch",
  group.by = "seurat_cluster.sketched"
)

p_cluster_tree <- plot_cluster_tree(
  mega_obj,
  reduction = embedding_reduction,
  ndims = max_dims,
  cluster_colors = cluster_colors
)

cluster_tree_combined <- p_projected | p_cluster_tree



# p_cluster_tree <- plot_cluster_tree(
#   object = mega_obj,
#   reduction = "harmony.sketch",
#   ndims = 30
# )
# 
# p_cluster_tree


save_plot(
  plot_object = cluster_tree_combined,
  file_stub = "cluster_tree",
  save_dir = plot_dir,
  width = 12,
  height = 7
)














#TODO make these exclude the top .5% so its a better visiuallization?
#additional plots via mega_umapper_plotter
plot_umap_features(
  object = mega_obj,
  features = c(
    paste0("nCount_", DefaultAssay(mega_obj)),
    paste0("nFeature_", DefaultAssay(mega_obj)),
    "percent.mt"
  ),
  save_dir = file.path(
    plot_dir,
    "QC"
  )
)


# plot_umap_metadata(
#   object = mega_obj,
#   variables = c(
#     "experimental_group_tissue",
#     "tissue"
#   ),
#   save_dir = file.path(
#     plot_dir,
#     "Metadata"
#   )
# )

# plot_umap_features(
#   object = mega_obj,
#   features = c(
#     "Epcam",
#     "Krt8",
#     "Krt19",
#     "Muc1"
#   ),
#   save_dir = file.path(
#     plot_dir,
#     "Markers"
#   )
# )

