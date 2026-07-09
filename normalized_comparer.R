
source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")

#loads objects
force_rebuild <- TRUE


object_sct <- load_normalized_object(
  sample_tissue = sample_tissue,
  raw_data_dir = raw_data_dir,
  analysis_mode = analysis_mode,
  bin_size = bin_size,
  normalization = "sct",
  min_counts = min_counts,
  max_counts = max_counts,
  min_features = min_features,
  max_features = max_features,
  max_percent_mt = max_percent_mt,
  force_rebuild = force_rebuild
)

object_log_norm <- load_normalized_object(
  sample_tissue = sample_tissue,
  raw_data_dir = raw_data_dir,
  analysis_mode = analysis_mode,
  bin_size = bin_size,
  normalization = "lognorm",
  min_counts = min_counts,
  max_counts = max_counts,
  min_features = min_features,
  max_features = max_features,
  max_percent_mt = max_percent_mt,
  force_rebuild = force_rebuild
)

# object <- load_visium_object(
#   analysis_mode = analysis_mode,
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   bin_size = bin_size,
#   stag = "qc",
#   force_rebuild = force_rebuild
# )
# 
# filtered_object <- load_visium_object(
#   sample_tissue = sample_tissue,
#   analysis_mode = analysis_mode,
#   raw_data_dir = raw_data_dir,
#   bin_size = bin_size,
#   stage = "filtered",
#   min_counts = min_counts,
#   max_counts = max_counts,
#   min_features = min_features,
#   max_features = max_features,
#   max_percent_mt = max_percent_mt,
#   force_rebuild = force_rebuild
# )


cat("object_sct:", DefaultAssay(object_sct), "\n")
cat("object_log_norm:", DefaultAssay(object_log_norm), "\n")
cat("object:", DefaultAssay(object), "\n")
cat("filtered_object:", DefaultAssay(filtered_object), "\n")



cat("number of cells in object:", ncol(object_log_norm))



# performs a sketched pca
# returns the object with the sketched pca attached
run_sketch_pca <- function(
    object,
    ncells = 50000,
    full_assay = DefaultAssay(object),
    sketch_assay = "sketch",
    sketch_reduction = "pca.sketch",
    dims = 1:50
) {
  
  # Ensure the correct assay is active
  DefaultAssay(object) <- full_assay
  
  # Pre-sketch processing
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)
  
  # Create sketch assay
  object <- SketchData(
    object = object,
    ncells = ncells,
    method = "LeverageScore",
    sketched.assay = sketch_assay
  )
  
  # Analyze sketch
  DefaultAssay(object) <- sketch_assay
  
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)
  
  
  #pca for the sketched data
  message("Running sketch PCA...")
  object <- RunPCA(
    object,
    reduction.name = sketch_reduction
  )
  
  return(object)
}


#gets the sketched pca's applied to objects
DefaultAssay(object_log_norm) <- "Spatial"
object_log_norm <- run_sketch_pca(object_log_norm)

DefaultAssay(object_sct) <- "sct"
object_sct <- run_sketch_pca(object_sct)




# used to make plots of PCA
plot_pca_diagnostics <- function(
    object,
    reduction = "pca.sketch",
    dims = 1:5,
    ndims = 50
) {
  
  # Visualize gene loadings
  loading_plot <- VizDimLoadings(
    object,
    dims = dims,
    reduction = reduction
  ) #x axis is the genes weight for the PC
  
  
  # Elbow plot
  elbow_plot <- ElbowPlot(
    object,
    reduction = reduction,
    ndims = ndims
  )
  
  return(list(
    loading_plot = loading_plot,
    elbow_plot = elbow_plot
  ))
  
}

# plot the PCA info for log_norm and sct

# decide on DIMS for analysis 

lognorm_pca <- plot_pca_diagnostics(object_log_norm)

lognorm_pca$loading_plot
lognorm_pca$elbow_plot


sct_pca <- plot_pca_diagnostics(object_sct)

sct_pca$loading_plot
sct_pca$elbow_plot


sct_pca$elbow_plot / lognorm_pca$elbow_plot






#clustering 
run_sketch_clustering <- function(
    object,
    full_assay = DefaultAssay(object),
    sketch_assay = "sketch",
    sketch_reduction = "pca.sketch",
    projected_reduction = "full.umap.sketch",
    dims = 1:30,
    sketch_umap = "umap.sketch",
    resolution = 3,
    cluster_name = "seurat_cluster.sketched"
) {
  
  
  # Cluster the sketch
  message("Finding sketch neighbors...")
  
  object <- FindNeighbors(
    object,
    assay = sketch_assay,
    reduction = sketch_reduction,
    dims = dims
  )
  
  message("Finding sketch clusters...")
  
  object <- FindClusters(
    object,
    cluster.name = cluster_name,
    resolution = resolution
  )
  
  
  # Compute UMAP on the sketch
  message("Running sketch UMAP...")
  
  object <- RunUMAP(
    object,
    reduction = sketch_reduction,
    reduction.name = sketch_umap,
    return.model = TRUE,
    dims = dims
  )
  
  
  # Determine normalization method
  if (inherits(object[[full_assay]], "SCTAssay")) {
    normalization_method <- "SCT"
  } else if (inherits(object[[full_assay]], "Assay5")) {
    normalization_method <- "LogNormalize"
  } else {
    stop(
      "Unsupported assay class: ",
      class(object[[full_assay]])[1],
      ". Expected an SCTAssay or Assay5."
    )
  }
  
  
  # Project PCA, clusters, and UMAP
  message("Projecting data to full dataset...")
  
  object <- ProjectData(
    object = object,
    assay = full_assay,
    sketched.assay = sketch_assay,
    sketched.reduction = sketch_reduction,
    full.reduction = projected_reduction,
    normalization.method = normalization_method,
    umap.model = sketch_umap,
    dims = dims
  )
  
  message("Finished projection.")
  
  return(object)
}


object_log_norm <- run_sketch_clustering(object_log_norm)

object_sct <- run_sketch_clustering(object_sct)




# ------------------------------------------------------------
# Plot sketch vs. projected UMAPs
#
# Creates a side-by-side comparison of:
#   1. The UMAP generated from the sketched cells.
#   2. The UMAP after projecting all cells into the sketch UMAP.
#
# This is useful for verifying that the projected embedding
# faithfully reproduces the structure learned from the sketch.
# ------------------------------------------------------------
plot_sketch_projection <- function(
    object,
    sketch_assay = "sketch",
    sketch_clusters = "seurat_cluster.sketched",
    projected_clusters = "seurat_cluster.projected",
    sketch_reduction = "umap.sketch",
    projected_reduction = "full.umap.sketch",
    sketch_title = "Sketch clustering",
    full_title = "Projected clustering"
) {
  
  # ----------------------------------------------------------
  # Save the original assay (typically "Spatial" or "sct")
  #
  # The plotting function temporarily switches between assays,
  # so we remember the original one to restore it before
  # plotting the projected UMAP, may be useful for other functions
  # ----------------------------------------------------------
  full_assay <- DefaultAssay(object)
  
  # ==========================================================
  # Plot the sketch UMAP
  # ==========================================================
  
  # Switch to the sketch assay
  DefaultAssay(object) <- sketch_assay
  
  # Color cells using the sketch cluster assignments
  Idents(object) <- sketch_clusters
  
  # Plot the UMAP computed directly on the sketched cells
  p1 <- DimPlot(
    object,
    reduction = sketch_reduction,
    label = FALSE
  ) +
    ggtitle(sketch_title) +
    theme(legend.position = "bottom")
  
  # ==========================================================
  # Plot the projected full-dataset UMAP
  # ==========================================================
  
  # Restore the original assay
  DefaultAssay(object) <- full_assay
  
  # Color cells using the projected cluster labels
  # assigned during ProjectData()
  Idents(object) <- projected_clusters
  
  # Plot the UMAP after projecting every cell into the
  # sketch-derived UMAP embedding
  p2 <- DimPlot(
    object,
    reduction = projected_reduction,
    label = FALSE
  ) +
    ggtitle(full_title) +
    theme(legend.position = "bottom")
  
  # ==========================================================
  # Return the two plots side-by-side
  # ==========================================================
  return(p1 | p2)
}


object_log_norm_umap <- plot_sketch_projection(object_log_norm)
object_log_norm_umap

object_sct_umap <- plot_sketch_projection(object_sct)
object_sct_umap


SpatialDimPlot(object_log_norm, label = T, repel = T, label.size = 4)



#bug fix stuff, ignore

# 
#  object <- Load10X_Spatial(
#    data.dir = file.path(raw_data_dir, sample_tissue, "outs"),
#    bin.size = bin_size,
#    slice = "slice1"
#  )
# 
# object <- load_visium_object(
#   sample_output_dir = sample_output_dir,
#   sample_name = sample_name,
#   analysis_mode = analysis_mode,
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   bin_size = bin_size,
#   stage = "raw",
#   force_rebuild = FALSE
# )
# 
# SpatialFeaturePlot(
#   object,
#   features = "Epcam",
#   slot = "counts"
# )












