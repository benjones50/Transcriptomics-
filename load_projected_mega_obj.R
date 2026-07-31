source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")

#loads an object that has already been projected, useful for analysis past mega object with multiple resolutions



sketch_label <- "sketch_33pct"


pca_dir <- get_pca_dir(
  root_dir = mega_dir,
  ndims = 50,
  sketch_label = sketch_label
)

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
  
  embedding_run_dir <- harmony_dir
  
} else {
  
  embedding_run_dir <- pca_dir
  
}


resolution <- 0.3

clustering_dir_normal <- get_clustering_dir(
  embedding_run_dir = embedding_run_dir,
  ndims = max_dims,
  resolution = resolution
)

clustering_dir <- get_nested_subdir(
  clustering_dir_normal,
  "multi_res_0.3"
)


umap_dir <- get_umap_dir(
  clustering_dir = clustering_dir
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





