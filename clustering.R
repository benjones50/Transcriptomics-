# NOTE --------------------------------------------------------------------
# These functions assume `run_sketch_pca()` has already been run.
# The sketch assay and `pca.sketch` reduction are treated as immutable
# downstream inputs and are never recreated here.
#
# -------------------------------------------------------------------------


# ==========================================================
# Stage 1: Graph construction + clustering only
#
# Responsibilities
# ----------------
# * Use the sketch assay
# * Build the neighbor graph
# * Find clusters
# * Do not compute UMAP
# * Do not project to the full dataset
#
# Output
# ------
# Adds:
#   * Neighbor graph
#   * seurat_cluster.sketched metadata
# ==========================================================
run_sketch_clustering <- function(
    object,
    sketch_assay = "sketch",
    sketch_reduction,
    dims,
    resolution, #this is weird bc updated to accept multiple resolutions while backwards compatible with using 1
    resolutions = resolution,
    selected_cluster_resolution = NULL,
    cluster_name = "seurat_cluster.sketched"#,
    #random_seed = 42
) {
  
  
  if (is.null(selected_cluster_resolution)) {
    selected_cluster_resolution <- resolution
  }
  
  if (!selected_cluster_resolution %in% resolutions) {
    stop(
      "selected_cluster_resolution must be one of the values in `resolutions`."
    )
  }
  
  # Convert to character for safe cluster-name construction
  resolutions_chr <- as.character(resolutions)
  selected_cluster_chr <- as.character(selected_cluster_resolution)
  
  message("--------------------------------------------")
  message("Sketch assay              : ", sketch_assay)
  message("PCA reduction             : ", sketch_reduction)
  message("Dimensions                : ", paste(range(dims), collapse = "-"))
  message("Scan resolutions          : ", paste(resolutions_chr, collapse = ", "))
  message("Selected cluster resolution: ", selected_cluster_chr)
  message("--------------------------------------------")
  
  # ==========================================================
  # Build the graph and cluster the sketch
  #
  # Clustering is performed ONLY on the sketch assay.
  # ==========================================================
  
  DefaultAssay(object) <- sketch_assay
  
  message("Finding nearest neighbors...")
  
  object <- FindNeighbors(
    object,
    assay = sketch_assay,
    reduction = sketch_reduction,
    dims = dims
  )
  
  message("Neighbor graph(s):")
  print(Graphs(object))
  
  message("Default assay:")
  print(DefaultAssay(object))
  
  message("Reductions:")
  print(Reductions(object))
  
  message("Graphs available:")
  print(Graphs(object))
  
  message("Finding clusters across resolutions...")
  
  #for multiple clusters
  cluster_prefix <- "cluster_res_"
  
  n_res <- length(resolutions)
  start_time <- Sys.time()
  
  for (i in seq_along(resolutions)) {
    
    res <- resolutions[i]
    
    message(
      sprintf(
        "[%d/%d] Clustering at resolution %.2f",
        i, n_res, res
      )
    )
    
    object <- FindClusters(
      object,
      cluster.name = paste0(cluster_prefix, res),
      resolution = res#,
      #random.seed = random_seed
    )
    
    
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    avg <- elapsed / i
    eta <- avg * (n_res - i)
    
    message(sprintf("Elapsed: %.1f min | ETA: %.1f min", elapsed, eta))
  }
  
  # ----------------------------------------------------------
  # Promote selected resolution to canonical downstream column
  # ----------------------------------------------------------
  selected_cluster_col <- paste0(cluster_prefix, selected_cluster_chr)
  
  if (!selected_cluster_col %in% colnames(object[[]])) {
    stop(
      "Selected cluster column not found: ",
      selected_cluster_col
    )
  }
  
  object[[cluster_name]] <- object[[selected_cluster_col]]
  Idents(object) <- cluster_name
  

  
  message("Sketch graph construction and clustering complete.")
  return(object)
}







# ==========================================================
# Stage 2: Sketch UMAP only
#
# Responsibilities
# ----------------
# * Use the sketch assay
# * Use the existing sketch PCA
# * Compute UMAP
# * Save the fitted model with return.model = TRUE
# * Do not rebuild graph/clusters
# * Do not project to the full dataset
#
# Output
# ------
# Adds:
#   * umap.sketch reduction
#   * fitted UMAP model
# ==========================================================
run_sketch_umap <- function(
    object,
    sketch_assay = "sketch",
    sketch_reduction,
    sketch_umap = "umap.sketch",
    dims = 1:30
) {
  

  
  message("--------------------------------------------")
  message("Sketch assay  : ", sketch_assay)
  message("PCA reduction : ", sketch_reduction)
  message("UMAP name     : ", sketch_umap)
  message("Dimensions    : ", paste(range(dims), collapse = "-"))
  message("--------------------------------------------")
  
  
  # ==========================================================
  # Compute the sketch UMAP
  #
  # return.model = TRUE is required because the fitted UMAP
  # model will later be used by ProjectData() to embed every
  # cell/bin in the full dataset.
  # ==========================================================
  
  DefaultAssay(object) <- sketch_assay
  
  message("Running sketch UMAP...")
  
  object <- RunUMAP(
    object,
    reduction = sketch_reduction,
    reduction.name = sketch_umap,
    return.model = TRUE,
    dims = dims
  )
  

  
  # NOTE: The fitted model is stored inside the UMAP reduction
  # object. Keeping this stage separate makes it easy to reuse
  # the same UMAP model for projection without rerunning the
  # graph or clustering steps.
  
  message("Sketch UMAP complete.")
  
  return(object)
}









# ==========================================================
# Stage 3: Project sketch analysis back to full dataset
#
# Responsibilities
# ----------------
# * Resolve the full assay automatically if needed
# * Detect normalization method from assay class
# * Project sketch PCA / clusters / UMAP back to the full set
# * Do not rebuild graph, clusters, or sketch UMAP
#
# Output
# ------
# Adds:
#   * projected.pca reduction
#   * seurat_cluster.projected metadata
#   * projected embedding in the sketch UMAP
# ==========================================================
run_sketch_projection <- function(
    object,
    normalization_method,
    full_assay = NULL,
    sketch_assay = "sketch",
    sketch_reduction,
    projected_reduction = "projected.pca",
    sketch_umap = "umap.sketch",
    ndims,
    cluster_name = "seurat_cluster.sketched"
) {
  
  # ==========================================================
  # Resolve the full assay
  #
  # This is the original assay that should receive the
  # projected analysis.
  # ==========================================================
  
  full_assay <- resolve_full_assay(
    object = object,
    sketch_assay = sketch_assay,
    full_assay = full_assay
  )
  

  
  message("--------------------------------------------")
  message("Full assay    : ", full_assay)
  message("Sketch assay  : ", sketch_assay)
  message("PCA reduction : ", sketch_reduction)
  message("UMAP model    : ", sketch_umap)
  message("Projection    : ", projected_reduction)
  message("Dimensions    : ", ndims)
  message("--------------------------------------------")
  
  
  # ==========================================================
  # Determine the normalization method
  #
  # ProjectData() needs to know whether the original assay
  # came from SCTransform or standard log normalization.
  #
  # This is determined automatically from the assay class,
  # so the function works across both workflows.
  # ==========================================================
  
  message("Normalization method: ", normalization_method)
  
  
  # ==========================================================
  # Project the sketch analysis back onto every cell/bin
  #
  # This computes:
  #   * projected.pca
  #   * projected clusters
  #   * projected UMAP/model-based embedding
  #
  # for every observation in the original assay.
  # ==========================================================
  
  DefaultAssay(object) <- full_assay
  
  message("Projecting sketch analysis to full dataset...")
  
  object <- ProjectData(
    object = object,
    assay = full_assay,
    sketched.assay = sketch_assay,
    sketched.reduction = sketch_reduction,
    full.reduction = projected_reduction,
    normalization.method = normalization_method,
    refdata = list(
      seurat_cluster.projected = cluster_name
    ),
    umap.model = sketch_umap,
    dims = 1:ndims
  )
  
  
  # NOTE: to keep the full assay
  # active for the next downstream analysis step, keeping
  # `DefaultAssay(object) <- full_assay` here is the cleanest
  # choice. It makes the returned object immediately ready for
  # full-assay-based plotting or differential expression.
  
  message("Sketch projection complete.")
  
  return(object)
}










