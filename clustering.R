# NOTE --------------------------------------------------------------------
# These functions assume `run_sketch_pca()` has already been run.
# The sketch assay and `pca.sketch` reduction are treated as immutable
# downstream inputs and are never recreated here.
#
# -------------------------------------------------------------------------

# ==========================================================
# Internal helper: resolve the full assay
#
# Problem:
# Downstream projection needs the original assay, but the
# default assay may already be set to the sketch assay.
#
# If the user does not explicitly specify `full_assay`,
# automatically choose the first assay that is NOT the sketch.
#
# This works for:
#   Spatial.016um
#   Spatial.008um
#   Spatial
#   RNA
#   SCT
#   segmented-cell assays
# ==========================================================
resolve_full_assay <- function(object, sketch_assay, full_assay = NULL) {
  
  if (is.null(full_assay)) {
    
    full_assay <- setdiff(
      Assays(object),
      sketch_assay
    )[1]
  }
  
  
  if (is.na(full_assay)) {
    stop(
      "Could not determine the full assay automatically.",
      call. = FALSE
    )
  }
  
  
  if (identical(full_assay, sketch_assay)) {
    stop(
      "full_assay cannot be the same as sketch_assay ('", sketch_assay, "').",
      call. = FALSE
    )
  }
  
  return(full_assay)
}


# ==========================================================
# Internal helper: resolve the normalization method
#
# ProjectData() needs to know whether the original assay
# came from SCTransform or standard log normalization.
#
# This is determined from the assay class so the projection
# stage can work across both workflows.
# ==========================================================
resolve_normalization_method <- function(object, full_assay) {
  
  full_assay_object <- object[[full_assay]]
  
  if (inherits(full_assay_object, "SCTAssay")) {
    
    return("SCT")
    
  } else if (inherits(full_assay_object, "Assay5")) {
    
    return("LogNormalize")
    
  } else {
    
    stop(
      "Unsupported assay class for full_assay '",
      full_assay,
      "': ",
      class(full_assay_object)[1],
      call. = FALSE
    )
    
  }
}






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
    sketch_reduction = "pca.sketch",
    dims = 1:30,
    resolution = 3,
    cluster_name = "seurat_cluster.sketched"
) {
  

  

  
  message("--------------------------------------------")
  message("Sketch assay  : ", sketch_assay)
  message("PCA reduction : ", sketch_reduction)
  message("Dimensions    : ", paste(range(dims), collapse = "-"))
  message("Resolution    : ", resolution)
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
  
  message("Finding clusters...")
  
  object <- FindClusters(
    object,
    cluster.name = cluster_name,
    resolution = resolution
  )
  
  
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
    sketch_reduction = "pca.sketch",
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
    full_assay = NULL,
    sketch_assay = "sketch",
    sketch_reduction = "pca.sketch",
    projected_reduction = "projected.pca",
    sketch_umap = "umap.sketch",
    dims = 1:30,
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
  message("Dimensions    : ", paste(range(dims), collapse = "-"))
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
  
  normalization_method <- resolve_normalization_method(
    object = object,
    full_assay = full_assay
  )
  
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
    dims = dims
  )
  
  
  # NOTE: to keep the full assay
  # active for the next downstream analysis step, keeping
  # `DefaultAssay(object) <- full_assay` here is the cleanest
  # choice. It makes the returned object immediately ready for
  # full-assay-based plotting or differential expression.
  
  message("Sketch projection complete.")
  
  return(object)
}










