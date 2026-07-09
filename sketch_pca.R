# ==========================================================
# Sketch PCA
#
# Responsibilities
# ----------------
# * Perform sketch-based PCA on a Seurat object.
# * Create the sketch assay.
# * Compute PCA on the sketched data.
#
# This function intentionally does NOT:
#   * Save objects
#   * Load cached objects
#   * Generate figures
#   * Manage directories
#
# Input
# -----
# A Seurat object containing a normalized assay.
#
# Output
# ------
# Returns the Seurat object with:
#   * sketch assay
#   * pca.sketch dimensional reduction
#
# Typical usage
# -------------
# object <- load_or_build_object(
#   directory = pca_dir,
#   build_function = function() {
#     run_sketch_pca(object)
#   }
# )
# ==========================================================

run_sketch_pca <- function(
    object,
    ncells = 50000,
    full_assay = NULL,
    sketch_assay = "sketch",
    sketch_reduction = "pca.sketch"
) {
  
  # ==========================================================
  # Determine assay to sketch
  # ==========================================================
  
  if (is.null(full_assay)) {
    full_assay <- DefaultAssay(object)
  }
  
  DefaultAssay(object) <- full_assay
  
  # ==========================================================
  # Prepare the full assay
  #
  # Variable features are used during leverage-score sketching.
  # ==========================================================
  
  message("Finding variable features...")
  
  object <- FindVariableFeatures(object)
  
  message("Scaling full assay...")
  
  object <- ScaleData(object)
  
  # ==========================================================
  # Construct the sketch assay
  # ==========================================================
  
  message(
    "Creating sketch assay (",
    ncells,
    " cells)..."
  )
  
  object <- SketchData(
    object = object,
    ncells = ncells,
    method = "LeverageScore",
    sketched.assay = sketch_assay
  )
  
  # ==========================================================
  # Analyze the sketch assay
  # ==========================================================
  
  DefaultAssay(object) <- sketch_assay
  
  message("Finding sketch variable features...")
  
  object <- FindVariableFeatures(object)
  
  message("Scaling sketch assay...")
  
  object <- ScaleData(object)
  
  # ==========================================================
  # Compute PCA
  # ==========================================================
  
  message("Running sketch PCA...")
  
  object <- RunPCA(
    object = object,
    reduction.name = sketch_reduction
  )
  
  message("Sketch PCA complete.")
  
  return(object)
  
}