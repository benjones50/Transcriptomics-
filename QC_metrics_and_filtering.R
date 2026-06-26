# ============================================================
# VISIUM HELPER FUNCTIONS
# ============================================================
#
# This file is meant to be sourced by your main loader script.
# It contains small reusable functions for:
#   1) adding QC metadata
#   2) filtering a Seurat object using QC thresholds
#
# Example use:
#   source(file.path(code_dir, "QC_metrics_and_filtering.R"))
# ============================================================




# ------------------------------------------------------------
# Add QC metrics to a Seurat object
# ------------------------------------------------------------
#
# What this does:
#   - Finds the column names for counts and features
#     based on the active assay name
#   - Adds percent.mt using mitochondrial genes
#
# Why this is useful:
#   - You only need to calculate these QC columns once per raw
#     object
#   - Later steps can reuse the saved QC-annotated object
# ------------------------------------------------------------

add_visium_qc_metrics <- function(object, mt_pattern = "^mt-") {
  
  
  
  # Get the name of the active assay, for example "Spatial"
  assay_name <- DefaultAssay(object)
  
  
  # Build the standard Seurat metadata column names
  # Example:
  #   nCount_Spatial
  #   nFeature_Spatial
  count_col <- paste0("nCount_", assay_name)
  feature_col <- paste0("nFeature_", assay_name)
  

  
  #remover
  object <- subset(object, subset = nCount_Spatial > 0)
  
  
  
  # Calculate mitochondrial percentage for each bin/spot
  # adds a new metadata column called percent.mt
  object[["percent.mt"]] <- PercentageFeatureSet(
    object,
    pattern = mt_pattern
  )
  

  # Return the updated object
  object
}




# ------------------------------------------------------------
# Filter a Seurat object using QC thresholds
# ------------------------------------------------------------
#
# What this does:
#   - Looks at the metadata columns for counts, features, and
#     percent.mt
#   - Keeps only bins/spots/cells that pass all thresholds
#   - Returns a subsetted Seurat object
#
# ------------------------------------------------------------
filter_visium_object_by_qc <- function(
    object,
    min_counts,
    max_counts,
    min_features,
    max_features,
    max_percent_mt
) {
  
  # Get the active assay name so we can build the right metadata column names
  assay_name <- DefaultAssay(object)
  
  # Standard Seurat QC column names
  count_col <- paste0("nCount_", assay_name)
  feature_col <- paste0("nFeature_", assay_name)
  
  # Pull the metadata table for easier filtering
  meta <- object@meta.data
  
  # Identify bins/spots/cells that pass every QC threshold
  #
  # Keep only rows where:
  #   - total counts are above min_counts and below max_counts
  #   - feature counts are above min_features and below max_features
  #   - mitochondrial percentage is below max_percent_mt
  keep_bins <- rownames(meta)[
    meta[[count_col]] > min_counts &
      meta[[count_col]] < max_counts &
      meta[[feature_col]] > min_features &
      meta[[feature_col]] < max_features &
      meta$percent.mt < max_percent_mt
  ]
  
  # Return a Seurat object containing only the kept bins
  subset(object, cells = keep_bins)
}


