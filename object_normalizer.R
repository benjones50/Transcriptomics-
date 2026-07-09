# ============================================================
# NORMALIZED SEURAT OBJECT LOADER
# ============================================================
#
# Goal:
#   Load or create normalized Seurat objects from filtered
#   Visium objects.
#
# Behavior:
#   - Search for an existing normalized object on disk.
#   - If found, load and return it.
#   - If not found, load the filtered object, normalize it,
#     save it to a normalized subfolder, and return it.
#
# Supported normalization methods:
#   - "lognorm"  : Seurat NormalizeData(LogNormalize)
#   - "sct"      : Seurat SCTransform
#
# ============================================================

#object_normal <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)

#object_SCT <- SCTransform(object,assay = "Spatial", new.assay.name = "SCT", vars.to.regress = "percent.mt", verbose = FALSE)



# ------------------------------------------------------------
# Normalize a Seurat object using the chosen method
# ------------------------------------------------------------
normalize_visium_object <- function(
    object,
    normalization = c("lognorm", "sct"),
    vars.to.regress = NULL, #"percent.mt",
    sct_assay_name = "SCT",
    scale.factor = 10000,
    verbose = FALSE
) {
  normalization <- match.arg(normalization)
  
  
  # ----------------------------------------------------------
  # Log-normalization
  # ----------------------------------------------------------
  if (normalization == "lognorm") {
    message(
      "Running LogNormalize"
    )
    
    object <- NormalizeData(
      object = object,
      normalization.method = "LogNormalize",
      scale.factor = scale.factor,
      verbose = verbose
    )
    
    message("Log-normalization complete.")
    return(object)
  }
  
  # ----------------------------------------------------------
  # SCTransform normalization
  # ----------------------------------------------------------
  if (normalization == "sct") {
    message(
      "Running SCTransform()"
    )
    
    assay_name <- DefaultAssay(object)
    
    object <- SCTransform(
      object = object,
      assay = assay_name,
      new.assay.name = sct_assay_name,
      vars.to.regress = vars.to.regress,
      verbose = verbose
    )
      
      
    # Make SCT the default assay for downstream analysis.
    DefaultAssay(object) <- sct_assay_name
    
    message(
      "SCTransform complete. Default assay set to '", sct_assay_name, "'."
    )
    
    return(object)
  }
  
  stop("Unknown normalization method requested.")
}







# ------------------------------------------------------------
# Main normalized-object loader
#
# This mirrors the behavior of existing seurat_object loader:
#   1) check for cached normalized object
#   2) if present, load it
#   3) if not, load filtered object, normalize, save, return
# ------------------------------------------------------------
load_normalized_object <- function( #TODO add vars to regress as a variable here
    sample_tissue,
    raw_data_dir,
    analysis_mode = c("binned", "segmented_cells"),
    bin_size = NULL,
    normalization = c("lognorm", "sct"),
    min_counts = NULL,
    max_counts = NULL,
    min_features = NULL,
    max_features = NULL,
    max_percent_mt = NULL,
    mt_pattern = "^mt-",
    force_rebuild = FALSE,
    verbose = FALSE
) {
  analysis_mode <- match.arg(analysis_mode)
  normalization <- match.arg(normalization)
  
  # ----------------------------------------------------------
  # Print a clear header so the console shows what is happening.
  # ----------------------------------------------------------
  message(
    "\n==============================\n",
    "Requested normalized object\n",
    "analysis_mode = ", analysis_mode, "\n",
    "normalization = ", normalization,
    if (!is.null(bin_size)) paste0("\nbin_size      = ", bin_size, "um") else "",
    "\n=============================="
  )
  
  # ----------------------------------------------------------
  # Validate analysis mode / bin size combinations.
  # ----------------------------------------------------------
  if (analysis_mode == "binned" && is.null(bin_size)) {
    stop("bin_size must be supplied when analysis_mode = 'binned'")
  }
  
  if (analysis_mode == "segmented_cells" && !is.null(bin_size)) {
    stop("bin_size must be NULL when analysis_mode = 'segmented_cells'")
  }
  
  # ----------------------------------------------------------
  # Build the exact normalized object path for these settings.
  # ----------------------------------------------------------
  object_file <- build_object_path(
    sample_tissue = sample_tissue,
    stage = normalization,
    analysis_mode = analysis_mode,
    bin_size = bin_size,
    min_counts = min_counts,
    max_counts = max_counts,
    min_features = min_features,
    max_features = max_features,
    max_percent_mt = max_percent_mt
  )
  
  # ----------------------------------------------------------
  # If the normalized object already exists, load it and return.
  # ----------------------------------------------------------
  if (file.exists(object_file) && !force_rebuild) {
    message(
      "Loading cached ",
      normalization,
      " object:\n  ",
      object_file
    )
    return(readRDS(object_file))
  }
  
  # ----------------------------------------------------------
  # Otherwise, build the normalized object from the filtered one.
  # ----------------------------------------------------------
  message(
    "No cached normalized object found. Building a new ",
    normalization,
    " object..."
  )
  
  message(
    "\n==============================\n",
    "Filtered object settings\n",
    "min_counts     = ", min_counts, "\n",
    "max_counts     = ", max_counts, "\n",
    "min_features   = ", min_features, "\n",
    "max_features   = ", max_features, "\n",
    "max_percent_mt = ", max_percent_mt, "\n",
    "=============================="
  )
  
  # Load the filtered object from existing loader.
  # This keeps the normalized loader separate from the raw/qc/filter pipeline.
  #uses load_visium_object from the seurat object loaders
  filtered_object <- load_visium_object(
    sample_tissue = sample_tissue,
    raw_data_dir = raw_data_dir,
    bin_size = bin_size,
    analysis_mode = analysis_mode,
    stage = "filtered",
    mt_pattern = mt_pattern,
    min_counts = min_counts,
    max_counts = max_counts,
    min_features = min_features,
    max_features = max_features,
    max_percent_mt = max_percent_mt,
    force_rebuild = force_rebuild
  )
  
  message(
    "Filtered object loaded successfully. Beginning normalization..."
  )
  
  # ----------------------------------------------------------
  # Run the chosen normalization method.
  # ----------------------------------------------------------
  object <- normalize_visium_object(
        object = filtered_object,
        normalization = normalization,
        vars.to.regress = NULL, #"percent.mt",
        sct_assay_name = "sct",
        scale.factor = 10000,
        verbose = verbose
      )
    
  
  # ----------------------------------------------------------
  # Save the normalized object for fast reloading later.
  # ----------------------------------------------------------
  
  message(
    "Saving normalized object..."
  )
  
  saveRDS(object, object_file)
  
  message(
    "Normalized object saved to:\n  ",
    object_file
  )
  
  return(object)
}






















