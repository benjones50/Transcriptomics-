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
# Build the filename for a normalized object
#
# This encodes the key settings used to generate the object,
# so different normalization choices and QC settings can coexist.
# ------------------------------------------------------------
normalized_object_path <- function(
    sample_output_dir,
    sample_name,
    normalization = c("lognorm", "sct"),
    analysis_mode = c("binned", "segmented_cells"),
    bin_size = NULL,
    min_counts = NULL,
    max_counts = NULL,
    min_features = NULL,
    max_features = NULL,
    max_percent_mt = NULL
) {
  normalization <- match.arg(normalization)
  analysis_mode <- match.arg(analysis_mode)
  
  # Put all normalized objects in a dedicated subfolder.
  norm_dir <- file.path(sample_output_dir, "objects", "normalized")
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build a readable filename that captures the object type and settings.
  if (analysis_mode == "binned") {
    
    if (is.null(bin_size)) {
      stop("bin_size must be supplied when analysis_mode = 'binned'")
    }
    
    base_name <- paste0(
      "object_normal_",
      sample_name,
      "_",
      bin_size,
      "um"
    )
    
  } else if (analysis_mode == "segmented_cells") {
    
    base_name <- paste0(
      "cell_seg_normal_",
      sample_name
    )
    
  } else {
    stop("Unknown analysis_mode supplied to normalized_object_path().")
  }
  
  # If this normalized object is based on a filtered object, encode
  # the QC thresholds so the cache key is unambiguous.
  if (
    is.null(min_counts) || is.null(max_counts) ||
    is.null(min_features) || is.null(max_features) ||
    is.null(max_percent_mt)
  ) {
    stop(
      "For normalized objects, you must provide min_counts, max_counts, ",
      "min_features, max_features, and max_percent_mt so the cache key ",
      "matches the filtered object it was built from."
    )
  }
  
  base_name <- paste0(
    base_name,
    "_minC", min_counts,
    "_maxC", max_counts,
    "_minF", min_features,
    "_maxF", max_features,
    "_maxMT", max_percent_mt
  )
  
  # Add the normalization method to the filename.
  if (normalization == "lognorm") {
    base_name <- paste0(base_name, "_lognorm")
  } else if (normalization == "sct") {
    base_name <- paste0(base_name, "_sct")
  }
  
  file.path(norm_dir, paste0(base_name, ".rds"))
}






# ------------------------------------------------------------
# Normalize a Seurat object using the chosen method
# ------------------------------------------------------------
normalize_visium_object <- function(
    object,
    normalization = c("lognorm", "sct"),
    vars.to.regress = "percent.mt",
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
load_normalized_object <- function(
    sample_output_dir,
    sample_name,
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
  
  # Make sure the output directory exists.
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # Build the exact normalized object path for these settings.
  # ----------------------------------------------------------
  object_file <- normalized_object_path(
    sample_output_dir = sample_output_dir,
    sample_name = sample_name,
    normalization = normalization,
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
  filtered_object <- load_visium_object(
    sample_output_dir = sample_output_dir,
    sample_name = sample_name,
    sample_tissue = sample_tissue,
    raw_data_dir = raw_data_dir,
    bin_size = bin_size,
    analysis_mode = analysis_mode,
    version = "filtered",
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
        vars.to.regress = "percent.mt",
        sct_assay_name = "SCT",
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






















