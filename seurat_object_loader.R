# ============================================================
# SEURAT OBJECT LOADER
# ============================================================
#
# Goal:
#   Load or create Visium Seurat objects at different stages:
#     - raw
#     - qc
#     - filtered
#
# Design:
#   - Filenames encode the settings used.
#   - If the exact file already exists, load it.
#   - If not, build it, save it, and return it.
#   - No separate manifest file is needed because the filename
#     itself records the settings.
#
# ============================================================



# loads helper functions for object loader
# functions to add qc and filter
source(file.path(code_dir, "QC_metrics_and_filtering.R"))


# ------------------------------------------------------------
# Build the filename for a Visium object version
#
# This encodes the settings for saved visium objects
# ------------------------------------------------------------
visium_object_path <- function(
    sample_output_dir,
    sample_name,
    version = c("raw", "qc", "filtered"),
    bin_size,
    min_counts = NULL,
    max_counts = NULL,
    min_features = NULL,
    max_features = NULL,
    max_percent_mt = NULL
) {
  version <- match.arg(version)
  
  # Start with the common pieces that every filename should have.
  # This makes the object easy to identify later.
  base_name <- paste0(
    "visium_seurat_",
    version,
    "_",
    sample_name,
    "_",
    bin_size,
    "um"
  )
  
  # Add QC thresholds only for filtered objects.
  # These settings become part of the filename so that different
  # filter versions can coexist safely.
  if (version == "filtered") {
    if (is.null(min_counts) || is.null(max_counts) ||
        is.null(min_features) || is.null(max_features) ||
        is.null(max_percent_mt)) {
      stop(
        "For version = 'filtered', you must provide min_counts, max_counts, ",
        "min_features, max_features, and max_percent_mt."
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
  }
  
  # Create a dedicated folder for Seurat objects
  object_dir <- file.path(
    sample_output_dir,
    "objects"
  )
  
  dir.create(
    object_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Return the full file path for the object
  file.path(
    object_dir,
    paste0(base_name, ".rds")
  )

}





# ------------------------------------------------------------
# Main loader
# ------------------------------------------------------------
load_visium_object <- function(
    sample_output_dir,         # where cached objects are stored
    sample_name,               # name used in filenames
    sample_tissue,             # folder name inside raw_data_dir
    raw_data_dir,              # parent directory of Space Ranger outputs
    bin_size,                  # Visium HD bin size, e.g. 2, 8, 16
    version = c("raw", "qc", "filtered"),
    mt_pattern = "^mt-",       # mitochondrial gene prefix
    min_counts = NULL,         # QC filter lower UMI bound
    max_counts = NULL,         # QC filter upper UMI bound
    min_features = NULL,       # QC filter lower feature bound
    max_features = NULL,       # QC filter upper feature bound
    max_percent_mt = NULL,     # QC filter mitochondrial cutoff
    force_rebuild = FALSE      # if TRUE, ignore cache and remake
) {
  version <- match.arg(version)
  
  # Make sure the output directory exists before trying to save there.
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # Determine the expected file path for this exact object
  # ----------------------------------------------------------
  object_file <- visium_object_path(
    sample_output_dir = sample_output_dir,
    sample_name = sample_name,
    version = version,
    bin_size = bin_size,
    min_counts = min_counts,
    max_counts = max_counts,
    min_features = min_features,
    max_features = max_features,
    max_percent_mt = max_percent_mt
  )
  
  # ----------------------------------------------------------
  # If the exact file already exists, load it immediately
  # ----------------------------------------------------------
  if (file.exists(object_file) && !force_rebuild) {
    message("Loading cached ", version, " object:\n  ", object_file)
    return(readRDS(object_file))
  }
  
  # ----------------------------------------------------------
  # Build the raw object
  # ----------------------------------------------------------
  if (version == "raw") {
    message("Creating raw Seurat object from Space Ranger output...")
    
    object <- Load10X_Spatial(
      data.dir = file.path(raw_data_dir, sample_tissue, "outs"),
      bin.size = bin_size
    )
    
    saveRDS(object, object_file)
    message("Saved raw object to:\n  ", object_file)
    return(object)
  }
  
  # ----------------------------------------------------------
  # QC-annotated object
  # ----------------------------------------------------------
  if (version == "qc") {
    message("Loading raw object, then adding QC metadata...")
    
    # First load the raw object.
    raw_object <- load_visium_object(
      sample_output_dir = sample_output_dir,
      sample_name = sample_name,
      sample_tissue = sample_tissue,
      raw_data_dir = raw_data_dir,
      bin_size = bin_size,
      version = "raw",
      mt_pattern = mt_pattern,
      force_rebuild = force_rebuild
    )
    
    # ------------------------------------------------------
    # Uses QC_metrics_and_filtering
    #
    # Adds:
    # - count_col,
    # - feature_col
    # - percent.mt
    # ------------------------------------------------------
    
    
    object <- add_visium_qc_metrics(
      object = raw_object,
      mt_pattern = mt_pattern
    )
    
    saveRDS(object, object_file)
    message("Saved QC-annotated object to:\n  ", object_file)
    return(object)
  }
  
  # ----------------------------------------------------------
  # Filtered object
  # ----------------------------------------------------------
  if (version == "filtered") {
    message("Loading QC object, then applying filters...")
    
    # First load the QC-annotated object.
    qc_object <- load_visium_object(
      sample_output_dir = sample_output_dir,
      sample_name = sample_name,
      sample_tissue = sample_tissue,
      raw_data_dir = raw_data_dir,
      bin_size = bin_size,
      version = "qc",
      mt_pattern = mt_pattern,
      force_rebuild = force_rebuild
    )
    
    # Make sure thresholds were provided.
    if (is.null(min_counts) || is.null(max_counts) ||
        is.null(min_features) || is.null(max_features) ||
        is.null(max_percent_mt)) {
      stop(
        "For version = 'filtered', you must supply all QC thresholds: ",
        "min_counts, max_counts, min_features, max_features, max_percent_mt."
      )
    }
    
    # ------------------------------------------------------
    # filtering helper
    #
    # applies defined filters to object
    # returns only bins matching criteria
    # 
    # ------------------------------------------------------
    
    object <- filter_visium_object_by_qc(
      object = qc_object,
      min_counts = min_counts,
      max_counts = max_counts,
      min_features = min_features,
      max_features = max_features,
      max_percent_mt = max_percent_mt
    )
    
    saveRDS(object, object_file)
    message("Saved filtered object to:\n  ", object_file)
    return(object)
  }
  
  # This should never happen because match.arg() checks version.
  stop("Unknown version requested.")
}




# ============================================================
# EXAMPLE USAGE
# ============================================================
#
#
# raw_object <- load_visium_object(
#   sample_output_dir = sample_output_dir,
#   sample_name = sample_name,
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   bin_size = bin_size,
#   version = "raw"
# )
#
# qc_object <- load_visium_object(
#   sample_output_dir = sample_output_dir,
#   sample_name = sample_name,
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   bin_size = bin_size,
#   version = "qc"
# )
#
# filtered_object <- load_visium_object(
#   sample_output_dir = sample_output_dir,
#   sample_name = sample_name,
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   bin_size = bin_size,
#   version = "filtered",
#   min_counts = min_counts,
#   max_counts = max_counts,
#   min_features = min_features,
#   max_features = max_features,
#   max_percent_mt = max_percent_mt
# )
# ============================================================