# ------------------------------------------------------------
# Build the cache path for a Seurat object
#
# Given a sample and processing stage, determine where the
# object should live on disk.
#
# This function centralizes ALL object naming so every loader
# (raw, filtered, normalized, etc.) follows the same directory
# structure.
# ------------------------------------------------------------
build_object_path <- function(
    sample_tissue,
    stage = c(
      "raw",
      "qc",
      "filtered",
      "lognorm",
      "sct"
    ),
    analysis_mode = c(
      "binned",
      "segmented_cells"
    ),
    bin_size = NULL,
    min_counts = NULL,
    max_counts = NULL,
    min_features = NULL,
    max_features = NULL,
    max_percent_mt = NULL
) {
  
  stage <- match.arg(stage)
  analysis_mode <- match.arg(analysis_mode)
  
  # ==========================================================
  # Determine sample name and cache directory
  # ==========================================================
  
  if (analysis_mode == "binned") {
    
    sample_name <- paste0(
      bin_size,
      "um_",
      sample_tissue
    )
    
    sample_object_dir <- file.path(
      output_dir,
      paste0(bin_size, "um"),
      sample_name
    )
    
    base_name <- paste0(
      "object_",
      stage,
      "_",
      sample_name,
      "_",
      bin_size,
      "um"
    )
    
  } else {
    
    sample_name <- paste0(
      "seg_",
      sample_tissue
    )
    
    sample_object_dir <- file.path(
      output_dir,
      "segmented_cells",
      sample_name
    )
    
    base_name <- paste0(
      "cell_segmented_",
      stage,
      "_",
      sample_name
    )
    
  }
  
  # ==========================================================
  # Filtered objects encode QC thresholds
  # ==========================================================
  
  if (stage %in% c("filtered","lognorm","sct")) {
    
    if (
      is.null(min_counts) ||
      is.null(max_counts) ||
      is.null(min_features) ||
      is.null(max_features) ||
      is.null(max_percent_mt)
    ) {
      stop(...)
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
    
  
  
  # ==========================================================
  # Create cache directory if needed
  # ==========================================================
  
  object_dir <- file.path(
    sample_object_dir,
    "objects"
  )
  
  dir.create(
    object_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ==========================================================
  # Return full cache path
  # ==========================================================
  
  file.path(
    object_dir,
    paste0(base_name, ".rds")
  )

}


