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


#TODO add project names to everything so its easier to merge seurats together, bulk analysis, etc


# loads helper functions for object loader
# functions to add qc and filter
source(file.path(code_dir, "QC_metrics_and_filtering.R"))
source(file.path(code_dir,"cell_segmentation_custom.R"))
source(file.path(code_dir,"object_path_builder.R"))







# ------------------------------------------------------------
# Main loader
load_visium_object <- function(
    sample_tissue,
    raw_data_dir,
    bin_size = NULL,
    analysis_mode = c("binned", "segmented_cells"),
    stage = c("raw", "qc", "filtered"),
    mt_pattern = "^mt-",
    min_counts = NULL,
    max_counts = NULL,
    min_features = NULL,
    max_features = NULL,
    max_percent_mt = NULL,
    force_rebuild = FALSE
) {
  stage <- match.arg(stage)
  analysis_mode <- match.arg(analysis_mode)
  
  
  
  # tell user exactly what object is being requested
  message(
    "\n==============================\n",
    "Requested object\n",
    "analysis_mode = ", analysis_mode, "\n",
    "stage       = ", stage,
    if (!is.null(bin_size))
      paste0("\nbin_size      = ", bin_size, "um")
    else "",
    "\n=============================="
  )
  
  
  
  #double check that a bin size is given
  if (
    analysis_mode == "binned" &&
    is.null(bin_size)
  ) {
    stop(
      "bin_size must be supplied when analysis_mode = 'binned'"
    )
  }
  
  #error message if someone put segmented cells with a bin size
  if (
    analysis_mode == "segmented_cells" &&
    !is.null(bin_size)
  ) {
    stop(
      "bin_size must be null for segmented cells"
    )
  }
  
  #builds the storage/output path of the visium object
  object_file <- build_object_path(
    sample_tissue = sample_tissue,
    stage = stage,
    analysis_mode = analysis_mode,
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
    #changes message
    if (analysis_mode == "binned") {
      message(
        "Loading cached ",
        stage,
        " ",
        bin_size,
        "um object:\n  ",
        object_file
      )
    } else {
      message(
        "Loading cached segmented-cell ",
        stage,
        " object:\n  ",
        object_file
      )
    }
    
    return(readRDS(object_file))
  }
  
  
  # ----------------------------------------------------------
  # Binned workflow
  # ----------------------------------------------------------
  if (analysis_mode == "binned") {
    
    if (stage == "raw") {
      
      message(
        "Building raw ",
        analysis_mode,
        " object..."
      )
      
      object <- Load10X_Spatial(
        data.dir = file.path(raw_data_dir, sample_tissue, "outs"),
        bin.size = bin_size
      )
      
      # ----------------------------------------------------------
      # Build binned object manually
      #
      # We avoid Load10X_Spatial() because Seurat currently
      # fails when constructing VisiumV2 objects in this
      # environment.
      # ----------------------------------------------------------
      
      # object <- build_binned_visium_object(
      #   sample_id = sample_tissue,
      #   bin_size = bin_size
      # )
      
      
      
      
      
      saveRDS(object, object_file)
      return(object)
    }
    
    
    
    if (stage == "qc") {
      
      message(
        "Building qc ",
        analysis_mode,
        " object..."
      )
      
      raw_object <- load_visium_object(
        sample_tissue = sample_tissue,
        raw_data_dir = raw_data_dir,
        bin_size = bin_size,
        analysis_mode = "binned",
        stage = "raw",
        mt_pattern = mt_pattern,
        force_rebuild = force_rebuild
      )
      
      object <- add_visium_qc_metrics(
        object = raw_object,
        mt_pattern = mt_pattern
      )
      
      saveRDS(object, object_file)
      return(object)
    }
    
    
    
    if (stage == "filtered") {
      
      message(
        "Building filtered ",
        analysis_mode,
        " object..."
      )
      
      qc_object <- load_visium_object(
        sample_tissue = sample_tissue,
        raw_data_dir = raw_data_dir,
        bin_size = bin_size,
        analysis_mode = "binned",
        stage = "qc",
        mt_pattern = mt_pattern,
        force_rebuild = force_rebuild
      )
      
      object <- filter_visium_object_by_qc(
        object = qc_object,
        min_counts = min_counts,
        max_counts = max_counts,
        min_features = min_features,
        max_features = max_features,
        max_percent_mt = max_percent_mt
      )
      
      saveRDS(object, object_file)
      return(object)
    }
  }
  
  # ----------------------------------------------------------
  # Segmented-cell workflow
  # ----------------------------------------------------------
  if (analysis_mode == "segmented_cells") {
    
    if (stage == "raw") {
      
      message(
        "Building raw ",
        analysis_mode,
        " object..."
      )
      
      # calls function built in
      # cell_segmentation_custom.R
      object <- build_segmented_visium_object(
        sample_id = sample_tissue
      )
      
      saveRDS(object, object_file)
      return(object)
    }
    
    
    
    if (stage == "qc") {
      
      message(
        "Building qc ",
        analysis_mode,
        " object..."
      )
      
      raw_object <- load_visium_object(
        sample_tissue = sample_tissue,
        raw_data_dir = raw_data_dir,
        analysis_mode = "segmented_cells",
        stage = "raw",
        force_rebuild = force_rebuild
      )
      
      object <- add_visium_qc_metrics(
        object = raw_object,
        mt_pattern = mt_pattern
      )
      
      saveRDS(object, object_file)
      return(object)
    }
    
    
    
    if (stage == "filtered") {
      
      message(
        "Building filtered ",
        analysis_mode,
        " object..."
      )
      
      qc_object <- load_visium_object(
        sample_tissue = sample_tissue,
        raw_data_dir = raw_data_dir,
        analysis_mode = "segmented_cells",
        stage = "qc",
        mt_pattern = mt_pattern,
        force_rebuild = force_rebuild
      )
      
      object <- filter_visium_object_by_qc(
        object = qc_object,
        min_counts = min_counts,
        max_counts = max_counts,
        min_features = min_features,
        max_features = max_features,
        max_percent_mt = max_percent_mt
      )
      
      saveRDS(object, object_file)
      return(object)
    }
  }
  #should never happen
  stop("Unknown analysis_mode or stage requested.")
}
    
    



# ============================================================
# 5) Example calls
# ============================================================
#
# Binned data:
#   analysis_mode = "binned"
#   bin_size = 2 / 8 / 16
#
# Segmented cells:
#   analysis_mode = "segmented_cells"
#   bin_size can be left NULL
# ============================================================
# 
# # Example: segmented cells
# object <- load_visium_object(
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   analysis_mode = analysis_mode,
#   stage = "raw",
#   force_rebuild = FALSE
# )
# 
# # Example: segmented-cell QC object
#   qc_object <- load_visium_object(
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   analysis_mode = analysis_mode,
#   stage = "qc",
#   force_rebuild = FALSE
# )
# 
# # Example: segmented-cell filtered object
# filtered_object <- load_visium_object(
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   analysis_mode = analysis_mode,
#   stage = "filtered",
#   min_counts = min_counts,
#   max_counts = max_counts,
#   min_features = min_features,
#   max_features = max_features,
#   max_percent_mt = max_percent_mt,
#   force_rebuild = FALSE
# )


