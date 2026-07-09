# mega/
#     2um/
#         SCT/
#             QC_15_1500_15_1250_MT20/
# 
#                 object.rds
# 
#                 metadata/
#                     run_info.yaml
#                     session_info.txt
#                     object_summary.csv
# 
#                 dimension_reduction/
#                     pca/
#                         dims50/
#                             object.rds          
#                             elbow.pdf
#                             loadings.pdf
# 
#                 clustering/
#                     dims30_res0.8/
#                         object.rds              
#                         umap.pdf
#                         cluster_markers.csv
#                         cluster_summary.csv
# 
#                 annotation/
#                     module_scores/
#                     cell_types/
# 
#                 differential_expression/
# 
#                 figures/
# 
#                 exports/
# 
# 




# ------------------------------------------------------------
# Build the root output directory for a merged ("mega") analysis
#
# This function determines where the root directory for a merged
# Seurat analysis should live on disk.
#
# Every downstream helper function (PCA, clustering, annotation,
# exports, etc.) should build upon the directory returned here.
#
# Example:
#
# mega/
#   2um/
#     SCT/
#       QC_15_1500_15_1250_MT20/
#
# This function is ONLY responsible for determining that root
# directory and ensuring it exists.
#
# It does NOT:
#
#   * save objects
#   * load objects
#   * create analysis subdirectories
#   * know anything about PCA, clustering, Harmony, etc.
#
# It simply returns the canonical analysis directory.
# ------------------------------------------------------------
build_mega_object_dir <- function(
    analysis_mode = c(
        "binned",
        "segmented_cells"
    ),
    normalization = c(
      "lognorm","sct"
    ),
    bin_size = NULL,
    min_counts,
    max_counts,
    min_features,
    max_features,
    max_percent_mt
) {

    # ==========================================================
    # Validate arguments
    # ==========================================================

    analysis_mode <- match.arg(analysis_mode)
    normalization <- match.arg(normalization)

    # Binned analyses require a bin size
    if (analysis_mode == "binned") {

        if (is.null(bin_size)) {
            stop(
                "bin_size must be supplied when ",
                "analysis_mode = 'binned'.",
                call. = FALSE
            )
        }

    }

    # Segmented-cell analyses should not use bin_size
    if (
        analysis_mode == "segmented_cells" &&
        !is.null(bin_size)
    ) {
        warning(
            "bin_size is ignored when ",
            "analysis_mode = 'segmented_cells'."
        )
    }

    # QC thresholds are required because every downstream
    # analysis is uniquely defined by these filtering values.
    required_args <- list(
        min_counts     = min_counts,
        max_counts     = max_counts,
        min_features   = min_features,
        max_features   = max_features,
        max_percent_mt = max_percent_mt
    )

    missing_args <- names(required_args)[
        vapply(required_args, is.null, logical(1))
    ]

    if (length(missing_args) > 0) {

        stop(
            "Missing required QC parameter(s): ",
            paste(missing_args, collapse = ", "),
            call. = FALSE
        )

    }

    # Verify the project output directory exists
    if (!dir.exists(output_dir)) {

        stop(
            "Output directory does not exist:\n",
            output_dir,
            call. = FALSE
        )

    }

    # ==========================================================
    # Build directory hierarchy
    # ==========================================================

    if (analysis_mode == "binned") {

        analysis_dir <- file.path(
            output_dir,
            "mega",
            paste0(bin_size, "um"),
            paste0(
                "QC_",
                min_counts, "_",
                max_counts, "_",
                min_features, "_",
                max_features, "_MT",
                max_percent_mt
            ),
            normalization
        )

    } else {

        analysis_dir <- file.path(
            output_dir,
            "mega",
            "segmented_cells",
            paste0(
                "QC_",
                min_counts, "_",
                max_counts, "_",
                min_features, "_",
                max_features, "_MT",
                max_percent_mt
            ),
            normalization,
        )

    }

    # ==========================================================
    # Create the directory if necessary
    # ==========================================================

    dir.create(
        analysis_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    # ==========================================================
    # Return the analysis directory
    # ==========================================================

    return(analysis_dir)

}








# ------------------------------------------------------------
# Get an analysis subdirectory
#
# Returns (and creates if necessary) a named analysis
# subdirectory beneath a mega analysis root directory.
#
# This is the generic helper used by all directory accessors
# (metadata, exports, figures, clustering, etc.).
#
# Example
# -------
# get_nested_subdir(
#     root_dir,
#     "metadata"
# )
#
# returns
#
# root_dir/
#     metadata/
#
# Responsibilities
# ----------------
# * validate root_dir
# * validate subdirectory_name
# * create the directory if necessary
# * return the full path
#
# This function knows NOTHING about PCA, clustering,
# metadata, figures, exports, etc.
# ------------------------------------------------------------
# ------------------------------------------------------------
# Get a nested analysis subdirectory
#
# Returns (and creates if necessary) a nested directory beneath
# an existing parent directory.
#
# This helper is intended for deeper directory hierarchies such
# as:
#
# dimension_reduction/
#     pca/
#         dims50/
#
# or
#
# dimension_reduction/
#     harmony/
#         dims50_theta2/
#
# Example
# -------
# get_nested_subdir(
#     parent_dir,
#     "pca",
#     "dims50"
# )
#
# returns
#
# parent_dir/
#     pca/
#         dims50/
#
# Responsibilities
# ----------------
# * validate parent_dir
# * validate subdirectory names
# * build the nested directory path
# * create the directory if necessary
# * return the full directory path
#
# This function knows NOTHING about PCA, Harmony, clustering,
# annotation, Seurat, or any analysis.
# ------------------------------------------------------------

get_nested_subdir <- function(
    parent_dir,
    ...
) {
  
  # ==========================================================
  # Validate parent directory
  # ==========================================================
  
  if (missing(parent_dir) || is.null(parent_dir)) {
    stop(
      "parent_dir must be supplied.",
      call. = FALSE
    )
  }
  
  if (!is.character(parent_dir) || length(parent_dir) != 1) {
    stop(
      "parent_dir must be a single character string.",
      call. = FALSE
    )
  }
  
  if (!dir.exists(parent_dir)) {
    stop(
      "The supplied parent_dir does not exist:\n",
      parent_dir,
      call. = FALSE
    )
  }
  
  # ==========================================================
  # Collect nested directory names
  # ==========================================================
  
  subdirs <- c(...)
  
  if (length(subdirs) == 0) {
    stop(
      "At least one subdirectory name must be supplied.",
      call. = FALSE
    )
  }
  
  if (!all(vapply(subdirs, is.character, logical(1)))) {
    stop(
      "All subdirectory names must be character strings.",
      call. = FALSE
    )
  }
  
  if (!all(vapply(subdirs, length, integer(1)) == 1)) {
    stop(
      "Each subdirectory name must be a single character string.",
      call. = FALSE
    )
  }
  
  if (any(nchar(subdirs) == 0)) {
    stop(
      "Subdirectory names cannot be empty strings.",
      call. = FALSE
    )
  }
  
  # ==========================================================
  # Build nested directory
  # ==========================================================
  
  nested_dir <- do.call(
    file.path,
    c(list(parent_dir), as.list(subdirs))
  )
  
  # ==========================================================
  # Create directory if needed
  # ==========================================================
  
  dir.create(
    nested_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # ==========================================================
  # Return directory
  # ==========================================================
  
  return(nested_dir)
  
}







get_metadata_dir <- function(root_dir) {
  
  get_nested_subdir(
    parent_dir = root_dir,
    "metadata"
  )
  
}

get_figures_dir <- function(root_dir) {
  
  get_nested_subdir(
    parent_dir = root_dir,
    "figures"
  )
  
}
get_exports_dir <- function(root_dir) {
  
  get_nested_subdir(
    parent_dir = root_dir,
    "exports"
  )
  
}











#next layer:

get_pca_dir <- function(
    root_dir,
    dims
) {
  
  
  
  get_nested_subdir(
    root_dir,
    "pca",
    paste0("dims", dims)
  )
  
}





# ------------------------------------------------------------
# Get the clustering analysis directory
#
# Returns (and creates if necessary) the directory used to
# store a sketch clustering analysis.
#
# Example
# -------
# clustering/
#     dims30_res3/
# ------------------------------------------------------------
get_clustering_dir <- function(
    root_dir,
    dims,
    resolution
) {
  
  get_nested_subdir(
    parent_dir = root_dir,
    "clustering",
    paste0(
      "dims",
      max(dims),
      "_res",
      resolution
    )
  )
  
}


# ------------------------------------------------------------
# Get the sketch projection analysis directory
#
# Returns (and creates if necessary) the directory used to
# store a projected sketch analysis.
#
# Example
# -------
# projection/
#     dims30/
# ------------------------------------------------------------
get_projection_dir <- function(
    root_dir,
    dims
) {
  
  get_nested_subdir(
    parent_dir = root_dir,
    "projection",
    paste0(
      "dims",
      max(dims)
    )
  )
  
}



# ------------------------------------------------------------
# Get the Harmony analysis directory
#
# Returns (and creates if necessary) the directory used to
# store a Harmony integration run.
#
# Example
# -------
# dimension_reduction/
#     harmony/
#         dims50_theta2/
# ------------------------------------------------------------
get_harmony_dir <- function(
    root_dir,
    dims,
    theta
) {
  
  
  get_nested_subdir(
    root_dir,
    "harmony",
    paste0(
      "dims",
      dims,
      "_theta",
      theta
    )
  )
  
}



# ------------------------------------------------------------
# Get the UMAP analysis directory
#
# Returns (and creates if necessary) the directory used to
# store a UMAP computed from a specified number of dimensions.
#
# Example
# -------
# dimension_reduction/
#     umap/
#         dims30/
# ------------------------------------------------------------
get_umap_dir <- function(
    root_dir,
    dims
) {
  
  get_nested_subdir(
    root_dir,
    "umap",
    paste0("dims", max(dims))
  )
  
}




# ------------------------------------------------------------
# Get the Spatially Variable Gene (SVI) analysis directory
#
# Returns (and creates if necessary) the directory used to
# store results from a spatially variable gene analysis.
#
# Supported methods include:
#   - moran
#   - spark
#   - hotspot
#
# Example
# -------
# dimension_reduction/
#     svi/
#         moran/
# ------------------------------------------------------------
get_svi_dir <- function(
    root_dir,
    method = c(
      "moran",
      "spark",
      "hotspot"
    )
) {
  
  method <- match.arg(method)
  
  get_nested_subdir(
    root_dir,
    "svi",
    method
  )
  
}





# ------------------------------------------------------------
# Build the path to a figure
#
# Returns the canonical path for a figure associated with an
# analysis directory.
#
# If necessary, the figures directory is created.
#
# Example
# -------
# build_figure_path(
#     directory = pca_dir,
#     filename = "elbow.pdf"
# )
#
# returns
#
# pca/
#     dims50/
#         figures/
#             elbow.pdf
#
# Responsibilities
# ----------------
# * validate filename
# * ensure the figures directory exists
# * return the complete figure path
#
# This function does NOT:
#
# * save the figure
# * load the figure
# * know anything about PCA, Harmony, clustering, etc.
# ------------------------------------------------------------

#basically paste0 when actually making a file name, but makes sure dir is there
build_figure_path <- function(
    directory,
    filename = NULL
) {
  
  # ==========================================================
  # Validate filename
  # ==========================================================
  
  
  if (!is.character(filename) || length(filename) != 1) {
    stop(
      "filename must be a single character string.",
      call. = FALSE
    )
  }
  
  if (nchar(filename) == 0) {
    stop(
      "filename cannot be an empty string.",
      call. = FALSE
    )
  }
  
  # ==========================================================
  # Get (or create) the figures directory
  # ==========================================================
  
  figures_dir <- get_nested_subdir(
    directory,
    "figures"
  )
  
  # ==========================================================
  # Return the full figure path
  # ==========================================================
  
  file.path(
    figures_dir,
    filename
  )
  
}





# ------------------------------------------------------------
# Build the canonical object path for an analysis directory
#
# Given an analysis directory, return the location of the
# cached Seurat object.
#
# Example
# -------
# clustering/
#     dims30_res0.8/
#
# becomes
#
# clustering/
#     dims30_res0.8/
#         object.rds
#
# This function NEVER creates directories, loads objects,
# or saves objects.
# ------------------------------------------------------------
build_object_path_new <- function(directory) {
  
  if (missing(directory) || is.null(directory)) {
    stop(
      "directory must be supplied.",
      call. = FALSE
    )
  }
  
  if (!is.character(directory) || length(directory) != 1) {
    stop(
      "directory must be a single character string.",
      call. = FALSE
    )
  }
  
  file.path(
    directory,
    "object.rds"
  )
  
}





# ------------------------------------------------------------
# Save a cached Seurat object
#
# Saves a Seurat object to the canonical cache location for an
# analysis directory.
#
# This function NEVER determines where objects belong.
# It only saves to the location returned by build_object_path_new().
# ------------------------------------------------------------
save_mega_object <- function(
    object,
    directory
) {
  
  object_path <- build_object_path_new(directory)
  
  message("----------------------------------------")
  message("Saving cached Seurat object")
  message(object_path)
  
  tryCatch(
    
    {
      
      saveRDS(
        object = object,
        file = object_path
      )
      
      message("Finished saving object.")
      
    },
    
    error = function(e) {
      
      stop(
        "Unable to save Seurat object.\n",
        "Path: ", object_path, "\n",
        e$message,
        call. = FALSE
      )
      
    }
    
  )
  
  invisible(object_path)
  
}




# ------------------------------------------------------------
# Load a cached Seurat object
#
# Loads a Seurat object from the canonical cache location for
# an analysis directory.
#
# This function NEVER determines where objects belong.
# It only loads from the location returned by
# build_object_path_new().
# ------------------------------------------------------------
load_mega_object <- function(directory) {
  
  object_path <- build_object_path_new(directory)
  
  if (!file.exists(object_path)) {
    
    stop(
      "Cached object not found:\n",
      object_path,
      call. = FALSE
    )
    
  }
  
  message("----------------------------------------")
  message("Loading cached Seurat object")
  message(object_path)
  
  object <- tryCatch(
    
    {
      
      readRDS(object_path)
      
    },
    
    error = function(e) {
      
      stop(
        "Unable to load cached Seurat object.\n",
        "Path: ", object_path, "\n",
        e$message,
        call. = FALSE
      )
      
    }
    
  )
  
  message("Finished loading object.")
  
  object
  
}








# ------------------------------------------------------------
# Check whether a cached object exists
#
# Determines whether an analysis directory contains a cached
# Seurat object at the canonical cache location.
#
# Returns
# -------
# TRUE  : object.rds exists
# FALSE : object.rds does not exist
#
# This function NEVER loads the object.
# It simply checks whether it exists.
# ------------------------------------------------------------
cached_object_exists <- function(directory) {
  
  object_path <- build_object_path_new(directory)
  
  file.exists(object_path)
  
}






# ------------------------------------------------------------
# Load a cached object if it exists, otherwise build it.
#
# This is the orchestration helper for stage-based pipelines.
#
# Typical workflow
# ----------------
# 1) If a cached object exists, load and return it.
# 2) Otherwise, call build_function().
# 3) Save the newly created object.
# 4) Return the newly built object.
#
# Parameters
# ----------
# directory:
#   Analysis directory for this pipeline stage.
#
# build_function:
#   A zero-argument function that performs the analysis and
#   returns the completed object.
#
#   Example:
#
#   load_or_build_object(
#     directory = pca_dir,
#     build_function = function() {
#       run_sketch_pca(
#         object = mega_obj,
#         ncells = 50000,
#         dims = 1:50
#       )
#     }
#   )
#
# force_rebuild:
#   If TRUE, always rebuild and overwrite any existing cache.
#
# Notes
# -----
# * Assumes the analysis directory already exists.
# * Does not create directories.
# * Does not know anything about PCA, clustering, UMAP, etc.
# * Simply decides whether to load or build.
# ------------------------------------------------------------
load_or_build_object <- function(
    directory,
    build_function,
    force_rebuild = FALSE
) {
  
  # ==========================================================
  # Validate analysis directory
  # ==========================================================
  
  if (!dir.exists(directory)) {
    stop(
      "Analysis directory does not exist:\n",
      directory,
      call. = FALSE
    )
  }
  
  if (!is.function(build_function)) {
    stop(
      "build_function must be a function.",
      call. = FALSE
    )
  }
  
  object_path <- build_object_path_new(directory)
  
  # ==========================================================
  # Load existing cached object
  # ==========================================================
  
  if (!force_rebuild && cached_object_exists(directory)) {
    
    message("----------------------------------------")
    message("Cached object found")
    message("Loading:")
    message("  ", object_path)
    
    return(
      load_mega_object(directory)
    )
    
  }
  
  # ==========================================================
  # Build a new object
  # ==========================================================
  
  if (force_rebuild && cached_object_exists(directory)) {
    
    message("----------------------------------------")
    message("force_rebuild = TRUE")
    message("Rebuilding cached object")
    message("Overwriting:")
    message("  ", object_path)
    
  } else {
    
    message("----------------------------------------")
    message("Cached object not found")
    message("Building new object")
    message("Saving to:")
    message("  ", object_path)
    
  }
  
  object <- tryCatch(
    
    build_function(),
    
    error = function(e) {
      
      stop(
        "Unable to build object.\n\n",
        "Directory:\n",
        directory,
        "\n\n",
        "Reason:\n",
        e$message,
        call. = FALSE
      )
      
    }
    
  )
  
  # ==========================================================
  # Save newly built object
  # ==========================================================
  
  tryCatch(
    
    save_mega_object(
      object = object,
      directory = directory
    ),
    
    error = function(e) {
      
      stop(
        "Object was built successfully but could not be saved.\n\n",
        "Directory:\n",
        directory,
        "\n\n",
        "Reason:\n",
        e$message,
        call. = FALSE
      )
      
    }
    
  )
  
  message("Build complete.")
  
  return(object)
  
}




