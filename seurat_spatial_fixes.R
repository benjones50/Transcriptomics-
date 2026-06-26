#the normal Read10X_Segmentations() doesn't work. 
# missmatch between visium.v2 object 

Read10X_Segmentations_safe <- function(
    image.dir,
    data.dir,
    image.name = "tissue_lowres_image.png",
    assay = "Spatial",
    slice = "slice1"
) {
  
  message(
    "Trying native Seurat Read10X_Segmentations()..."
  )
  
  native_obj <- try(
    Seurat:::Read10X_Segmentations(
      image.dir = image.dir,
      data.dir = data.dir,
      image.name = image.name,
      assay = assay,
      slice = slice
    ),
    silent = TRUE
  )
  
  if (!inherits(native_obj, "try-error")) {
    
    message(
      "✓ Using native Seurat Read10X_Segmentations()"
    )
    
    return(native_obj)
  }
  
  
  
  message(
    "✗ Native Read10X_Segmentations() failed"
  )
  
  message(
    "✓ Using patched VisiumV2 workaround"
  )
  
  
  
  
  
  sf.obj <- Seurat:::Read10X_HD_GeoJson(
    data.dir = data.dir,
    segmentation.type = "cell"
  )
  
  segmentations <- CreateSegmentation(
    sf.obj,
    compact = TRUE
  )
  
  centroids <- CreateCentroids(
    sf.obj,
    nsides = Inf,
    radius = NULL,
    theta = 0
  )
  
  boundaries <- list(
    segmentations = segmentations,
    centroids = centroids
  )
  
  image <- png::readPNG(
    source = file.path(
      image.dir,
      image.name
    )
  )
  
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(
      image.dir,
      "scalefactors_json.json"
    )
  )
  
  
  #only section that is changed
  
  # #this is the default
  # 
  # visium.v2 <- new(
  #   Class = "VisiumV2",
  #   boundaries = boundaries,
  #   assay = assay,
  #   key = key,
  #   image = image,
  #   scale.factors = scale.factors,
  #   coords_x_orientation = "horizontal"
  # )
  
  
  # Seurat code expects:
  #   VisiumV2
  # ├─ image
  # ├─ scale.factors
  # ├─ molecules
  # ├─ boundaries
  # ├─ assay
  # ├─ key
  # └─ coords_x_orientation
  # 
  # loaded class gives:
  #   VisiumV2
  # ├─ image
  # ├─ scale.factors
  # ├─ molecules
  # ├─ boundaries
  # ├─ assay
  # └─ key
  
  #no coords_x_orentation
  
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = boundaries,
    assay = assay,
    key = Key(slice, quiet = TRUE),
    image = image,
    scale.factors = scale.factors
  )
  
  return(visium.v2)
}









Read10X_Image_safe <- function(
    image.dir,
    image.name = "tissue_lowres_image.png",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE,
    image.type = "VisiumV2"
) {
  
  # ==========================================================
  # First try Seurat's native implementation
  # ==========================================================
  #
  # If Seurat has been fixed in a future version, we want to
  # use the official implementation rather than maintaining
  # our own copy forever.
  # ==========================================================
  
  message(
    "Trying native Seurat Read10X_Image()..."
  )
  
  native_obj <- try(
    Seurat:::Read10X_Image(
      image.dir = image.dir,
      image.name = image.name,
      assay = assay,
      slice = slice,
      filter.matrix = filter.matrix,
      image.type = image.type
    ),
    silent = TRUE
  )
  
  if (!inherits(native_obj, "try-error")) {
    
    message(
      "✓ Using native Seurat Read10X_Image()"
    )
    
    return(native_obj)
  }
  
  message(
    "✗ Native Read10X_Image() failed"
  )
  
  message(
    "✓ Using patched VisiumV2 workaround"
  )
  
  
  # ==========================================================
  # Validate image type
  # ==========================================================
  #
  # Seurat currently supports VisiumV1 and VisiumV2 images.
  # Visium HD uses VisiumV2.
  # ==========================================================
  
  image.type <- match.arg(
    image.type,
    choices = c(
      "VisiumV1",
      "VisiumV2"
    )
  )
  
  
  # ==========================================================
  # Locate image file
  # ==========================================================
  #
  # Seurat first looks inside the requested image directory.
  # If the image is missing there, it falls back to the
  # older Space Ranger directory structure.
  # ==========================================================
  
  primary.path <- file.path(
    image.dir,
    image.name
  )
  
  fallback.path <- file.path(
    dirname(dirname(dirname(image.dir))),
    "spatial",
    image.name
  )
  
  
  # ==========================================================
  # Read tissue image
  # ==========================================================
  #
  # This image is used as the background for spatial plots.
  # ==========================================================
  
  image <- tryCatch(
    {
      png::readPNG(primary.path)
    },
    error = function(e) {
      if (file.exists(fallback.path)) {
        png::readPNG(fallback.path)
      } else {
        stop(
          "Neither primary nor fallback image could be read:\n",
          primary.path,
          "\n",
          fallback.path
        )
      }
    }
  )
  
  
  # ==========================================================
  # Read image scaling factors
  # ==========================================================
  #
  # These convert between pixel coordinates on the image and
  # spot/cell coordinates stored by Space Ranger.
  # ==========================================================
  
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(
      image.dir,
      "scalefactors_json.json"
    )
  )
  
  
  # ==========================================================
  # Read spot coordinates
  # ==========================================================
  #
  # For binned Visium data, these are the coordinates of each
  # spatial barcode/spot on the tissue image.
  # ==========================================================
  
  coordinates <- Read10X_Coordinates(
    filename = Sys.glob(
      file.path(
        image.dir,
        "*tissue_positions*"
      )
    ),
    filter.matrix = filter.matrix
  )
  
  
  # ==========================================================
  # Generate Seurat image key
  # ==========================================================
  #
  # The key is used internally by Seurat to identify this
  # image/FOV.
  # ==========================================================
  
  key <- Key(
    slice,
    quiet = TRUE
  )
  
  
  # ==========================================================
  # Legacy VisiumV1 support
  # ==========================================================
  #
  # Included to match Seurat's original implementation.
  # Not expected to be used for Visium HD.
  # ==========================================================
  
  if (image.type == "VisiumV1") {
    
    visium.v1 <- new(
      Class = image.type,
      assay = assay,
      key = key,
      coordinates = coordinates,
      scale.factors = scale.factors,
      image = image
    )
    
    visium.v1@spot.radius <- Radius(
      visium.v1
    )
    
    return(visium.v1)
  }
  
  
  # ==========================================================
  # Build a Field Of View (FOV)
  # ==========================================================
  #
  # Seurat stores spatial coordinates inside an FOV object.
  # This creates centroid locations for each spot.
  # ==========================================================
  
  fov <- CreateFOV(
    coordinates[, c("imagecol", "imagerow")],
    type = "centroids",
    radius = scale.factors[["spot"]],
    assay = assay,
    key = key
  )
  
  
  # ==========================================================
  # Create VisiumV2 image object
  # ==========================================================
  #
  # ONLY DIFFERENCE FROM SEURAT:
  #
  # Native Seurat code attempts:
  #
  #   coords_x_orientation = "horizontal"
  #
  # when constructing VisiumV2.
  #
  # In this environment:
  #
  #   slotNames("VisiumV2")
  #
  # does NOT contain coords_x_orientation, causing:
  #
  #   invalid name for slot of class "VisiumV2"
  #
  # Therefore we construct the same VisiumV2 object but omit
  # the unsupported slot.
  # ==========================================================
  
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = image,
    scale.factors = scale.factors
  )
  
  return(visium.v2)
}






build_binned_visium_object <- function(
    sample_id,
    bin_size = 8
) {
  
  # ==========================================================
  # Determines the folder where originally binned data comes from
  # ==========================================================
  
  bin_folder <- switch(
    as.character(bin_size),
    "2"  = "square_002um",
    "8"  = "square_008um",
    "16" = "square_016um",
    stop("Unsupported bin size")
  )
  
  bin_dir <- file.path(
    raw_data_dir,
    sample_id,
    "outs",
    "binned_outputs",
    bin_folder
  )
  
  # ==========================================================
  # Read matrix
  # ==========================================================
  
  mat <- Read10X_h5(
    file.path(
      bin_dir,
      "filtered_feature_bc_matrix.h5"
    )
  )
  
  # ==========================================================
  # Create Seurat object
  # ==========================================================
  
  object <- CreateSeuratObject(
    counts = mat,
    assay = "Spatial",
    project = sample_id
  )
  
  # ==========================================================
  # Read image using safe constructor
  # ==========================================================
  
  image <- Read10X_Image_safe(
    image.dir = file.path(
      bin_dir,
      "spatial"
    ),
    assay = "Spatial",
    slice = paste0(
      "slice1.",
      sample_id
    )
  )
  
  
  
  # ==========================================================
  # Synchronize image cells and assay cells
  # ==========================================================
  #
  # The image/FOV object and expression matrix must contain
  # the same spots in the same order.
  #
  # Native Load10X_Spatial() handles this internally.
  # Since we are building the object manually, we need to
  # enforce the alignment ourselves.
  #
  # If the image and assay are not aligned:
  #
  #   identical(
  #     Cells(image),
  #     colnames(object)
  #   )
  #
  # returns FALSE and functions such as
  # SpatialFeaturePlot() may fail.
  # ==========================================================
  
  common_cells <- intersect(
    Cells(image),
    colnames(object)
  )
  
  # Keep only spots present in both objects
  
  image <- subset(
    image,
    cells = common_cells
  )
  
  object <- subset(
    object,
    cells = common_cells
  )
  
  # ==========================================================
  # Reorder image cells to match assay cell order
  # ==========================================================
  #
  # Even if both contain the same spots, Seurat requires
  # them to appear in the exact same order.
  # ==========================================================
  
  image <- subset(
    image,
    cells = colnames(object)
  )
  
  # Verify alignment before attaching
  
  stopifnot(
    identical(
      Cells(image),
      colnames(object)
    )
  )
  
  # ==========================================================
  # Attach image to Seurat object
  # ==========================================================
  
  object[[paste0(
    "slice1.",
    sample_id
  )]] <- image
  
  return(object)
}



  
  