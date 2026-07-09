# ============================================================
# Load segmented-cell count matrix
# ============================================================
#
# Space Ranger produced a cell-level count matrix where:
#   rows    = genes
#   columns = segmented cells
#
# These are NOT 2um/8um/16um bins.
# Each column is an actual segmented cell from Space Ranger.
# ============================================================

mat <- Read10X_h5(
  file.path(
    raw_data_dir,
    sample_tissue,
    "outs",
    "segmented_outputs",
    "filtered_feature_cell_matrix.h5"
  )
)

# Check matrix dimensions
dim(mat)


# ============================================================
# Create Seurat object from segmented-cell matrix
# ============================================================
#
# This gives us the gene-expression side of the object.
# At this point, the object has counts, but no spatial image yet.
# ============================================================

object <- CreateSeuratObject(
  counts = mat,
  assay = "Spatial",
  project = sample_name
)

object


# ============================================================
# Build the polygon / image object manually
# ============================================================
#
# IMPORTANT:
# Seurat:::Read10X_Segmentations() exists in Seurat 5.5.0,
# but in this Seurat / SeuratObject combination it still
# errors when it tries to create a VisiumV2 object with the
# coords_x_orientation slot.
#
# The workaround is to build the same object manually, but
# without that unsupported slot.
# ============================================================

# Read the cell segmentation GeoJSON from Space Ranger output
sf.obj <- Seurat:::Read10X_HD_GeoJson(
  data.dir = file.path(
    raw_data_dir,
    sample_tissue,
    "outs"
  ),
  segmentation.type = "cell"
)

# Convert the GeoJSON polygons into Seurat segmentation objects
segmentations <- CreateSegmentation(
  sf.obj,
  compact = TRUE
)

# Create centroid coordinates for the segmented cells
centroids <- CreateCentroids(
  sf.obj,
  nsides = Inf,
  radius = NULL,
  theta = 0
)

# Store segmentation boundaries and centroids together
boundaries <- list(
  segmentations = segmentations,
  centroids = centroids
)

# Read the tissue image used for spatial plotting
image <- png::readPNG(
  source = file.path(
    raw_data_dir,
    sample_tissue,
    "outs",
    "segmented_outputs",
    "spatial",
    "tissue_lowres_image.png"
  )
)

# Read the spatial scale factors from Space Ranger output
scale.factors <- Read10X_ScaleFactors(
  filename = file.path(
    raw_data_dir,
    sample_tissue,
    "outs",
    "segmented_outputs",
    "spatial",
    "scalefactors_json.json"
  )
)

# Create the VisiumV2 object manually
# This works around the coords_x_orientation slot problem
pol_obj <- new(
  Class = "VisiumV2",
  boundaries = boundaries,
  assay = "Spatial",
  key = Key(
    paste0("slice1.", sample_tissue),
    quiet = TRUE
  ),
  image = image,
  scale.factors = scale.factors
)

pol_obj


# ============================================================
# Update the polygon object before attaching it
# ============================================================
#
# This helps Seurat bring the image object closer to the newer
# internal structure before it is attached to the Seurat object.
# ============================================================

pol_obj <- UpdateSeuratObject(pol_obj)


# ============================================================
# Match cells between counts and polygons
# ============================================================
#
# Keep only cells present in BOTH objects.
# This avoids mismatches later.
# ============================================================

common_cells <- intersect(
  Cells(pol_obj),
  colnames(object)
)

length(common_cells)


# ============================================================
# Remove unmatched cells
# ============================================================
#
# Now both objects contain the same cells.
# ============================================================

pol_obj <- subset(
  pol_obj,
  cells = common_cells
)

object <- subset(
  object,
  cells = common_cells
)


# ============================================================
# Verify ordering
# ============================================================
#
# Cell order must be identical in both objects.
# ============================================================

if (!identical(
  Cells(pol_obj),
  colnames(object)
)) {
  stop(
    "Cell order mismatch between polygons and counts."
  )
}


# ============================================================
# Attach polygons to Seurat object
# ============================================================
#
# This is the key step that populates object@images.
# ============================================================

object[[paste0(
  "slice1.",
  sample_tissue
)]] <- pol_obj


# ============================================================
# Update the full Seurat object
# ============================================================
#
# This is needed because the image was manually constructed.
# ============================================================

object <- UpdateSeuratObject(object)


# ============================================================
# Confirm the image/segmentation was attached
# ============================================================

Images(object)


# ============================================================
# Add mitochondrial percentage
# ============================================================

object[["percent.mt"]] <- PercentageFeatureSet(
  object,
  pattern = "^mt-"
)


# ============================================================
# Try the spatial plot
# ============================================================
#
# NOTE:
# In this package combination, SpatialFeaturePlot() may still
# complain about object/image alignment even after UpdateSeuratObject().
# If that happens, the segmentation attachment is still usable,
# but Seurat's spatial plotting validator is being strict.
# ============================================================

spatial_test <- try(
  SpatialFeaturePlot(
    object,
    features = "percent.mt",
    images = paste0("slice1.", sample_tissue)
  ),
  silent = TRUE
)

if (inherits(spatial_test, "try-error")) {
  message(
    "SpatialFeaturePlot still fails for this object/image combo. ",
    "The object is attached, but Seurat's validator is still unhappy."
  )
  print(
    ImageDimPlot(object)
  )
} else {
  print(spatial_test)
}


# ============================================================
# Alternative visualization
# ============================================================
#
# This should still show the segmentation boundaries.
# ============================================================

ImageDimPlot(object)


# ============================================================
# Save object
# ============================================================

# saveRDS(
#   object,
#   file.path(
#     sample_output_dir,
#     paste0(
#       sample_name,
#       "_segmented_cells.rds"
#     )
#   )
# )