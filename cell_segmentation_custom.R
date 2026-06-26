library(Seurat)
library(BPCells)
library(tidyverse)


#options(future.globals.maxSize = 1e10)

# ============================================================
# Build segmented-cell matrix cache
# ============================================================
#
# uses on-disk BPCells matrices.
# This block creates them the first time, then reuses them.
# ============================================================

bp_root_dir <- file.path(data_dir, "ondisk_bp_cells_mtx_segmented")
dir.create(bp_root_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Helper: build polygon / image object
# ============================================================
#
# Read10X_Segmentations isn't playing nice in this version of seurat. 
# manually copied and edited the function to work myself
#
# Uses native Seurat implementation when possible.
# Falls back to patched VisiumV2 constructor if Seurat's
# implementation fails because of version mismatches.
#
# ============================================================

build_visium_polygons <- function(sample_id) {
  
  image_dir <- file.path(
    raw_data_dir,
    sample_id,
    "outs",
    "segmented_outputs",
    "spatial"
  )
  
  data_dir_out <- file.path(
    raw_data_dir,
    sample_id,
    "outs"
  )
  #uses custom 10x reader that works even with visium v2 object bug
  Read10X_Segmentations_safe(
    image.dir = image_dir,
    data.dir = data_dir_out,
    assay = "Spatial",
    slice = paste0(
      "slice1.",
      sample_id
    )
  )
}







# ============================================================
# build a segmented visium object for one sample. used bc 
# ============================================================

build_segmented_visium_object <- function(sample_id) {
  
  message("\n==============================")
  message("Processing sample: ", sample_id)
  message("==============================")
  
  # ----------------------------------------------------------
  # Load or create BPCells matrix
  # ----------------------------------------------------------
  
  h5_fp <- file.path(
  raw_data_dir,
  sample_id,
  "outs",
  "segmented_outputs",
  "filtered_feature_cell_matrix.h5"
)
  
  bp_dir <- file.path(
    bp_root_dir,
    sample_id
  )
  
  if (!dir.exists(bp_dir)) {
    message("Creating BPCells on-disk matrix for ", sample_id)
    
    h5_tmp <- open_matrix_10x_hdf5(path = h5_fp)
    
    write_matrix_dir(
      mat = h5_tmp,
      dir = bp_dir,
      overwrite = TRUE
    )
  }
  
  matrix_temp <- open_matrix_dir(dir = bp_dir)
  
  
  # ============================================================
  # Convert Ensembl IDs to gene symbols
  # ============================================================
  #
  # The BPCells matrix currently uses Ensembl IDs as row names:
  #
  #   ENSMUSG00000051951
  # We want gene symbols instead:
  #
  #   Xkr4
  #   Rp1
  #   Sox17
  #
  # Space Ranger stores the mapping between:
  #
  #   Ensembl ID <--> Gene Symbol
  #
  # inside:
  #
  #   features.tsv.gz
  #
  # ============================================================
  
  feat_fp <- file.path(
    raw_data_dir,
    sample_id,
    "outs",
    "segmented_outputs",
    "filtered_feature_cell_matrix",
    "features.tsv.gz"
  )
  
  # Read the feature annotation table.
  #
  # X1 = Ensembl ID
  # X2 = Gene Symbol
  # X3 = Feature Type
  #
  # Example:
  #
  # ENSMUSG00000051951   Xkr4    Gene Expression
  #
  feature_temp <- readr::read_delim(
    feat_fp,
    delim = "\t",
    col_names = FALSE
  ) %>%
    dplyr::filter(
      X3 == "Gene Expression"
    ) %>%
    dplyr::select(
      ensembl = X1,
      symbol = X2
    )
  
  # ============================================================
  # Match every matrix row to the annotation table
  # ============================================================
  #
  # match() returns the row number in feature_temp corresponding
  # to each gene in the matrix.
  #
  # Example:
  #
  # matrix gene               match result
  # -----------------------   ------------
  # ENSMUSG00000051951        1
  # ENSMUSG00000025900        2
  #
  # These indices let us pull out the gene symbols.
  # ============================================================
  
  idx <- match(
    rownames(matrix_temp),
    feature_temp$ensembl
  )
  
  # Replace Ensembl IDs with gene symbols.
  #
  # Example:
  #
  # ENSMUSG00000051951 -> Xkr4
  # ENSMUSG00000025900 -> Rp1
  #
  gene_symbol <- feature_temp$symbol[idx]
  
  # ============================================================
  # Handle missing gene symbols
  # ============================================================
  #
  # Occasionally a gene may not have an assigned symbol.
  #
  # Instead of creating NA rownames, keep the original
  # Ensembl ID.
  #
  # Example:
  #
  # ENSMUSG00000123456 -> ENSMUSG00000123456
  #
  # rather than:
  #
  # ENSMUSG00000123456 -> NA
  # ============================================================
  
  gene_symbol[
    is.na(gene_symbol) |
      gene_symbol == ""
  ] <- rownames(matrix_temp)[
    is.na(gene_symbol) |
      gene_symbol == ""
  ]
  
  # ============================================================
  # Handle duplicated gene symbols
  # ============================================================
  #
  # Multiple Ensembl IDs can map to the same gene symbol.
  #
  # Example:
  #
  # ENSMUSG_A -> GeneX
  # ENSMUSG_B -> GeneX
  #
  # Seurat requires unique feature names.
  #
  # ============================================================
  
  f <- factor(
    gene_symbol,
    levels = unique(gene_symbol)
  )
  
  # ============================================================
  # Build a sparse mapping matrix
  # ============================================================
  #
  # Multiplying this by the count matrix sums all rows
  # belonging to the same gene symbol.
  #
  # ============================================================
  
  map_matrix <- Matrix::sparseMatrix(
    i = as.integer(f),
    j = seq_along(f),
    x = 1,
    dims = c(
      nlevels(f),
      length(f)
    ),
    dimnames = list(
      levels(f),
      rownames(matrix_temp)
    )
  )
  
  # ============================================================
  # Collapse duplicated genes
  # ============================================================
  #
  # Before:
  #
  # ENSMUSG_A
  # ENSMUSG_B
  #
  # After:
  #
  # GeneX
  #
  #
  # Final result:
  #
  # rownames(matrix_temp) now contain gene symbols instead
  # of Ensembl IDs.
  #
  # ============================================================
  
  matrix_temp <- map_matrix %*% matrix_temp
  
  # ----------------------------------------------------------
  # Create Seurat object
  # ----------------------------------------------------------
  
  seurat_temp <- CreateSeuratObject(
    counts = matrix_temp,
    assay = "Spatial",
    project = sample_id
  )
  
  # ----------------------------------------------------------
  # Build polygon/image object
  # ----------------------------------------------------------
  
  #uses helper function from beginning
  polygon_object <- build_visium_polygons(sample_id)
  
  # ----------------------------------------------------------
  # Remove cells that do not exist in both objects, polygon and seurat
  # ----------------------------------------------------------
  
  common_temp <- intersect(
    Cells(polygon_object),
    colnames(seurat_temp)
  )
  
  #removal based on whats in common
  polygon_object <- subset(
    polygon_object,
    cells = common_temp
  )
  
  seurat_temp <- subset(
    seurat_temp,
    cells = common_temp
  )
  
  # ----------------------------------------------------------
  # Make sure cell order is identical
  # ----------------------------------------------------------
  
  if (!identical(
    Cells(polygon_object),
    colnames(seurat_temp)
  )) {
    stop(
      "Cell IDs in Seurat object and polygon object are not the same or not in the same order for sample ",
      sample_id
    )
  }
  
  # ----------------------------------------------------------
  # Attach polygon/image object to Seurat object
  # ----------------------------------------------------------
  
  seurat_temp[[paste0("slice1.", sample_id)]] <- polygon_object
  
  # # ----------------------------------------------------------
  # # Add mitochondrial percentage
  # # ----------------------------------------------------------
  # 
  # seurat_temp[["percent.mt"]] <- PercentageFeatureSet(
  #   seurat_temp,
  #   pattern = "^mt-"
  # )
  

  
  return(seurat_temp)
}









