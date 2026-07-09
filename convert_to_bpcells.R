# ============================================================
# CONVERT SEURAT OBJECT TO BPCELLS
# ============================================================
#
# Goal:
#   Convert the counts layer of a Seurat v5 object from an
#   in-memory sparse matrix into an on-disk BPCells matrix.
#
#   Stores only the counts layer as a BPCells matrix because it is by far the largest
#   component of the Seurat object and provides the greatest memory savings. The
#   normalized and scaled layers are much smaller.
#
# Why?
#   BPCells stores count matrices on disk rather than in RAM,
#   allowing much larger datasets to be analyzed with lower
#   memory usage.
#
# Behavior:
#   - Creates a BPCells directory next to the Seurat object.
#   - Writes the counts layer to disk.
#   - Replaces the in-memory counts layer with the BPCells
#     disk-backed matrix.
#   - Returns the updated Seurat object.
#
#
#
# ============================================================

convert_to_bpcells <- function(
    object,
    object_file
) {
  
  message(
    "\n==============================\n",
    "Converting object to BPCells\n",
    "=============================="
  )
  
  # ----------------------------------------------------------
  # Determine which assay contains the counts layer.
  # ----------------------------------------------------------
  
  assay_name <- DefaultAssay(object)
  
  message(
    "Assay: ",
    assay_name
  )
  
  
  # ----------------------------------------------------------
  # Build the BPCells directory name from the Seurat object
  # filename.
  # ----------------------------------------------------------
  
  bp_path <- sub(
    "\\.rds$",
    "_bp",
    object_file
  )
  
  
  message(
    "BPCells directory:\n  ",
    bp_path
  )
  
  
  if (dir.exists(bp_path)) {
    
    message(
      "Existing BPCells directory found. Removing..."
    )
    
    unlink(
      bp_path,
      recursive = TRUE,
      force = TRUE
    )
    
  }
  
  
  
  
  # ----------------------------------------------------------
  # Retrieve the counts layer from the Seurat object.
  # ----------------------------------------------------------
  
  message(
    "Loading counts layer..."
  )
  
  counts <- LayerData(
    object,
    assay = assay_name,
    layer = "counts"
  )
  
  # ----------------------------------------------------------
  # Write the counts matrix to disk in BPCells format.
  # ----------------------------------------------------------
  
  message(
    "Writing BPCells matrix to disk..."
  )
  
  write_matrix_dir(
    counts,
    dir = bp_path
  )
  
  # ----------------------------------------------------------
  # Replace the in-memory counts layer with the disk-backed
  # BPCells matrix.
  # ----------------------------------------------------------
  
  message(
    "Replacing counts layer with BPCells matrix..."
  )
  
  LayerData(object, assay = assay_name, layer = "counts") <- open_matrix_dir(bp_path)
  
  message(
    "BPCells conversion complete."
  )
  
  message(
    "Counts layer is now disk-backed.\n"
  )
  
  return(object)
}



