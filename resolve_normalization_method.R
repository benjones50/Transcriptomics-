# ==========================================================
# Internal helper: resolve the normalization method
#
# ProjectData() needs to know whether the original assay
# came from SCTransform or standard log normalization.
#
# This is determined from the assay class so the projection
# stage can work across both workflows.
# ==========================================================
resolve_normalization_method <- function(object, full_assay) {
  
  full_assay_object <- object[[full_assay]]
  
  if (inherits(full_assay_object, "SCTAssay")) {
    
    return("SCT")
    
  } else if (inherits(full_assay_object, "Assay5")) {
    
    return("LogNormalize")
    
  } else {
    
    stop(
      "Unsupported assay class for full_assay '",
      full_assay,
      "': ",
      class(full_assay_object)[1],
      call. = FALSE
    )
    
  }
}



# ==========================================================
# Internal helper: resolve the full assay
#
# Problem:
# Downstream projection needs the original assay, but the
# default assay may already be set to the sketch assay.
#
# If the user does not explicitly specify `full_assay`,
# automatically choose the first assay that is NOT the sketch.
#
# This works for:
#   Spatial.016um
#   Spatial.008um
#   Spatial
#   RNA
#   SCT
#   segmented-cell assays
# ==========================================================
resolve_full_assay <- function(object, sketch_assay, full_assay = NULL) {
  
  if (is.null(full_assay)) {
    
    full_assay <- setdiff(
      Assays(object),
      sketch_assay
    )[1]
  }
  
  
  if (is.na(full_assay)) {
    stop(
      "Could not determine the full assay automatically.",
      call. = FALSE
    )
  }
  
  
  if (identical(full_assay, sketch_assay)) {
    stop(
      "full_assay cannot be the same as sketch_assay ('", sketch_assay, "').",
      call. = FALSE
    )
  }
  
  return(full_assay)
}