#used to change the data in mega_obj@meta.data$tissue to be more accurate bio

# ------------------------------------------------------------
# Rename metadata values using a lookup table
#
# Replaces values in a metadata column according to a named
# lookup vector.
#
# Parameters
# ----------
# object:
#     Seurat object.
#
# column:
#     Metadata column to read from.
#
# lookup:
#     Named character vector where:
#
#         names(lookup) = existing values
#         lookup values = replacement values
#
# new_column:
#     Metadata column to write the renamed values to.
#
# Notes
# -----
# If new_column is identical to column, the original metadata
# column is overwritten. Otherwise, a new metadata column is
# created and the original is preserved.
#
# Example
# -------
# lookup <- c(
#     "ABC" = "Woohoo",
#     "XYZ" = "Banana"
# )
#
# object <- rename_metadata_values(
#     object = object,
#     column = "tissue",
#     lookup = lookup,
#     new_column = "tissue_name"
# )
# ------------------------------------------------------------
rename_metadata_values <- function(
    object,
    column,
    lookup,
    new_column = column
) {
  
  # ==========================================================
  # Validate inputs
  # ==========================================================

  if (!column %in% colnames(object[[]])) {
    stop(
      "Metadata column '",
      column,
      "' not found.",
      call. = FALSE
    )
  }
  
  
  # ==========================================================
  # Rename matching values
  # ==========================================================
  
  values <- object[[column]][, 1]
  
  matches <- values %in% names(lookup)
  
  values[matches] <- lookup[values[matches]]
  
  # ==========================================================
  # Save to requested metadata column
  # ==========================================================
  
  object[[new_column]] <- values
  
  return(object)
  
}





#splits info in a column into a list of other columns based on a delimiter

split_metadata_column <- function(
    object,
    column,
    into,
    delimiter = "__"
) {
  
  values <- object[[column]][, 1]
  
  split_values <- strsplit(
    values,
    split = delimiter,
    fixed = TRUE
  )
  
  lengths <- lengths(split_values)
  
  if (!all(lengths == length(into))) {
    stop(
      "Not every value split into ",
      length(into),
      " fields.",
      call. = FALSE
    )
  }
  
  split_matrix <- do.call(
    rbind,
    split_values
  )
  
  split_df <- as.data.frame(
    split_matrix,
    stringsAsFactors = FALSE
  )
  
  colnames(split_df) <- into
  
  object <- AddMetaData(
    object,
    metadata = split_df
  )
  
  return(object)
  
}
