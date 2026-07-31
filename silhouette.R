# makes approximate silhouette

# ------------------------------------------------------------
# Internal helper: validate a vector of reduced dimensions
# ------------------------------------------------------------
.silhouette_validate_dims <- function(dims, n_available_dims) {
  dims <- unique(as.integer(dims))
  
  if (length(dims) == 0L || any(is.na(dims))) {
    stop("`dims` must be a non-empty integer vector.")
  }
  
  if (any(dims < 1L)) {
    stop("`dims` must contain positive integers only.")
  }
  
  if (max(dims) > n_available_dims) {
    stop(
      sprintf(
        "`dims` requests dimension %s, but the reduction only has %s dimensions.",
        max(dims),
        n_available_dims
      )
    )
  }
  
  dims
}

# ------------------------------------------------------------
# Function: add_silhouette()
# ------------------------------------------------------------
add_silhouette <- function(
    object,
    reduction,
    dims,
    cluster_column = NULL,
    metadata_column = "silhouette"
) {
  
  if (!reduction %in% Reductions(object)) {
    stop(sprintf("Reduction '%s' was not found in the object.", reduction))
  }
  
  emb <- Embeddings(object, reduction = reduction)
  
  dims <- .silhouette_validate_dims(
    dims = dims,
    n_available_dims = ncol(emb)
  )
  
  if (is.null(cluster_column)) {
    
    cluster_values <- Idents(object)
    cluster_source <- "Idents(object)"
    
  } else {
    
    cluster_values <- object[[cluster_column]][, 1]
    cluster_source <- cluster_column
    
  }
  
  names(cluster_values) <- colnames(object)
  
  # ------------------------------------------------------------
  # Align clusters to the reduction
  # ------------------------------------------------------------
  cells_used <- rownames(emb)
  
  cluster_values <- cluster_values[cells_used]
  
  keep_cells <- !is.na(cluster_values)
  
  cells_used <- cells_used[keep_cells]
  cluster_use <- factor(cluster_values[keep_cells])
  
  if (length(cells_used) < 2L) {
    stop("At least two cells with cluster assignments are required.")
  }
  
  if (nlevels(cluster_use) < 2L) {
    stop("Silhouette requires at least two clusters.")
  }
  
  x <- emb[cells_used, dims, drop = FALSE]
  
  # ------------------------------------------------------------
  # Compute approximate silhouette
  # ------------------------------------------------------------
  approx_df <- as.data.frame(
    bluster::approxSilhouette(
      x = x,
      clusters = cluster_use
    ),
    stringsAsFactors = FALSE
  )
  
  sil_df <- data.frame(
    cell = rownames(approx_df),
    silhouette = as.numeric(approx_df$width),
    stringsAsFactors = FALSE
  )
  
  # ------------------------------------------------------------
  # Store silhouette values
  # ------------------------------------------------------------
  silhouette_values <- rep(NA_real_, ncol(object))
  names(silhouette_values) <- colnames(object)
  
  silhouette_values[sil_df$cell] <- sil_df$silhouette
  
  object[[metadata_column]] <- silhouette_values
  
  # ------------------------------------------------------------
  # Console summary
  # ------------------------------------------------------------
  sil_vals <- sil_df$silhouette

  cat("\n===== Silhouette summary =====\n")
  cat("Metadata column:", metadata_column, "\n")
  cat("Cluster source:", cluster_source, "\n")
  cat("Reduction:", reduction, "\n")
  cat("Dims:", paste(dims, collapse = ", "), "\n")
  cat("Cluster source:", cluster_source, "\n\n")
  cat("Cells used:", format(length(cells_used), big.mark = ","), "\n")
  cat("Clusters used:", nlevels(cluster_use), "\n\n")
  cat(sprintf("Mean: %.4f\n", mean(sil_vals, na.rm = TRUE)))
  cat(sprintf("Median: %.4f\n", stats::median(sil_vals, na.rm = TRUE)))
  cat(sprintf("SD: %.4f\n", stats::sd(sil_vals, na.rm = TRUE)))
  cat(sprintf("Minimum: %.4f\n", min(sil_vals, na.rm = TRUE)))
  cat(sprintf("Maximum: %.4f\n", max(sil_vals, na.rm = TRUE)))
  cat(sprintf("Percent negative: %.2f%%\n", 100 * mean(sil_vals < 0, na.rm = TRUE)))
  cat(sprintf("Percent < 0.25: %.2f%%\n", 100 * mean(sil_vals < 0.25, na.rm = TRUE)))
  cat(sprintf("Percent > 0.75: %.2f%%\n", 100 * mean(sil_vals > 0.75, na.rm = TRUE)))
  
  object
}