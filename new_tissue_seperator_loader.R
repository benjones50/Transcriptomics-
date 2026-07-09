

# ============================================================
# TISSUE POLYGON SEPARATOR
# ============================================================
#
# Purpose:
#   Assign each Visium spot to the left or right tissue using
#   a manually drawn polygon.
#
# Workflow
#   1) Read polygon vertices from CSV
#   2) Convert vertices to an sf polygon
#   3) Determine which spots fall inside the polygon
#   4) Add tissue metadata to the Seurat object
#   5) Save a QC plot
#   6) Store polygon information in object@misc
#
# Benjamin Jones
# ============================================================

assign_tissues_by_polygon <- function(
    object,
    sample_tissue,
    polygon_file,
    output_dir = NULL,
    
    tissue_col = "tissue",
    tissue_id_col = "tissue_id",
    inside_col = "inside_polygon",
    
    inside_label = "left",
    outside_label = "right",
    
    show_plot = FALSE,
    verbose = TRUE
) {
  
  # ----------------------------------------------------------
  # Print header
  # ----------------------------------------------------------
  
  if (verbose) {
    
    message(
      "\n",
      "========================================\n",
      "Assign tissues by polygon\n",
      "========================================"
    )
    
  }
  
  # ----------------------------------------------------------
  # Validate object
  # ----------------------------------------------------------
  
  if (!inherits(object, "Seurat")) {
    
    stop(
      "`object` must be a Seurat object."
    )
    
  }
  
  # ----------------------------------------------------------
  # Validate polygon file
  # ----------------------------------------------------------
  
  if (!file.exists(polygon_file)) {
    
    stop(
      "Polygon file does not exist:\n",
      polygon_file
    )
    
  }
  

  # ----------------------------------------------------------
  # Create output directory if needed
  # ----------------------------------------------------------

    if (is.null(output_dir)) {
      
      stop(
        "output_dir must be supplied to save plot"
      )
      
    }
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
  
  
  # ----------------------------------------------------------
  # Extract spatial coordinates
  # ----------------------------------------------------------
  
  coords <- GetTissueCoordinates(object)
  
  coords <- coords[, c("x", "y")]
  
  coords$barcode <- rownames(coords)
  
  pts <- sf::st_as_sf(
    coords,
    coords = c("x", "y")
  )
  
  
  # ----------------------------------------------------------
  # Convert spots to sf points
  # ----------------------------------------------------------
  
  pts <- sf::st_as_sf(
    coords,
    coords = c("x","y")
  )
  
  # ----------------------------------------------------------
  # Read polygon CSV
  # ----------------------------------------------------------
  
  polygon_df <- read.csv(
    polygon_file,
    check.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Locate columns belonging to this sample
  # ----------------------------------------------------------
  polygon_sample <- sub("_prime$", "", sample_tissue)
  
  column_base <- sub("_prime", "", colnames(polygon_df))
  
  polygon_columns <- colnames(polygon_df)[
    column_base %in% c(
      paste0(polygon_sample, "_x"),
      paste0(polygon_sample, "_y")
    )
  ]
  
  if (length(polygon_columns) == 0) {
    
      stop(
        "No polygon columns found for sample: ",
        sample_tissue,
        "\n\nAvailable columns:\n",
        paste(colnames(polygon_df), collapse = "\n")
      )
    
    
  }
  
  if (length(polygon_columns) != 2) {
    
    stop(
      "Expected exactly two polygon columns.\n",
      "Found:\n",
      paste(
        polygon_columns,
        collapse = ", "
      )
    )
    
  }
  
  if (verbose) {
    
    message(
      "Polygon columns:\n",
      "  ",
      paste(
        polygon_columns,
        collapse = ", "
      )
    )
    
  }
  
  # ----------------------------------------------------------
  # Extract vertices
  # ----------------------------------------------------------
  
  vertices <- polygon_df[
    ,
    polygon_columns,
    drop = FALSE
  ]
  
  # Remove NA padding
  
  vertices <- vertices[
    rowSums(
      is.na(vertices)
    ) != ncol(vertices),
    ,
    drop = FALSE
  ]
  
  if (nrow(vertices) < 3) {
    
    stop(
      "Polygon contains fewer than three vertices."
    )
    
  }
  
  # ----------------------------------------------------------
  # Convert to numeric matrix
  # ----------------------------------------------------------
  
  vertex_matrix <- as.matrix(vertices)
  
  storage.mode(vertex_matrix) <- "double"
  
  # ----------------------------------------------------------
  # Close polygon
  # ----------------------------------------------------------
  
  if (
    !all(
      vertex_matrix[1, ] ==
      vertex_matrix[
        nrow(vertex_matrix),
      ]
    )
  ) {
    
    vertex_matrix <- rbind(
      vertex_matrix,
      vertex_matrix[
        1,
        ,
        drop = FALSE
      ]
    )
    
  }
  
  if (verbose) {
    
    message(
      "Vertices used: ",
      nrow(vertex_matrix) - 1
    )
    
  }
  
  # ----------------------------------------------------------
  # Build polygon
  # ----------------------------------------------------------
  
  polygon <- sf::st_polygon(
    list(vertex_matrix)
  )
  
  polygon_sf <- sf::st_sf(
    
    sample_tissue = sample_tissue,
    
    geometry =
      sf::st_sfc(
        polygon
      )
    
  )
  
  # ----------------------------------------------------------
  # Determine which spots lie inside
  # ----------------------------------------------------------
  
  inside <- sf::st_within(
    pts,
    polygon_sf,
    sparse = FALSE
  )[,1]
  
  if (
    length(inside) !=
    nrow(coords)
  ) {
    
    stop(
      "Internal error during polygon assignment."
    )
    
  }
  
  # ----------------------------------------------------------
  # Add metadata
  # ----------------------------------------------------------
  
  object[[inside_col]] <- inside
  
  object[["side"]] <- ifelse(
    inside,
    inside_label,
    outside_label
  )
  
  object$tissue <- paste(
    sample_tissue,
    object$side,
    sep = "_"
  )
  
  
  
  
  # ----------------------------------------------------------
  # Store assignment information in object@misc
  # ----------------------------------------------------------
  
  object@misc$tissue_assignment <- list(
    sample_tissue = sample_tissue,
    polygon_file = polygon_file,
    polygon_columns = polygon_columns,
    polygon = polygon_sf,
    vertices = vertex_matrix,
    inside_label = inside_label,
    outside_label = outside_label
  )
  
  
  
  # ----------------------------------------------------------
  # Assignment summary
  # ----------------------------------------------------------
  
  total_spots <- length(inside)
  
  left_spots <- sum(
    inside,
    na.rm = TRUE
  )
  
  right_spots <- total_spots - left_spots
  
  left_percent <-
    round(
      100 * left_spots / total_spots,
      2
    )
  
  right_percent <-
    round(
      100 * right_spots / total_spots,
      2
    )
  
  
  # ----------------------------------------------------------
  # Print summary
  # ----------------------------------------------------------
  
  if (verbose) {
    
    message(
      "\n",
      "========================================\n",
      "Tissue assignment summary\n",
      "========================================\n",
      "Sample: ", sample_tissue, "\n",
      "Polygon vertices: ", nrow(vertex_matrix) - 1, "\n",
      "Total spots: ", format(total_spots, big.mark=","), "\n",
      "Left tissue: ", format(left_spots, big.mark=","), " (", left_percent, "%)\n",
      "Right tissue: ", format(right_spots, big.mark=","), " (", right_percent, "%)\n",
      "========================================"
    )
    
  }

    
  # ----------------------------------------------------------
  # QC plot
  # ----------------------------------------------------------
  
  plot_df <- data.frame(
    x = coords$x,
    y = coords$y,
    tissue = factor(
      ifelse(
        inside,
        inside_label,
        outside_label
      ),
      levels = c(
        inside_label,
        outside_label
      )
    )
  )
  
  polygon_df <- data.frame(
    x = vertex_matrix[, 1],
    y = vertex_matrix[, 2]
  )
  
  qc_plot <- ggplot() +
    geom_point(
      data = plot_df,
      aes(x, y),
      color = "grey85",
      size = pt_size
    ) +
    geom_point(
      data = plot_df,
      aes(
        x,
        y,
        color = tissue
      ),
      size = pt_size
    ) +
    geom_polygon(
      data = polygon_df,
      aes(x, y),
      fill = "red",
      alpha = 0.08,
      color = "red",
      linewidth = 0.7
    ) +
    coord_equal() +
    theme_classic() +
    labs(
      title = paste0(
        sample_tissue,
        " Tissue Assignment"
      ),
      subtitle = paste0(
        "Left: ",
        format(left_spots, big.mark = ","),
        "   Right: ",
        format(right_spots, big.mark = ",")
      ),
      color = "Tissue"
    )
    
  
  # ----------------------------------------------------------
  # Save QC plot
  # ----------------------------------------------------------
  
  
    save_plot(
      plot_object = qc_plot,
      file_stub = paste0(
        sample_tissue,
        "_tissue_assignment"
      ),
      save_dir = sample_output_dir,
      width = 8,
      height = 8
    )
  
  
  if (show_plot) print(qc_plot)
  
  
  # ----------------------------------------------------------
  # Finished
  # ----------------------------------------------------------
  
  if (verbose) {
    message(
      "\n",
      "Tissue assignment complete.\n"
    )
  }
  return(
    list(
      object = object,
      plot = qc_plot,
      polygon = polygon_sf,
      vertices = vertex_matrix,
      summary = list(
        total_spots = total_spots,
        left_spots = left_spots,
        right_spots = right_spots,
        left_percent = left_percent,
        right_percent = right_percent
      )
    )
  )
  
}


# 
# 
# 
# 
# 
# # Load one filtered object
# object <- load_visium_object(
#   sample_tissue = sample_tissue,
#   raw_data_dir = raw_data_dir,
#   analysis_mode = analysis_mode,
#   stage = "filtered",
#   bin_size = bin_size,
#   mt_pattern = "^mt-",
#   min_counts = min_counts,
#   max_counts = max_counts,
#   min_features = min_features,
#   max_features = max_features,
#   max_percent_mt = max_percent_mt,
#   force_rebuild = TRUE
# )
# 
# 
# polygon_file <- "/projects/nagy_lab_projects/projects_benjones/share/data/tissue_separation_polygon_coords.csv"
# 
# result <- assign_tissues_by_polygon(
#   object = object,
#   sample_tissue = sample_tissue,
#   polygon_file = polygon_file,
#   output_dir = sample_output_dir,
#   show_plot = TRUE,
#   verbose = TRUE
# )
# 
# result
# 
# object <- result$object
# qc_plot <- result$plot
# summary <- result$summary
# 
# 
# qc_plot
# 
# 
# 
# head(object@meta.data)
# 
# table(object$tissue_id)
# str(object@misc$tissue_assignment)
# 
# 
# 
# 
# 
# 
# 
# left_object <- subset(
#   object,
#   subset = tissue == "left"
# )
# 
# right_object <- subset(
#   object,
#   subset = tissue == "right"
# )
# 
# SpatialPlot(
#   left_object,
#   image.alpha = 0,
#   pt.size.factor = pt_size,
#   cols = "forestgreen"
# )
# 
# SpatialPlot(
#   right_object,
#   image.alpha = 0,
#   pt.size.factor = pt_size,
#   cols = "forestgreen"
# )
# 
# 
# 
# 
# 
# 
# 












  
  