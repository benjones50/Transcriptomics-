library(Seurat)
library(sf)

#gets directories, bin size etc
source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")


#TODO could wrap all of this into a couple functions? load gates, apply them, generate new gates, view gates etc





#loads seurat object
source(file.path(code_dir, "seurat_object_loader.R"))

#does this first time, otherwise saves and loads rds file
# object <- Load10X_Spatial(
#   data.dir = paste0(raw_data_dir, sample_tissue, "/outs"),
#   bin.size = bin_size
# )





# ------------------------------------------------------------
# Coordinates
# ------------------------------------------------------------
coords <- GetTissueCoordinates(object)

xy <- data.frame(
  x = coords$x,
  y = coords$y
)

pts <- st_as_sf(
  xy,
  coords = c("x", "y")
)

# ------------------------------------------------------------
# Polygon file
# Save/load your tissue polygons here instead of manual selection
# ------------------------------------------------------------
poly_file <- file.path(data_dir, sample_tissue, "tissue_polygons.rds")

# Start with all spots unassigned
object$tissue <- "Unassigned"

# ------------------------------------------------------------
# Define polygons if no saved polygon file exists
# Edit these coordinates for your sample
# ------------------------------------------------------------
if (file.exists(poly_file)) {
  polygons <- readRDS(poly_file)
} else {
  # Example polygon for left tissue
  poly_left <- st_polygon(list(
    matrix(
      c(
        0, 0,
        0, 17000,
        8000, 15700,
        9000, 15000,
        10000, 14000,
        10000, 0,
        0, 0
      ),
      ncol = 2,
      byrow = TRUE
    )
  ))
  
  # Example polygon for right tissue
  # Replace these coordinates with your actual right-side boundary
  poly_right <- st_polygon(list(
    matrix(
      c(
        10000, 0,
        10000, 14000,
        9000, 15000,
        8000, 15700,
        4000, 22000,
        16000,22000,
        16000,0,
        10000, 0
      ),
      ncol = 2,
      byrow = TRUE
    )
  ))
  
  # can add more polygons to the list in the same way as before
  polygons <- list(
    left_tissue = poly_left,
    right_tissue = poly_right
  )
  
  saveRDS(polygons, poly_file)
}

# ------------------------------------------------------------
# Convert polygons to sf objects and assign tissue labels
# Using st_intersects includes points on the boundary too
# ------------------------------------------------------------
for (tissue_name in names(polygons)) {
  poly_sf <- st_sfc(polygons[[tissue_name]])
  
  inside <- st_intersects(pts, poly_sf, sparse = FALSE)[, 1]
  
  object$tissue[inside] <- tissue_name
}

# ------------------------------------------------------------
# Tissue separation visualization
# ------------------------------------------------------------


#code to plot the polygons in the future

#has capabilities for more than just the left and right polygon mapping
#but this is unused currently
plot_tissue_polygons <- function(
    xy,
    pts,
    polygons,
    sample_name,
    colors = c("red", "blue", "green", "purple", "orange")
) {
  
  # default color = unassigned
  bin_color <- rep("black", nrow(xy))
  
  # store counts for legend
  tissue_counts <- c()
  
  # color bins and draw polygon outlines
  for (i in seq_along(polygons)) {
    
    tissue_name <- names(polygons)[i]
    
    inside <- st_intersects(
      pts,
      st_sfc(polygons[[i]]),
      sparse = FALSE
    )[, 1]
    
    bin_color[inside] <- colors[i]
    
    tissue_counts[tissue_name] <- sum(inside)
  }
  
  # plot bins
  plot(
    xy$x,
    xy$y,
    col = bin_color,
    pch = 16,
    cex = 0.2,
    asp = 1,
    xlab = "X coordinate",
    ylab = "Y coordinate",
    main = paste0(sample_name, "\nTissue Separation")
  )
  
  # draw outlines
  for (i in seq_along(polygons)) {
    
    plot(
      st_sfc(polygons[[i]]),
      add = TRUE,
      border = colors[i],
      lwd = 3
    )
  }
  
  # legend
  legend(
    "topright",
    legend = c(
      paste0(
        names(tissue_counts),
        " (n=",
        tissue_counts,
        ")"
      ),
      "Unassigned"
    ),
    col = c(
      colors[seq_along(tissue_counts)],
      "black"
    ),
    pch = 16,
    bty = "n"
  )
}




#code to save plot of tissue seperation
pdf(
  file.path(sample_output_dir, "tissue_sep_visualization.pdf"),
  width = 8,
  height = 8
)

plot_tissue_polygons(
  xy = xy,
  pts = pts,
  polygons = polygons,
  sample_name = sample_name
)

dev.off()
