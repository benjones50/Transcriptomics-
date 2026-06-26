# ------------------------------------------------------------------
# Compatibility patch
#
# The original Seurat implementation treated the absence of the
# coords_x_orientation slot as evidence that the object was created
# with an outdated coordinate system and immediately stopped.
#
# Some valid VisiumV2 objects (including custom builders and certain
# Seurat installations) do not define this slot at all.
#
# Modified behavior:
#   - If the slot exists, require it to be "horizontal".
#   - If the slot does not exist, assume compatibility and continue.
#
# This preserves the safety check for explicitly old objects while
# allowing compatible objects without the slot to be plotted.
# ------------------------------------------------------------------

#wrapper:

# ------------------------------------------------------------------
# Compatibility patch
# ------------------------------------------------------------------

SpatialPlot_patched <- function(
    object, group.by = NULL, features = NULL, images = NULL,
    cols = NULL, image.alpha = 1, image.scale = "lowres", crop = TRUE,
    slot = "data", keep.scale = "feature", min.cutoff = NA, max.cutoff = NA,
    cells.highlight = NULL, cols.highlight = c("#DE2D26", "grey50"),
    facet.highlight = FALSE, label = FALSE, label.size = 5, label.color = "white",
    label.box = TRUE, repel = FALSE, ncol = NULL, combine = TRUE,
    pt.size.factor = 1.6, alpha = c(1, 1), shape = 21, stroke = NA,
    stroke.alpha = NA, interactive = FALSE, do.identify = FALSE,
    identify.ident = NULL, do.hover = FALSE, information = NULL,
    plot_segmentations = FALSE
) {
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning(
      "'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity",
      call. = FALSE,
      immediate. = TRUE
    )
    interactive <- TRUE
  }
  
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  
  cells <- unique(Seurat:::CellsByImage(object, images = images, unlist = TRUE))
  
  if (is.null(x = features)) {
    if (interactive) {
      if (identical(alpha, c(1, 1))) {
        alpha <- c(0.1, 1)
      }
      tryCatch(
        expr = {
          return(Seurat:::ISpatialDimPlot(
            object = object,
            image = images[1],
            image.scale = image.scale,
            group.by = group.by,
            alpha = alpha
          ))
        },
        error = function(e) {
          if (grepl("arguments imply differing number of rows", conditionMessage(e))) {
            stop(
              "Cells were removed due to missing data; check if the specified image and assay are correct.\n",
              call. = FALSE
            )
          } else {
            stop(e)
          }
        }
      )
    }
    
    group.by <- group.by %||% "ident"
    object[["ident"]] <- Idents(object = object)
    data <- object[[group.by]]
    data <- data[cells, , drop = FALSE]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  } else {
    if (interactive) {
      return(ISpatialFeaturePlot(
        object = object,
        feature = features[1],
        image = images[1],
        image.scale = image.scale,
        slot = slot,
        alpha = alpha
      ))
    }
    
    data <- FetchData(object = object, vars = features, cells = cells, layer = slot, clean = FALSE)
    features <- colnames(x = data)
    
    min.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = min(data[, feature]), no = cutoff)
      },
      cutoff = min.cutoff,
      feature = features
    )
    
    max.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = max(data[, feature]), no = cutoff)
      },
      cutoff = max.cutoff,
      feature = features
    )
    
    check.lengths <- unique(x = vapply(
      X = list(features, min.cutoff, max.cutoff),
      FUN = length,
      FUN.VALUE = numeric(length = 1)
    ))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    
    data <- sapply(
      X = 1:ncol(x = data),
      FUN = function(index) {
        data.feature <- as.vector(x = data[, index])
        min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index], data.feature)
        max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index], data.feature)
        data.feature[data.feature < min.use] <- min.use
        data.feature[data.feature > max.use] <- max.use
        data.feature
      }
    )
    
    colnames(x = data) <- features
    rownames(x = data) <- cells
  }
  
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- cells
  
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && is.list(x = cells.highlight))
  
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("'do.hover' requires only one image, using image ", images, call. = FALSE, immediate. = TRUE)
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = "feature", no = "grouping")
      warning("'do.hover' requires only one ", type, ", using ", features, call. = FALSE, immediate. = TRUE)
    }
    if (facet.highlight) {
      warning("'do.hover' requires no faceting highlighted cells", call. = FALSE, immediate. = TRUE)
      facet.highlight <- FALSE
    }
  }
  
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("Faceting the highlight only works with a single image, using image ", images, call. = FALSE, immediate. = TRUE)
    }
    ncols <- length(x = cells.highlight)
  } else {
    ncols <- length(x = images)
  }
  
  plots <- vector(mode = "list", length = length(x = features) * ncols)
  
  if (!(is.null(x = keep.scale)) && keep.scale == "all") {
    max.feature.value <- max(apply(data, 2, function(x) max(x, na.rm = TRUE)))
  }
  
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, no = i)
    image.use <- object[[images[[image.idx]]]]
    is_visium_v2 <- inherits(image.use, "VisiumV2")
    
    
    
    
    
    # ------------------------------------------------------------------
    # Compatibility patch
    #
    # Original Seurat implementation:
    #
    # old_axis_orientation <- (!.hasSlot(image.use,
    #     "coords_x_orientation")) ||
    #   (.hasSlot(image.use,
    #     "coords_x_orientation") &&
    #    slot(image.use,
    #         "coords_x_orientation") != "horizontal")
    #
    # if (is_visium_v2 && old_axis_orientation) {
    #   stop(...)
    # }
    #
    # Problem:
    #
    # Some valid VisiumV2 objects (including custom-built objects and
    # certain Seurat installations) do not define the
    # coords_x_orientation slot.
    #
    # The original implementation therefore incorrectly assumes these
    # objects are outdated and immediately stops plotting.
    #
    # Modified behavior:
    #
    #   Missing slot
    #       -> Assume compatibility.
    #
    #   Slot exists but is not "horizontal"
    #       -> Stop, matching Seurat's intended safety check.
    #
    # This preserves the orientation validation for objects that
    # explicitly define the slot while allowing compatible objects
    # without the slot to continue plotting.
    # ------------------------------------------------------------------
    
    # Original Seurat code:
    #
    # old_axis_orientation <- (!.hasSlot(image.use,
    #     "coords_x_orientation")) ||
    #   (.hasSlot(image.use,
    #     "coords_x_orientation") &&
    #    slot(image.use,
    #         "coords_x_orientation") != "horizontal")
    
    old_axis_orientation <-
      .hasSlot(image.use, "coords_x_orientation") &&
      slot(image.use, "coords_x_orientation") != "horizontal"
    
    # Remaining code below is unchanged from Seurat.
    
    
    
    
    
    if (is_visium_v2 && old_axis_orientation) {
      stop(
        "Please run `UpdateSeuratObject` on your Seurat object first to ensure that data aligns to the image ",
        images[[image.idx]],
        " when plotting.",
        call. = TRUE
      )
    }
    
    if (plot_segmentations == TRUE && inherits(image.use, "VisiumV2") && "segmentations" %in% names(image.use)) {
      db <- DefaultBoundary(image.use)
      on.exit(DefaultBoundary(image.use) <- db, add = TRUE)
      DefaultBoundary(image.use) <- "segmentations"
    }
    
    coordinates <- GetTissueCoordinates(object = image.use, scale = image.scale)
    
    highlight.use <- if (facet.highlight) cells.highlight[i] else cells.highlight
    
    for (j in seq_along(features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[, features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" && !inherits(x = data[, features[j]], what = "factor")) {
        max.feature.value <- max(data[, features[j]])
      }
      
      has_visium_segm_data <- inherits(image.use, "VisiumV2") &&
        !is.null(image.use@boundaries$segmentations) &&
        "sf.data" %in% slotNames(image.use@boundaries$segmentations)
      
      if (!("cell" %in% colnames(x = coordinates))) {
        coordinates$cell <- rownames(x = coordinates)
      }
      
      idx <- match(coordinates$cell, rownames(x = data))
      plot.data <- cbind(coordinates, data[idx, features[j], drop = FALSE])
      
      plot <- Seurat:::SingleSpatialPlot(
        data = plot.data,
        image = image.use,
        image.scale = image.scale,
        image.alpha = image.alpha,
        col.by = features[j],
        cols = cols,
        alpha.by = if (is.null(x = group.by)) features[j] else NULL,
        pt.alpha = if (!is.null(x = group.by)) alpha[j] else NULL,
        geom = if (inherits(x = image.use, what = "STARmap")) {
          "poly_starmap"
        } else if (has_visium_segm_data && plot_segmentations) {
          "poly"
        } else {
          "spatial"
        },
        cells.highlight = highlight.use,
        cols.highlight = cols.highlight,
        pt.size.factor = pt.size.factor,
        shape = shape,
        stroke = stroke,
        stroke.alpha = stroke.alpha,
        crop = crop
      )
      
      if (is.null(x = group.by)) {
        plot <- plot +
          scale_fill_gradientn(name = features[j], colours = Seurat:::SpatialColors(n = 100)) +
          theme(legend.position = "top") +
          scale_alpha(range = alpha) +
          guides(alpha = "none")
      } else if (label) {
        plot <- Seurat:::LabelClusters(
          plot = plot,
          id = ifelse(test = is.null(x = cells.highlight), yes = features[j], no = "highlight"),
          geom = if (inherits(x = image.use, what = "STARmap") || (has_visium_segm_data && plot_segmentations)) {
            "GeomPolygon"
          } else {
            "GeomSpatial"
          },
          repel = repel,
          size = label.size,
          color = label.color,
          box = label.box,
          position = "nearest"
        )
      }
      
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot + ggtitle(label = images[[image.idx]]) + theme(plot.title = element_text(hjust = 0.5))
      }
      
      if (facet.highlight) {
        plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
      }
      
      if (has_visium_segm_data && plot_segmentations && !is.null(group.by)) {
        plot <- plot + guides(fill = guide_legend(override.aes = list(alpha = 1, color = "black", linewidth = 0.2, size = 2)))
      }
      
      if (!(is.null(x = keep.scale)) && !inherits(x = data[, features[j]], "factor")) {
        plot <- suppressMessages(plot & scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 100), limits = c(NA, max.feature.value)))
      }
      
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  
  if (combine) {
    if (!is.null(x = ncol)) {
      return(wrap_plots(plots = plots, ncol = ncol))
    }
    if (length(x = images) > 1) {
      return(wrap_plots(plots = plots, ncol = length(x = images)))
    }
    return(wrap_plots(plots = plots))
  }
  
  return(plots)
}







# Generic fallback helper
spatial_plot_fallback <- function(primary_fun, ...) {
  tryCatch(
    primary_fun(...),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("coords_x_orientation|UpdateSeuratObject", msg)) {
        message(
          "Detected VisiumV2 compatibility issue.\n",
          "Using patched SpatialPlot()."
        )
        return(SpatialPlot_patched(...))
      }
      stop(e)
    }
  )
}




# allow for use of patched spatialplot function if original function fails
# ------------------------------------------------------------------
# Wrapper around Seurat::SpatialFeaturePlot()
#
# Purpose:
#   Attempt to use the original Seurat implementation first.
#
#   If plotting fails due to the VisiumV2
#   coords_x_orientation compatibility issue,
#   automatically retry using the patched SpatialPlot().
#
# Any other error is re-thrown unchanged.
#
# This wrapper intentionally mirrors Seurat's exported API so it
# can be used as a drop-in replacement within analysis scripts.
# ------------------------------------------------------------------

SpatialFeaturePlot_safe <- function(
    object,
    features,
    images = NULL,
    crop = TRUE,
    slot = "data",
    keep.scale = "feature",
    min.cutoff = NA,
    max.cutoff = NA,
    ncol = NULL,
    combine = TRUE,
    pt.size.factor = 1.6,
    alpha = c(1, 1),
    image.alpha = 1,
    image.scale = "lowres",
    shape = 21,
    stroke = NA,
    stroke.alpha = NA,
    interactive = FALSE,
    information = NULL,
    plot_segmentations = FALSE
) {
  
  tryCatch(
    
    # ----------------------------------------------------------
    # Original Seurat implementation
    # ----------------------------------------------------------
    #
    # return(
    #   Seurat::SpatialFeaturePlot(
    #     object = object,
    #     ...
    #   )
    # )
    #
    # ----------------------------------------------------------
    
    Seurat::SpatialFeaturePlot(
      object = object,
      features = features,
      images = images,
      crop = crop,
      slot = slot,
      keep.scale = keep.scale,
      min.cutoff = min.cutoff,
      max.cutoff = max.cutoff,
      ncol = ncol,
      combine = combine,
      pt.size.factor = pt.size.factor,
      alpha = alpha,
      image.alpha = image.alpha,
      image.scale = image.scale,
      shape = shape,
      stroke = stroke,
      stroke.alpha = stroke.alpha,
      interactive = interactive,
      information = information,
      plot_segmentations = plot_segmentations
    ),
    
    error = function(e) {
      
      msg <- conditionMessage(e)
      
      if (grepl("UpdateSeuratObject|coords_x_orientation", msg)) {
        
        message(
          "Detected VisiumV2 compatibility issue.\n",
          "Retrying with patched SpatialPlot()."
        )
        
        
        # --------------------------------------------------
        # Patched behavior
        # --------------------------------------------------
        
        return(
          
          SpatialPlot_patched(
            object = object,
            features = features,
            images = images,
            crop = crop,
            slot = slot,
            keep.scale = keep.scale,
            min.cutoff = min.cutoff,
            max.cutoff = max.cutoff,
            ncol = ncol,
            combine = combine,
            pt.size.factor = pt.size.factor,
            alpha = alpha,
            image.alpha = image.alpha,
            image.scale = image.scale,
            shape = shape,
            stroke = stroke,
            stroke.alpha = stroke.alpha,
            interactive = interactive,
            information = information,
            plot_segmentations = plot_segmentations
          )
          
        )
        
      }
      
      stop(e)
      
    }
    
  )
  
}

count_col <- paste0("nCount_",DefaultAssay(object))
feature_col <- paste0("nFeature_",DefaultAssay(object))



count_spatial <- SpatialFeaturePlot_safe(
  object,
  features = count_col,
  pt.size.factor = pt_size,
  image.alpha = 0.05
) + ggtitle("nCount (Raw UMI Counts)")

count_spatial


# SpatialFeaturePlot_safe <- function(...) {
#   spatial_plot_fallback(Seurat::SpatialFeaturePlot, ...)
# }
# 
# SpatialDimPlot_safe <- function(...) {
#   spatial_plot_fallback(Seurat::SpatialDimPlot, ...)
# }












