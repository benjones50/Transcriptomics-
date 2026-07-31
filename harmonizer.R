library(harmony)
run_harmony_embedding <- function(
    object,
    ndims,
    sketch_reduction = "pca.sketch",
    harmony_reduction = "harmony",
    group.by.vars = "tissue",
    theta = NULL,
    lambda = NULL,
    sigma = 0.1,
    nclust = NULL,
    max_iter = 10L,
    early_stop = TRUE,
    reference_values = NULL,
    project.dim = TRUE,
    verbose = TRUE,
    ...
)
{
  
  # ==========================================================
  # Validate inputs
  # ==========================================================
  
  if (!sketch_reduction %in% Reductions(object = object)) {
    stop(
      "Sketch reduction not found in object: ",
      sketch_reduction,
      call. = FALSE
    )
  }
  
  if (harmony_reduction %in% Reductions(object = object)) {
    stop(
      "harmony_reduction already exists in object: ",
      harmony_reduction,
      call. = FALSE
    )
  }
  
  if (!group.by.vars %in% colnames(object[[]])) {
    stop(
      "Metadata column not found in object: ",
      group.by.vars,
      call. = FALSE
    )
  }
  
  ndims <- as.integer(x = ndims)
  
  sketch_reduction_object <- object[[sketch_reduction]]
  available_dims <- ncol(x = Embeddings(object = sketch_reduction_object))
  
  if (ndims > available_dims) {
    stop(
      "Requested dimensions exceed those available in reduction '",
      sketch_reduction,
      "'. Requested dims: ",
      ndims,
      "; available dims: ",
      available_dims,
      call. = FALSE
    )
  }
  
  if (ndims < 2L) {
    stop(
      "Harmony requires at least 2 dimensions.",
      call. = FALSE
    )
  }
  

  
  # ==========================================================
  # Run Harmony directly on the sketch PCA
  #
  # RunHarmony() takes the PCA reduction, the metadata variable
  # to remove, and the dimensions to use. This avoids needing a
  # temporary reduction or IntegrateLayers().
  # ==========================================================
  
  object <- tryCatch(
    expr = {
      harmony::RunHarmony(
        object = object,
        group.by.vars = group.by.vars,
        reduction.use = sketch_reduction,
        dims.use = seq_len(ndims),
        reduction.save = harmony_reduction,
        theta = theta,
        lambda = lambda,
        sigma = sigma,
        nclust = nclust,
        max_iter = max_iter,
        early_stop = early_stop,
        reference_values = reference_values,
        project.dim = project.dim,
        verbose = verbose,
        ...
      )
    },
    error = function(e) {
      stop(
        "Unable to run Harmony embedding.\n",
        "Sketch reduction: ", sketch_reduction, "\n",
        "Harmony reduction: ", harmony_reduction, "\n",
        "Group by variable: ", group.by.vars, "\n",
        "Reason: ", e$message,
        call. = FALSE
      )
    }
  )
  message("Harmony embedding complete.")
  
  return(object)
}