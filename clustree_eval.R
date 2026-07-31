source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")
library(clustree)

#this code was used to test to make sure my pipeline was creating the same clusters using the same seeds all the way through


options(future.globals.maxSize = 1000 * 1024^2)

# Decide on embedding
max_dims <- 11

cluster_resolutions <- c(0.2,0.4,0.5,0.6,0.8,1.0,1.2)

selected_cluster_resolution <- 0.5



run_pipeline <- function() {
  
  obj <- load_mega_object(
    directory = mega_dir
  )
  
  obj <- load_mega_object(directory = mega_dir)
  
  obj <- run_sketch_pca(
    object = obj,
    ncells = 50000
  )


  obj <- run_harmony_embedding(
    object = obj,
    ndims = max_dims,
    sketch_reduction = "pca.sketch",
    group.by.vars = "tissue"
  )
  
  
  obj <- run_sketch_clustering(
    object = obj,
    sketch_reduction = embedding_reduction,
    dims = 1:max_dims,
    resolutions = cluster_resolutions,
    selected_cluster_resolution = selected_cluster_resolution
  )
  
  obj
  
  #TODO move into above function so object saved and loaded is way smaller
  #takes only cells used in the sketch
  
  cells <- WhichCells(
    obj,
    expression = !is.na(seurat_cluster.sketched)
  )
  
  #takes only cells used in the sketch
  obj <- subset(obj, cells = cells)
  
}



#prints comparisons of 2 objects 
compare_objects <- function(obj1, obj2) {
  
  cat("====================================\n")
  cat("Comparing Seurat objects\n")
  cat("====================================\n\n")
  
  #----------------------------------------------------------
  # Sketch comparison
  #
  # Compare the sketched cells, make sure the same
  #----------------------------------------------------------
  cat(
    "Sketch identical: ",
    identical(
      colnames(obj1[["sketch"]]),
      colnames(obj2[["sketch"]])
    ),
    "\n\n"
  )
  
  #----------------------------------------------------------
  # PCA comparison
  #
  # Compare the PCA coordinates for every sketch cell.
  #
  #----------------------------------------------------------
  cat(
    "PCA embeddings identical: ",
    isTRUE(
      all.equal(
        Embeddings(obj1, "pca.sketch"),
        Embeddings(obj2, "pca.sketch")
      )
    ),
    "\n\n"
  )
  
  #----------------------------------------------------------
  # Harmony comparison
  #
  # Harmony is iterative and may introduce small numerical
  # differences?
  #----------------------------------------------------------
  if (
    "harmony.sketch" %in% Reductions(obj1) &&
    "harmony.sketch" %in% Reductions(obj2)
  ) {
    
    cat(
      "Harmony embeddings identical: ",
      isTRUE(
        all.equal(
          Embeddings(obj1, "harmony.sketch"),
          Embeddings(obj2, "harmony.sketch")
        )
      ),
      "\n\n"
    )
  }
  
  #----------------------------------------------------------
  # Cluster comparison
  #
  # Compare every clustering resolution present in both
  # objects.
  #
  # For each resolution we report:
  #
  #  • whether every cell has the exact same cluster
  #    assignment
  #
  #  • how many cells changed cluster
  #
  #  • what percentage of cells changed
  #
  #  • the total number of clusters produced in each run
  #
  # This distinguishes:
  #
  #   - identical clustering
  #   - same number of clusters but different assignments
  #   - completely different clustering structures
  #----------------------------------------------------------
  
  cluster_cols <- intersect(
    grep("^cluster_res_", colnames(obj1[[]]), value = TRUE),
    grep("^cluster_res_", colnames(obj2[[]]), value = TRUE)
  )
  
  cat("Cluster comparison:\n\n")
  
  for (col in cluster_cols) {
    
    # Cluster assignments from each object
    x <- obj1[[col]][, 1]
    y <- obj2[[col]][, 1]
    
    # Number of unique clusters produced
    n_clusters_1 <- length(unique(na.omit(x)))
    n_clusters_2 <- length(unique(na.omit(y)))
    
    # TRUE only if every cell has exactly the same label
    identical_assignments <- identical(x, y)
    
    # Compare only cells present in both objects
    valid_cells <- !is.na(x) & !is.na(y)
    
    # Cells assigned to different clusters
    changed <- valid_cells & (x != y)
    
    # Number of changed cells
    n_changed <- sum(changed)
    
    # Percentage of cells with different assignments
    pct_changed <- 100 * n_changed / sum(valid_cells)
    
    cat(
      sprintf(
        "%-18s identical=%-5s changed=%6d (%.3f%%) clusters=%2d vs %2d\n",
        col,
        identical_assignments,
        n_changed,
        pct_changed,
        n_clusters_1,
        n_clusters_2
      )
    )
  }
  
  cat("\n")
  
  #----------------------------------------------------------
  # Metadata comparison
  #
  # Display all metadata columns shared by both objects.
  #
  # This is useful for verifying that expected clustering
  # columns and annotations exist before comparing objects.
  #----------------------------------------------------------
  
  cat("Shared metadata columns:\n")
  
  print(
    intersect(
      colnames(obj1[[]]),
      colnames(obj2[[]])
    )
  )
  
  cat("\nComparison complete.\n")
}

obj1 <- run_pipeline()
obj2 <- run_pipeline()

compare_objects(obj1, obj2)
