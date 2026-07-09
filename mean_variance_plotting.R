
source(file.path(code_dir, "config.R"))


#loads objects
force_rebuild <- TRUE




object_sct <- load_normalized_object(
  sample_tissue = sample_tissue,
  raw_data_dir = raw_data_dir,
  analysis_mode = analysis_mode,
  bin_size = bin_size,
  normalization = "sct",
  min_counts = min_counts,
  max_counts = max_counts,
  min_features = min_features,
  max_features = max_features,
  max_percent_mt = max_percent_mt,
  force_rebuild = force_rebuild
)

object_log_norm <- load_normalized_object(
  sample_tissue = sample_tissue,
  raw_data_dir = raw_data_dir,
  analysis_mode = analysis_mode,
  bin_size = bin_size,
  normalization = "lognorm",
  min_counts = min_counts,
  max_counts = max_counts,
  min_features = min_features,
  max_features = max_features,
  max_percent_mt = max_percent_mt,
  force_rebuild = force_rebuild
)
  
  
  #mean variance plot
  plot_variable_features <- function(object, nfeatures = 2000, top_n = 10) {
    
    # Find variable features
    object <- FindVariableFeatures(
      object,
      selection.method = "vst",
      nfeatures = nfeatures
    )
    
    # Top variable genes
    top_genes <- head(VariableFeatures(object), top_n)
    
    # Create plot
    plot <- VariableFeaturePlot(object)
    
    # Label top genes
    plot <- LabelPoints(
      plot = plot,
      points = top_genes,
      repel = TRUE
    )
    
    # Return everything needed downstream
    list(
      object = object,
      plot = plot
    )
  }
  
  
  
  # mean variance plot for LogNormalize object
  lognorm_results <- plot_variable_features(object_log_norm)
  
  object_log_norm <- lognorm_results$object
  variable_features_plot_lognorm <- lognorm_results$plot
  
  # mean variance plot for SCT object
  sct_results <- plot_variable_features(object_sct)
  
  object_sct <- sct_results$object
  variable_features_plot_sct <- sct_results$plot
  
  # Combine plots
  variable_features_plot_all <-
    variable_features_plot_lognorm | variable_features_plot_sct
  
  
  #mean variance plot!!!!! very cool
  save_plot(
    variable_features_plot_all,
    paste0("variable_features_plot_all",sample_name),
    sample_output_dir,
    width = 16,
    height = 8,
    dpi = 400
  )
  
