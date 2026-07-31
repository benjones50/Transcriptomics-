
# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------


prep_de_df <- function(df) {
  
  # Convert to a standard data.frame
  df <- as.data.frame(df)
  
  #------------------------------------------------------------
  # Ensure a gene column exists
  #------------------------------------------------------------
  # Seurat stores gene names as row names rather than a dedicated
  # column. If no gene column exists, copy the row names into one.
  if (!"gene" %in% colnames(df)) {
    df$gene <- rownames(df)
  }
  
  
  #------------------------------------------------------------
  # Replace problematic values
  #------------------------------------------------------------
  
  # Replace NA or zero adjusted p-values with the smallest positive
  # floating-point number. This avoids issues when plotting
  # -log10(p-value), since log10(0) is undefined.
  df$p_val_adj[df$p_val_adj < 1e-300] <- 1e-300
  
  # Replace missing fold changes with zero so plotting functions
  # always receive valid numeric values.
  df$avg_log2FC[
    is.na(df$avg_log2FC)
  ] <- 0
  
  # Return the cleaned and standardized DE table.
  df
}




# ----------------------------------------------------------------------
# Volcano plotter
#
# Take one DE table
# ↓
# Clean and standardize the data
# ↓
# Identify significantly differentially expressed genes
# ↓
# Select a small set of genes to label
# ↓
# Determine appropriate x-axis limits
# ↓
# build the figure
# ↓
# Return the ggplot object
# ----------------------------------------------------------------------
make_volcano_plot <- function(df,
                              title,
                              p_cutoff = 0.05,
                              fc_cutoff = 0.25,
                              n_labels = 5) {
  
  # Standardize column names, convert data types, and replace
  # problematic values (e.g. NA or zero adjusted p-values).
  df <- prep_de_df(df)
  
  #------------------------------------------------------------
  # Identify significant genes
  #------------------------------------------------------------
  
  # Keep only genes that pass both the adjusted p-value cutoff
  # and the absolute log2 fold-change cutoff. These are the
  # genes considered biologically and statistically significant.
  # Keep significant genes
  sig <- df %>%
    filter(
      p_val_adj < p_cutoff,
      abs(avg_log2FC) >= fc_cutoff
    ) %>%
    arrange(desc(abs(avg_log2FC)))
  
  # Label the strongest effect sizes
  select_labels <- sig %>%
    slice_head(n = n_labels) %>%
    pull(gene)
  
  #------------------------------------------------------------
  # Choose genes to label
  #------------------------------------------------------------
  
  # Significant upregulated genes
  sig_up <- sig %>%
    filter(avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  # Significant downregulated genes
  sig_down <- sig %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC)
  
  # Label the strongest genes in each direction
  select_labels <- c(
    sig_up %>%
      slice_head(n = n_labels) %>%
      pull(gene),
    
    sig_down %>%
      slice_head(n = n_labels) %>%
      pull(gene)
  )
  
  select_labels <- unique(select_labels)
  
  # If no significant genes were found, label the genes with the
  # largest absolute fold changes so the plot still contains
  # informative annotations.
  if (length(select_labels) == 0) {
    
    select_labels <- df %>%
      arrange(desc(abs(avg_log2FC))) %>%
      slice_head(n = 2 * n_labels) %>%
      pull(gene) %>%
      unique()
    
  }
  
  
  
  #------------------------------------------------------------
  # Determine x-axis limits
  #------------------------------------------------------------
  
  # Find the largest absolute log2 fold change so the plot is
  # centered around zero with symmetric limits.
  xlim_max <- max(abs(df$avg_log2FC), na.rm = TRUE)
  
  # If every fold change is zero (or invalid), use a default
  # range so the plot can still be drawn.
  if (!is.finite(xlim_max) || xlim_max == 0) {
    xlim_max <- 1
  }
  
  
  #------------------------------------------------------------
  # Create volcano plot
  #------------------------------------------------------------
  
  # Flag significant genes
  df$significant <-
    df$p_val_adj < p_cutoff &
    abs(df$avg_log2FC) >= fc_cutoff
  
  # Symmetric limits for x-axis and color scale
  max_fc <- max(abs(df$avg_log2FC), na.rm = TRUE)
  
  p <- ggplot(
    df,
    aes(
      x = avg_log2FC,
      y = -log10(p_val_adj)
    )
  ) +
    
    # Points
    geom_point(
      aes(
        color = avg_log2FC,
        alpha = significant
      ),
      size = 1.6
    ) +
    
    # Continuous log2FC color scale
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-3, 3),
      oob = scales::squish,
      name = expression(Log[2]~FC)
    ) +
    
    # Fade non-significant genes
    scale_alpha_manual(
      values = c(
        "TRUE" = 0.90,
        "FALSE" = 0.20
      ),
      guide = "none"
    ) +
    
    # Threshold lines
    geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    # Gene labels
    geom_text_repel(
      data = subset(df, gene %in% select_labels),
      aes(label = gene),
      color = "black",
      size = 3.2,
      fontface = "bold",
      box.padding = 0.4,
      point.padding = 0.2,
      max.overlaps = Inf
    ) +
    
    # Axes
    scale_x_continuous(
      limits = c(-(max_fc + 0.5), max_fc + 0.5),
      breaks = scales::pretty_breaks(n = 11)
    ) +
    
    # Labels
    labs(
      title = title,
      x = expression(Log[2]~fold~change),
      y = expression(-Log[10]~adjusted~italic(P)),
      caption = paste0("n genes = ", nrow(df))
    ) +
    
    # Theme
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  p
}







# ------------------------------------------------------------
# Cluster vs rest, FindAllMarkers
#
# Idents(object) <- seurat_cluster.sketched
# ↓
# FindAllMarkers()
# ↓
# one Wilcoxon test for every gene / cluster
# ↓
# one large dataframe with gene, cluster, avg_log2FC, p_val_adj
# ------------------------------------------------------------
run_cluster_vs_rest_de <- function(object,
                                   ident_col,
                                   assay = NULL,
                                   min_pct = 0.01,
                                   test_use = "wilcox") {
  if (is.null(assay)) assay <- DefaultAssay(object)
  
  Idents(object) <- ident_col
  
  markers <- FindAllMarkers(
    object,
    assay = assay,
    only.pos = FALSE,
    logfc.threshold = 0.1,
    min.pct = min_pct,
    test.use = test_use,
    #return.thresh = 0.01,
    verbose = TRUE
  )
  
  markers <- prep_de_df(markers)
  
  markers
}







# ------------------------------------------------------------
# for each cluster:
# extract cluster info
# ↓
# make_volcano_plot()
# ↓
# store plot
# wrap_plots()
# ↓
# save_plot()
# ------------------------------------------------------------
plot_cluster_vs_rest_volcanoes <- function(cluster_vs_rest_df,
                                           plot_dir,
                                           file_stub = "cluster_vs_rest_marker_grid",
                                           width = 27,
                                           height = 18,
                                           p_cutoff = 0.05,
                                           fc_cutoff = 0.25,
                                           n_labels = 5) {
  
  cluster_vs_rest_df <- prep_de_df(cluster_vs_rest_df)
  
  if (!"cluster" %in% colnames(cluster_vs_rest_df)) {
    stop("cluster_vs_rest_df must contain a 'cluster' column from FindAllMarkers().")
  }
  
  clusters <- unique(as.character(cluster_vs_rest_df$cluster))
  
  plot_list <- lapply(clusters, function(cl) {
    df_cl <- cluster_vs_rest_df %>% filter(as.character(cluster) == cl)
    
    make_volcano_plot(
      df = df_cl,
      title = paste0("Cluster ", cl, " vs rest"),
      p_cutoff = p_cutoff,
      fc_cutoff = fc_cutoff,
      n_labels = n_labels
    )
  })
  
  grid_ncol <- ceiling(sqrt(length(plot_list)))
  cluster_vs_rest_grid <- wrap_plots(plotlist = plot_list, ncol = grid_ncol)
  
  save_plot(
    plot_object = cluster_vs_rest_grid,
    file_stub = file_stub,
    save_dir = file.path(plot_dir, "differential_expression", "cluster_vs_rest"),
    width = width,
    height = height
  )
  
  invisible(cluster_vs_rest_grid)
}

# ------------------------------------------------------------
#    
#  Pairwise DE, FindMarkers again

# 
# FindMarkers(0 vs 1)
# FindMarkers(0 vs 2)
# FindMarkers(0 vs 3)
# ...
# FindMarkers(7 vs 8)
# 
# Output: 
# a list pairwise_results
# ├── 0_vs_1
# ├── 0_vs_2
# ├── 0_vs_3
# ├── ...
# └── 7_vs_8
# 
# Each element is its own DE table.
# ------------------------------------------------------------

run_pairwise_de <- function(object,
                            ident_col,
                            assay = NULL,
                            min_pct = 0.01,
                            test_use = "wilcox") {
  if (is.null(assay)) assay <- DefaultAssay(object)
  
  Idents(object) <- ident_col
  cluster_levels <- levels(Idents(object))
  cluster_pairs <- combn(cluster_levels, 2, simplify = FALSE)
  
  pairwise_results <- list()
  
  for (pair in cluster_pairs) {
    cl1 <- pair[[1]]
    cl2 <- pair[[2]]
    comp_name <- paste0(cl1, "_vs_", cl2)
    
    message("Running ", comp_name)
    
    de_df <- FindMarkers(
      object,
      ident.1 = cl1,
      ident.2 = cl2,
      assay = assay,
      only.pos = FALSE,
      min.pct = min_pct,
      test.use = test_use,
      #return.thresh = 1,
      verbose = FALSE
    )
    
    pairwise_results[[comp_name]] <- prep_de_df(de_df)
  }
  
  pairwise_results
}






# ------------------------------------------------------------
# for every pair:
#
# pairwise_results[[pair]]
# ↓
# make_volcano_plot()
# ↓
# place in matrix
# ↓
# wrap_plots()
# ↓
# save_plot()
# ------------------------------------------------------------
plot_pairwise_volcano_grid <- function(pairwise_results,
                                       object,
                                       ident_col,
                                       plot_dir = plot_dir,
                                       file_stub = "pairwise_volcano_grid",
                                       width = 50,
                                       height = 50,
                                       p_cutoff = 0.05,
                                       fc_cutoff = 0.25,
                                       n_labels = 5) {
  Idents(object) <- ident_col
  cluster_levels <- levels(Idents(object))
  n_clust <- length(cluster_levels)
  
  plot_mat <- vector("list", n_clust * n_clust)
  k <- 1
  
  for (i in seq_len(n_clust)) {
    for (j in seq_len(n_clust)) {
      if (j <= i) {
        plot_mat[[k]] <- plot_spacer()
      } else {
        key <- paste0(cluster_levels[i], "_vs_", cluster_levels[j])
        
        if (!key %in% names(pairwise_results)) {
          alt_key <- paste0(cluster_levels[j], "_vs_", cluster_levels[i])
          if (!alt_key %in% names(pairwise_results)) {
            warning("Missing pairwise result for ", key, " and ", alt_key)
            plot_mat[[k]] <- plot_spacer()
            k <- k + 1
            next
          }
          
          df <- pairwise_results[[alt_key]]
          title <- paste0("Cluster ", cluster_levels[j], " vs Cluster ", cluster_levels[i])
        } else {
          df <- pairwise_results[[key]]
          title <- paste0("Cluster ", cluster_levels[i], " vs Cluster ", cluster_levels[j])
        }
        
        plot_mat[[k]] <- make_volcano_plot(
          df = df,
          title = title,
          p_cutoff = p_cutoff,
          fc_cutoff = fc_cutoff,
          n_labels = n_labels
        )
      }
      k <- k + 1
    }
  }
  
  pairwise_grid <- wrap_plots(plotlist = plot_mat, ncol = n_clust)
  
  save_plot(
    plot_object = pairwise_grid,
    file_stub = file_stub,
    save_dir = file.path(plot_dir, "differential_expression", "pairwise"),
    width = width,
    height = height,
    dpi = 300
  )
  
  invisible(pairwise_grid)
}





# ------------------------------------------------------------
# Example run
# ------------------------------------------------------------

#----------------------------------------------------------
# Run differential expression
#----------------------------------------------------------


cluster_vs_rest_markers <-
  run_cluster_vs_rest_de(
    object = mega_obj,
    ident_col = "seurat_cluster.sketched"
  )

# pairwise_results <-
#   run_pairwise_de(
#     object = mega_obj,
#     ident_col = "seurat_cluster.sketched"
#   )



# Plot dir
plot_dir <- get_figures_dir(projection_dir)

#----------------------------------------------------------
# Create figures
#----------------------------------------------------------

cluster_vs_rest_grid <-
  plot_cluster_vs_rest_volcanoes(
    cluster_vs_rest_df = cluster_vs_rest_markers,
    plot_dir = plot_dir
  )

# pairwise_grid <-
#   plot_pairwise_volcano_grid(
#     pairwise_results = pairwise_results,
#     object = mega_obj,
#     ident_col = "seurat_cluster.sketched",
#     plot_dir = plot_dir
#   )



# #code to test a single volcano, used for testing 
# cluster0_df <- cluster_vs_rest_markers %>%
#   filter(cluster == "0")
# 
# cluster0_plot <- make_volcano_plot(
#   df = cluster0_df,
#   title = "Cluster 0 vs Rest"
# )
# 
# save_plot(
#   plot_object = cluster0_plot,
#   file_stub = "test_cluster0_plot",
#   save_dir = file.path(plot_dir, "differential_expression"),
#   width = 12,
#   height = 6
# )







