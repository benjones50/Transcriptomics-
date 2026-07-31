# ------------------------------------------------------------
# Cluster vs rest, FindAllMarkers
#
# Idents(object) <- seurat_cluster.sketched
# ↓
# presto::wilcoxauc()
# ↓
# one Wilcoxon test for every gene / cluster
# ↓
# one large dataframe with gene, cluster, avg_log2FC, p_val_adj
# ------------------------------------------------------------
run_cluster_vs_rest_de <- function(object,
                                   ident_col,
                                   assay = NULL) {
  
  if (is.null(assay)) assay <- DefaultAssay(object)
  
  Idents(object) <- ident_col
  DefaultAssay(object) <- assay
  
  expr <- GetAssayData(
    object,
    assay = assay,
    layer = "data"
  )
  
  markers <- presto::wilcoxauc(
    X = expr,
    y = Idents(object)
  ) %>%
    dplyr::rename(
      gene = feature,
      cluster = group,
      avg_log2FC = logFC,
      pct.1 = pct_in,
      pct.2 = pct_out,
      p_val = pval,
      p_val_adj = padj,
      auc = auc
    ) %>%
    mutate(
      p_val_adj = pmax(p_val_adj, 1e-300),
      avg_log2FC = dplyr::coalesce(avg_log2FC, 0)
    ) %>%
    as.data.frame()
  
  markers
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
  
  
  # auc strength
  df$auc_strength <- abs(df$auc - 0.5) * 2
  
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
        color = auc_strength,
        alpha = abs(avg_log2FC)
      ),
      size = 1.6
    ) +
    
    # Continuous log2FC color scale
    scale_color_gradient(
      low = "grey90",
      high = "darkred",
      limits = c(0, 1),
      oob = scales::squish,
      name = "AUC strength"
    ) +
    scale_alpha_continuous(
      range = c(0.2, 1),
      name = expression("|Log"[2]*"FC|")
    )+
    
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
    save_dir = file.path(plot_dir,"differential_expression"),
    width = width,
    height = height
  )
  
  invisible(cluster_vs_rest_grid)
}













