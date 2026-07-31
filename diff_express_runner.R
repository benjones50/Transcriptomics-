source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")



#----------------------------------------------------------
# Run differential expression
#----------------------------------------------------------
cluster_vs_rest_markers <- #from volcano plot code
  run_cluster_vs_rest_de(
    object = mega_obj,
    ident_col = "seurat_cluster.sketched"
  )

# Plot dir
plot_dir <- get_figures_dir(projection_dir)

write.csv(
  cluster_vs_rest_markers,
  file.path(
    plot_dir,
    "differential_expression",
    "cluster_vs_rest_markers.csv"
  ),
  row.names = FALSE
)

#----------------------------------------------------------
# Create figures
#----------------------------------------------------------


cluster_vs_rest_grid <- #creates and saves a volcano plot of clusters vs all other cells
  plot_cluster_vs_rest_volcanoes(
    cluster_vs_rest_df = cluster_vs_rest_markers,
    plot_dir = plot_dir
  )

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

#----------------------------------------
# Get top marker genes
#----------------------------------------
# Keep only strong markers
# Remove non-significant markers (p_val_adj cutoff).
#
# Sort by:
#   highest log fold change first
#     if two genes have the same logFC, the one with higher auc comes first.
#----------------------------------------

get_top_markers <- function(
    processed_markers,
    auc_cutoff = 0.60,
    p_cutoff = 0.05,
    fc_cutoff = 0.25,
    top_n = NULL,
    rank_by = c("logFC", "auc", "auc_strength"),
    clusters = NULL
) {
  
  rank_by <- match.arg(rank_by)

  
  top_markers <- processed_markers %>%
    dplyr::filter(
      (auc >= auc_cutoff | auc <= (1 - auc_cutoff)), #keeps auc at >0.6 or <0.4
      p_val_adj < p_cutoff,
      abs(avg_log2FC) >= fc_cutoff
    )
  
  
  
  # Optionally keep only selected clusters
  if (!is.null(clusters)) {
    top_markers <- top_markers %>%
      dplyr::filter(cluster %in% clusters)
  }
  
  # Rank within each cluster
  top_markers <- top_markers %>%
    dplyr::group_by(cluster)
  
  if (rank_by == "logFC") {
    
    top_markers <- top_markers %>%
      dplyr::arrange(
        dplyr::desc(avg_log2FC),
        dplyr::desc(auc),
        .by_group = TRUE
      )
    
  } else if (rank_by == "auc") {
    
    top_markers <- top_markers %>%
      dplyr::arrange(
        dplyr::desc(auc),
        dplyr::desc(avg_log2FC),
        .by_group = TRUE
      )
    
  } else if (rank_by == "auc_strength") {
    
    top_markers <- top_markers %>%
      dplyr::arrange(
        dplyr::desc(auc_strength),
        dplyr::desc(avg_log2FC),
        .by_group = TRUE
      )
    
  }
  
  # Keep only the top N markers per cluster
  if (!is.null(top_n)) {
    top_markers <- top_markers %>%
      dplyr::slice_head(n = top_n)
  }
  
  top_markers <- top_markers %>%
    dplyr::ungroup()
  
  return(top_markers)
}




#----------------------------------------
# Add marker metrics and define valid markers
#----------------------------------------

auc_cutoff <- 0.60
p_cutoff <- 0.05
fc_cutoff <- 0.25

cluster_vs_rest_markers <- cluster_vs_rest_markers %>%
  dplyr::mutate(
    auc_strength = abs(auc - 0.5) * 2,
    valid_marker =
      (auc >= auc_cutoff | auc <= (1 - auc_cutoff)) &
      p_val_adj < p_cutoff &
      abs(avg_log2FC) >= fc_cutoff
  )

# Keep only valid markers for all summary plots
all_valid_markers <-
  cluster_vs_rest_markers %>%
  dplyr::filter(valid_marker)

#----------------------------------------
# Get valuable marker genes
#----------------------------------------
# top 10 markers for each cluster
top_markers <-
  get_top_markers(
    processed_markers = all_valid_markers,
    top_n = 10,
    rank_by = "auc"
  )

length(top_markers)
top_markers

length(all_valid_markers)
all_valid_markers

# pulls colors from seurat umap for coloring
cluster_colors <- extract_dimplot_colors(
  object = mega_obj,
  reduction = "umap.sketch",
  group.by = "seurat_cluster.sketched"
)

cluster_levels <- names(cluster_colors)

cluster_vs_rest_markers$cluster <- factor(
  cluster_vs_rest_markers$cluster,
  levels = cluster_levels
)

all_valid_markers$cluster <- factor(
  all_valid_markers$cluster,
  levels = cluster_levels
)

top_markers$cluster <- factor(
  top_markers$cluster,
  levels = cluster_levels
)

#----------------------------------------
# Distribution of cluster marker strength
#----------------------------------------

vln_cluster_marker_strength <- ggplot(
  cluster_vs_rest_markers,
  aes(
    x = factor(cluster),
    y = auc,
    fill = factor(cluster)
  )
) +
  geom_violin(
    color = "black",
    trim = FALSE,
    alpha = 0.8
  ) +
  geom_jitter(
    data = dplyr::filter(cluster_vs_rest_markers, !valid_marker),
    color = "grey70",
    width = 0.15,
    height = 0,
    size = 0.8,
    alpha = 0.25
  ) +
  geom_jitter(
    data = all_valid_markers,
    aes(color = factor(cluster)),
    width = 0.15,
    height = 0,
    size = 1.3,
    alpha = 1
  ) +
  geom_hline(
    yintercept = auc_cutoff,
    linetype = "dashed",
    linewidth = 0.6,
    color = "red"
  ) +
  geom_hline(
    yintercept = 1 - auc_cutoff,
    linetype = "dashed",
    linewidth = 0.6,
    color = "red"
  ) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Cluster",
    y = "AUC",
    title = "Distribution of marker strength by cluster"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# save_plot(
#   plot_object = vln_cluster_marker_strength,
#   file_stub = "vln_cluster_marker_strength",
#   save_dir = file.path(plot_dir, "differential_expression"),
#   width = 10,
#   height = 10
# )

#----------------------------------------
# Number of marker genes per cluster
#----------------------------------------

marker_counts <-
  all_valid_markers %>%
  dplyr::count(cluster, .drop = FALSE)

bar_cluster_marker_counts <- ggplot(
  marker_counts,
  aes(
    x = cluster,
    y = n,
    fill = cluster
  )
) +
  geom_col(color = "black") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = cluster_colors,
    drop = FALSE
  ) +
  labs(
    x = "Cluster",
    y = "Marker genes",
    title = "Number of marker genes per cluster"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# save_plot(
#   plot_object = bar_cluster_marker_counts,
#   file_stub = "bar_cluster_marker_counts",
#   save_dir = file.path(plot_dir, "differential_expression"),
#   width = 10,
#   height = 10
# )

scatter_cluster_marker_quality <- ggplot(
  cluster_vs_rest_markers,
  aes(
    x = avg_log2FC,
    y = auc
  )
) +
  geom_point(
    data = dplyr::filter(cluster_vs_rest_markers, !valid_marker),
    color = "black",
    alpha = 0.5,
    size = 2
  ) +
  geom_point(
    data = dplyr::filter(cluster_vs_rest_markers, valid_marker),
    aes(color = factor(cluster)),
    alpha = 1,
    size = 2
  ) +
  geom_hline(
    yintercept = c(1 - auc_cutoff, auc_cutoff),
    linetype = "dashed",
    linewidth = 0.6,
    color = "red"
  ) +
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    linetype = "dashed",
    linewidth = 0.6,
    color = "red"
  ) +
  facet_wrap(~cluster) +
  scale_color_manual(values = cluster_colors) +
  labs(
    x = "Average log2 fold change",
    y = "AUC",
    title = "Marker quality by cluster"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# save_plot(
#   plot_object = scatter_cluster_marker_quality,
#   file_stub = "scatter_cluster_marker_quality",
#   save_dir = file.path(plot_dir, "differential_expression"),
#   width = 12,
#   height = 10
# )

#----------------------------------------
# High-confidence marker logFC by cluster
#----------------------------------------

high_auc_markers <-
  all_valid_markers

high_auc_markers_plot <- ggplot(
  high_auc_markers,
  aes(
    x = cluster,
    y = avg_log2FC,
    color = cluster
  )
) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 1
  ) +
  geom_hline(
    yintercept = fc_cutoff,
    linetype = "dashed",
    color = "red"
  ) +
  geom_hline(
    yintercept = -fc_cutoff,
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = cluster_colors) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Cluster",
    y = "Average log2 fold change",
    title = "High-confidence marker logFC by cluster \nAUC >= 0.6 or <= 0.4, p_adj < 0.05"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# save_plot(
#   plot_object = high_auc_markers_plot,
#   file_stub = "high_auc_markers_plot",
#   save_dir = file.path(plot_dir, "differential_expression"),
#   width = 12,
#   height = 10
# )

#----------------------------------------
# Distribution of marker AUC by cluster
#----------------------------------------

vln_logfc_marker_strength <- ggplot(
  cluster_vs_rest_markers,
  aes(
    x = factor(cluster),
    y = auc,
    fill = factor(cluster)
  )
) +
  geom_violin(
    color = "black",
    trim = FALSE,
    alpha = 0.8
  ) +
  geom_jitter(
    data = dplyr::filter(cluster_vs_rest_markers, !valid_marker),
    color = "grey70",
    width = 0.15,
    height = 0,
    size = 0.8,
    alpha = 0.25
  ) +
  geom_jitter(
    data = all_valid_markers,
    aes(color = avg_log2FC),
    width = 0.15,
    height = 0,
    size = 2.8,
    alpha = 0.95
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1)
  ) +
  geom_hline(
    yintercept = auc_cutoff,
    linetype = "dashed",
    linewidth = 0.6,
    color = "red"
  ) +
  geom_hline(
    yintercept = 1 - auc_cutoff,
    linetype = "dashed",
    linewidth = 0.6,
    color = "red"
  ) +
  scale_fill_manual(
    values = cluster_colors,
    guide = "none"
  ) +
  scale_color_viridis_c(
    option = "inferno",
    name = "Average\nlog2FC"
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Cluster",
    y = "AUC",
    title = "Distribution of marker AUC by cluster"
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )

# save_plot(
#   plot_object = vln_logfc_marker_strength,
#   file_stub = "vln_logfc_marker_strength",
#   save_dir = file.path(plot_dir, "differential_expression"),
#   width = 10,
#   height = 10
# )

#----------------------------------------
# Percentage of valid markers by cluster
#----------------------------------------

#----------------------------------------
# Mean AUC and AUC strength by cluster
#----------------------------------------

#----------------------------------------
# Mean AUC and AUC strength of valid markers by cluster
#----------------------------------------

mean_auc_summary <-
  valid_marker_summary %>%
  dplyr::select(
    cluster,
    mean_auc,
    mean_auc_strength
  ) %>%
  tidyr::pivot_longer(
    cols = c(mean_auc, mean_auc_strength),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("mean_auc", "mean_auc_strength"),
      labels = c("Mean AUC", "Mean AUC Strength")
    )
  )

mean_auc_plot <- ggplot(
  mean_auc_summary,
  aes(
    x = cluster,
    y = value,
    fill = metric
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_text(
    aes(label = sprintf("%.3f", value)),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "Mean AUC" = "#4C78A8",
      "Mean AUC Strength" = "#F58518"
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Cluster",
    y = "Mean value",
    fill = NULL,
    title = "Mean AUC of valid markers by cluster"
  ) +
  theme_classic()
#----------------------------------------
# Combined marker quality figure
#----------------------------------------

combined_marker_summary <-
  (
    vln_cluster_marker_strength |
      bar_cluster_marker_counts
  ) /
  (
    scatter_cluster_marker_quality |
      high_auc_markers_plot
  ) /
  (
    vln_logfc_marker_strength |
      mean_auc_plot
  ) +
  plot_annotation(
    title = "Cluster Marker Summary"
  )

save_plot(
  plot_object = combined_marker_summary,
  file_stub = "combined_marker_summary",
  save_dir = file.path(plot_dir, "differential_expression"),
  width = 18,
  height = 30
)






















# code for dotplot based on dotplot from
# https://romanhaa.github.io/projects/scrnaseq_workflow/

#------------------------------------------------------------
# Top marker gene(s) for each cluster
#------------------------------------------------------------
# If the same gene is selected for multiple clusters,
# combine the cluster IDs into one label (e.g. 0_2_3_Tagln).
#------------------------------------------------------------
# Build unique gene labels
#------------------------------------------------------------

markers_for_plot <-
  top_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  group_by(gene) %>%
  summarise(
    first_cluster = min(as.integer(as.character(cluster))),
    max_logFC = max(avg_log2FC),
    clusters = paste(sort(cluster), collapse = "_"),
    .groups = "drop"
  ) %>%
  arrange(
    first_cluster,
    desc(max_logFC),
    gene
  ) %>%
  mutate(
    gene_label = paste0(clusters, "_", gene)
  )

genes_to_show <- markers_for_plot$gene
gene_labels   <- markers_for_plot$gene_label

# cluster_ids in the order they should appear on the plot.
cluster_ids <- levels(Idents(mega_obj))

# Expression matrix from the assay you want to plot.
expr <- GetAssayData(
  mega_obj,
  layer = "data"
)








#------------------------------------------------------------
# Mean expression per cluster
#------------------------------------------------------------

expression_levels_list <- lapply(cluster_ids, function(cl) {
  
  # Cells belonging to the current cluster
  cells <- WhichCells(
    mega_obj,
    idents = cl
  )
  
  # Subset the expression matrix to the selected genes and cells
  expr_subset <- expr[genes_to_show, cells, drop = FALSE]
  
  # Give the rows unique names so duplicated genes across clusters
  # do not cause problems downstream
  rownames(expr_subset) <- gene_labels
  
  # Mean expression for each gene in this cluster
  means <- Matrix::rowMeans(expr_subset)
  
  # Return a tidy data frame
  data.frame(
    cluster = cl,
    gene_label = names(means),
    expression = as.numeric(means),
    row.names = NULL,
    check.names = FALSE
  )
})

expression_levels <- dplyr::bind_rows(expression_levels_list)

#------------------------------------------------------------
# Percent of cells expressing each gene per cluster
#------------------------------------------------------------

percent_expressing_list <- lapply(cluster_ids, function(cl) {
  
  # Cells belonging to the current cluster
  cells <- WhichCells(
    mega_obj,
    idents = cl
  )
  
  # Subset the expression matrix to the selected genes and cells
  expr_subset <- expr[genes_to_show, cells, drop = FALSE]
  
  # Give the rows unique names so duplicated genes across clusters
  # do not cause problems downstream
  rownames(expr_subset) <- gene_labels
  
  # Fraction of cells expressing each gene (> 0)
  pct <- Matrix::rowMeans(expr_subset > 0)
  
  # Return a tidy data frame
  data.frame(
    cluster = cl,
    gene_label = names(pct),
    percent_cells = as.numeric(pct),
    row.names = NULL,
    check.names = FALSE
  )
})

percent_expressing <- dplyr::bind_rows(percent_expressing_list)

#------------------------------------------------------------
# Join expression and percent-expressing tables
#------------------------------------------------------------

plot_df <-
  dplyr::left_join(
    expression_levels,
    percent_expressing,
    by = c("cluster", "gene_label")
  ) %>%
  dplyr::mutate(
    cluster = factor(cluster, levels = cluster_ids),
    gene_label = factor(gene_label, levels = gene_labels)
  )














#------------------------------------------------------------
# Dot plot
#------------------------------------------------------------

dotplot_top_markers <- ggplot(
  plot_df,
  aes(
    gene_label,
    cluster
  )
) +
  geom_point(
    aes(
      color = expression,
      size = percent_cells
    )
  ) +
  scale_color_viridis_c(
    option = "magma",
    direction = -1,
    name = "Average\nexpression"
  ) +
  scale_size_continuous(
    name = "% expressing",
    labels = scales::percent
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


get_dotplot_dimensions <- function(n_x,
                                   n_y,
                                   cell_size = 0.35,
                                   min_width = 4,
                                   min_height = 4,
                                   pad_x = 3.5,
                                   pad_y = 5) {
  
  list(
    height = max(min_height, pad_y + n_y * cell_size),
    width = max(min_width,  pad_x + n_x * cell_size)
  )
}
dims <- get_dotplot_dimensions(
  n_x = length(levels(plot_df$gene_label)),
  n_y = length(levels(plot_df$cluster))
  
)

save_plot(
  plot_object = dotplot_top_markers,
  file_stub = "dotplot_top_markers",
  save_dir = file.path(plot_dir, "differential_expression"),
  width = dims$width,
  height = dims$height
)











#----------------------------------------------------
# Plot marker gene UMAPs
#----------------------------------------------------

plot_marker_umaps <- function(
    object,
    features = NULL,
    processed_markers = NULL,
    selection = c("manual", "top_auc", "top_auc_strength", "top_logFC"),
    top_n = 20,
    include_negative = FALSE,
    reduction = "full.umap.sketch",
    ncol = 4,
    pt.size = 0.05,
    order = TRUE
) {
  
  selection <- match.arg(selection)
  
  #--------------------------------------------------
  # Automatically select marker genes
  #--------------------------------------------------
  
  if (selection != "manual") {
    
    if (is.null(processed_markers)) {
      stop("processed_markers must be supplied when selection != 'manual'.")
    }
    
    markers <- processed_markers
    
    if (!include_negative) {
      markers <- markers %>%
        dplyr::filter(avg_log2FC > 0)
    }
    
    markers <-
      switch(
        selection,
        
        top_auc =
          markers %>%
          dplyr::arrange(desc(auc)),
        
        top_auc_strength =
          markers %>%
          dplyr::arrange(desc(auc_strength)),
        
        top_logFC =
          markers %>%
          dplyr::arrange(desc(avg_log2FC))
      )
    
    features <-
      markers %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(gene) %>%
      unique()
    
  }
  
  #--------------------------------------------------
  # Check that features exist
  #--------------------------------------------------
  
  if (is.null(features) || length(features) == 0) {
    stop("No features selected.")
  }
  
  features <- intersect(features, rownames(object))
  
  if (length(features) == 0) {
    stop("None of the selected features are present in the object.")
  }
  
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  
  FeaturePlot(
    object = object,
    features = features,
    reduction = reduction,
    combine = TRUE,
    ncol = ncol,
    pt.size = pt.size,
    order = order,
    keep.scale = "feature"
  ) &
    scale_color_viridis_c(
      option = "magma",
      direction = -1
    ) &
    theme_bw() &
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  
}

marker_umaps <- plot_marker_umaps(
  object = mega_obj,
  processed_markers = cluster_vs_rest_markers,
  selection = "top_auc_strength",
  top_n = 16
)

save_plot(
  plot_object = marker_umaps,
  file_stub = "marker_gene_umaps",
  save_dir = file.path(plot_dir, "differential_expression"),
  width = 20,
  height = 16
)









plot_marker_heatmaps <- function(
    object,
    markers,
    expr = GetAssayData(object, layer = "data"),
    cluster_ids = levels(Idents(object)),
    heatmap_dir,
    file_prefix = "marker",
    subtitle = NULL,
    height_scale = 0.26,
    auc_cutoff = 0.60,
    positive_outline = "forestgreen",
    negative_outline = "royalblue",
    tile_border_width = 0.35,
    marker_border_width = 2
) {
  
  #------------------------------------------------------------
  # Build gene list
  #------------------------------------------------------------
  
  markers_for_plot <-
    markers %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    group_by(gene) %>%
    summarise(
      first_cluster = min(as.integer(as.character(cluster))),
      max_logFC = max(avg_log2FC),
      clusters = paste(sort(cluster), collapse = "_"),
      .groups = "drop"
    ) %>%
    arrange(
      first_cluster,
      desc(max_logFC),
      gene
    ) %>%
    mutate(
      gene_label = paste0(clusters, "_", gene)
    )
  
  genes_to_show <- markers_for_plot$gene
  gene_labels   <- markers_for_plot$gene_label
  
  #------------------------------------------------------------
  # Mean expression
  #------------------------------------------------------------
  
  expression_levels_list <- lapply(cluster_ids, function(cl) {
    
    cells <- WhichCells(
      object,
      idents = cl
    )
    
    expr_subset <- expr[
      genes_to_show,
      cells,
      drop = FALSE
    ]
    
    rownames(expr_subset) <- gene_labels
    
    means <- Matrix::rowMeans(expr_subset)
    
    data.frame(
      cluster = cl,
      gene_label = names(means),
      expression = as.numeric(means),
      row.names = NULL,
      check.names = FALSE
    )
    
  })
  
  expression_levels <- bind_rows(expression_levels_list)
  
  #------------------------------------------------------------
  # Percent expressing
  #------------------------------------------------------------
  
  percent_expressing_list <- lapply(cluster_ids, function(cl) {
    
    cells <- WhichCells(
      object,
      idents = cl
    )
    
    expr_subset <- expr[
      genes_to_show,
      cells,
      drop = FALSE
    ]
    
    rownames(expr_subset) <- gene_labels
    
    pct <- Matrix::rowMeans(expr_subset > 0)
    
    data.frame(
      cluster = cl,
      gene_label = names(pct),
      percent_cells = as.numeric(pct),
      row.names = NULL,
      check.names = FALSE
    )
    
  })
  
  percent_expressing <- bind_rows(percent_expressing_list)
  
  #------------------------------------------------------------
  # Build heatmap dataframe
  #------------------------------------------------------------
  
  heatmap_df <-
    expression_levels %>%
    left_join(
      percent_expressing,
      by = c("cluster", "gene_label")
    ) %>%
    group_by(gene_label) %>%
    mutate(
      z_expression = as.numeric(scale(expression)),
      z_percent = as.numeric(scale(percent_cells))
    ) %>%
    ungroup() %>%
    mutate(
      cluster = factor(cluster, levels = cluster_ids),
      gene_label = factor(gene_label, levels = gene_labels)
    )
  

  #------------------------------------------------------------
  # Marker outline annotation
  #------------------------------------------------------------
  
  marker_outline_df <-
    markers %>%
    mutate(
      gene_label = markers_for_plot$gene_label[
        match(gene, markers_for_plot$gene)
      ],
      outline = dplyr::case_when(
        auc >= auc_cutoff ~ "Positive",
        auc <= (1 - auc_cutoff) ~ "Negative",
        TRUE ~ NA_character_
      ),
      auc_label = sprintf("%.2f", auc)
    ) %>%
    dplyr::filter(!is.na(outline)) %>%
    mutate(
      outline = factor(outline, levels = c("Positive", "Negative")),
      cluster = factor(cluster, levels = cluster_ids),
      gene_label = factor(gene_label, levels = gene_labels)
    )
  
  outline_colors <- c(
    "Positive" = "forestgreen",
    "Negative" = "royalblue"
  )
  #------------------------------------------------------------
  # Shared theme
  #------------------------------------------------------------
  
  heatmap_theme <-
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8),
      axis.title = element_blank(),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      legend.title = element_text(face = "bold")
    )
  
  

  #------------------------------------------------------------
  # Heatmaps
  #------------------------------------------------------------
  
#  Average Expression
  p_expression <-
    ggplot(
      heatmap_df,
      aes(cluster, gene_label, fill = expression)
    ) +
    geom_tile(
      color = "white",
      linewidth = tile_border_width
    ) +
    geom_tile(
      data = marker_outline_df,
      aes(cluster, gene_label, color = outline),
      inherit.aes = FALSE,
      fill = NA,
      linewidth = marker_border_width,
      show.legend = TRUE
    ) +
    geom_text(
      data = marker_outline_df,
      aes(cluster, gene_label, label = auc_label),
      inherit.aes = FALSE,
      color = "white",
      fontface = "bold",
      size = 2.6,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = outline_colors,
      name = "Marker type"
    ) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      name = "Average\nExpression"
    ) +
    labs(
      title = "Average Expression",
      subtitle = subtitle
    ) +
    heatmap_theme
  
  
  
  
#  Expression Z-score
  p_expression_z <-
    ggplot(
      heatmap_df,
      aes(cluster, gene_label, fill = z_expression)
    ) +
    geom_tile(
      color = "white",
      linewidth = tile_border_width
    ) +
    geom_tile(
      data = marker_outline_df,
      aes(cluster, gene_label, color = outline),
      inherit.aes = FALSE,
      fill = NA,
      linewidth = marker_border_width,
      show.legend = TRUE
    ) +
    geom_text(
      data = marker_outline_df,
      aes(cluster, gene_label, label = auc_label),
      inherit.aes = FALSE,
      color = "white",
      fontface = "bold",
      size = 2.6,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = outline_colors,
      name = "Marker type"
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "grey",
      high = "red",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = "Expression\nZ-score"
    ) +
    labs(
      title = "Relative Expression (Z-score)",
      subtitle = subtitle
    ) +
    heatmap_theme
  
  
  
#  Percent Expressing
  p_percent <-
    ggplot(
      heatmap_df,
      aes(cluster, gene_label, fill = percent_cells)
    ) +
    geom_tile(
      color = "white",
      linewidth = tile_border_width
    ) +
    geom_tile(
      data = marker_outline_df,
      aes(cluster, gene_label, color = outline),
      inherit.aes = FALSE,
      fill = NA,
      linewidth = marker_border_width,
      show.legend = TRUE
    ) +
    geom_text(
      data = marker_outline_df,
      aes(cluster, gene_label, label = auc_label),
      inherit.aes = FALSE,
      color = "white",
      fontface = "bold",
      size = 2.6,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = outline_colors,
      name = "Marker type"
    ) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      labels = scales::percent,
      name = "%\nExpressing"
    ) +
    labs(
      title = "Percent Expressing",
      subtitle = subtitle
    ) +
    heatmap_theme
# Percent Expressing Z-score
  p_percent_z <-
    ggplot(
      heatmap_df,
      aes(cluster, gene_label, fill = z_percent)
    ) +
    geom_tile(
      color = "white",
      linewidth = tile_border_width
    ) +
    geom_tile(
      data = marker_outline_df,
      aes(cluster, gene_label, color = outline),
      inherit.aes = FALSE,
      fill = NA,
      linewidth = marker_border_width,
      show.legend = TRUE
    ) +
    geom_text(
      data = marker_outline_df,
      aes(cluster, gene_label, label = auc_label),
      inherit.aes = FALSE,
      color = "white",
      fontface = "bold",
      size = 2.6,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = outline_colors,
      name = "Marker type"
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "grey",
      high = "red",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = "% Expr.\nZ-score"
    ) +
    labs(
      title = "Percent Expressing (Z-score)",
      subtitle = subtitle
    ) +
    heatmap_theme
  
  #------------------------------------------------------------
  # Save
  #------------------------------------------------------------
  
  save_plot(
    p_expression,
    paste0(file_prefix, "_heatmap_expression"),
    heatmap_dir,
    width = max(4, length(gene_labels) * height_scale / 2.5),
    height = max(8, length(gene_labels) * height_scale)
  )
  
  save_plot(
    p_expression_z,
    paste0(file_prefix, "_heatmap_expression_zscore"),
    heatmap_dir,
    width = max(4, length(gene_labels) * height_scale / 2.5),
    height = max(8, length(gene_labels) * height_scale)
  )
  
  save_plot(
    p_percent,
    paste0(file_prefix, "_heatmap_percent_expressing"),
    heatmap_dir,
    width = max(4, length(gene_labels) * height_scale / 2.5),
    height = max(8, length(gene_labels) * height_scale)
  )
  
  save_plot(
    p_percent_z,
    paste0(file_prefix, "_heatmap_percent_expressing_zscore"),
    heatmap_dir,
    width = max(4, length(gene_labels) * height_scale / 2.5),
    height = max(8, length(gene_labels) * height_scale)
  )
  
  invisible(
    list(
      expression = p_expression,
      expression_z = p_expression_z,
      percent = p_percent,
      percent_z = p_percent_z,
      data = heatmap_df
    )
  )
  
}




plot_marker_heatmaps(
  object = mega_obj,
  markers = top_markers,
  heatmap_dir = heatmap_dir,
  file_prefix = "top_marker",
  subtitle = "Top marker genes passing differential expression thresholds",
  height_scale = 0.22
)

plot_marker_heatmaps(
  object = mega_obj,
  markers = all_valid_markers,
  heatmap_dir = heatmap_dir,
  file_prefix = "all_marker",
  subtitle = "All marker genes passing differential expression thresholds",
  height_scale = 0.25
)





























