source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")
library(clustree)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(grid)


options(future.globals.maxSize = 1000 * 1024^2)

#start by loading in object

sketch_label <- "sketch_33pct"



# Decide on embedding
max_dims <- 11

pca_dir <- get_pca_dir(
  root_dir = mega_dir,
  ndims = 50,
  sketch_label = sketch_label
)

#   mega_obj <- load_or_build_object(  #would need if not using harmony for analysis
#   directory = pca_dir,
#   build_function = function() {
#     run_sketch_pca(
#       object = mega_obj,
#       ncells = ncells_sketch
#     )
#   },
#   force_rebuild = FALSE
# )

if (embedding == "harmony") {
  
  group.by.vars <- "tissue"
  
  harmony_dir <- get_harmony_dir(
    root_dir = mega_dir,
    ndims = max_dims,
    sketch_label = "sketch_33pct",
    group.by.vars = group.by.vars
  )
  
  mega_obj <- load_or_build_object(
    directory = harmony_dir,
    build_function = function() {
      run_harmony_embedding(
        object = mega_obj,
        ndims = max_dims,
        sketch_reduction = "pca.sketch",
        group.by.vars = "tissue"
      )
    },
    force_rebuild = FALSE
  )

  embedding_run_dir <- harmony_dir
  
} else {
  
  
  embedding_run_dir <- pca_dir
  
}



cluster_resolutions <- c(0.3,0.4,0.5,0.6,0.7,0.8,1,1.2)

selected_cluster_resolution <- 0.3



clustering_dir_normal <- get_clustering_dir(
  embedding_run_dir = embedding_run_dir,
  ndims = max_dims,
  resolution = selected_cluster_resolution
)


clustering_dir <- get_nested_subdir(
  clustering_dir_normal,
  paste0("multi_res_",selected_cluster_resolution)
)

mega_obj <- load_or_build_object(
  directory = clustering_dir,
  build_function = function() {
    
    mega_obj <- run_sketch_clustering(
      object = mega_obj,
      sketch_reduction = embedding_reduction,
      dims = 1:max_dims,
      resolutions = cluster_resolutions,
      selected_cluster_resolution = selected_cluster_resolution
    )
    
  },
  force_rebuild = FALSE
)


#takes only cells used in the sketch
cells <- WhichCells(
  mega_obj,
  expression = !is.na(seurat_cluster.sketched)
)


#takes only cells used in the sketch
mega_obj_sketch <- subset(mega_obj, cells = cells)


colnames(mega_obj_sketch@meta.data)

# moves cluster resolutions to the end for better clarity
cluster_cols <- grep(
  "^cluster_res_",
  colnames(mega_obj_sketch@meta.data),
  value = TRUE
)
mega_obj_sketch@meta.data <- mega_obj_sketch@meta.data[, c(
  setdiff(colnames(mega_obj_sketch@meta.data), cluster_cols),
  cluster_cols
)]

colnames(mega_obj_sketch@meta.data)


# citation("clustree")
# 
# Zappia L, Oshlack A. Clustering trees: a visualization for
# evaluating clusterings at multiple resolutions. GigaScience.
# 2018;7. DOI:gigascience/giy083
#
# https://lazappi.github.io/clustree/


# Saving/loading plot directory
plot_dir <- get_figures_dir(paste0(clustering_dir,"/clustree"))








# Basic clustering tree
p_clustree <- clustree(
  mega_obj_sketch,
  prefix = "cluster_res_",
  layout = "sugiyama"
)

save_plot(
  plot_object = p_clustree,
  file_stub = "clustree_basic",
  save_dir = plot_dir,
  width = 12,
  height = 8
)

# Mean leverage score
p_clustree_leverage <- clustree(
  mega_obj_sketch,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_colour = "leverage.score",
  node_colour_aggr = "mean"
)

save_plot(
  plot_object = p_clustree_leverage,
  file_stub = "clustree_leverage_score",
  save_dir = plot_dir,
  width = 12,
  height = 8
)

# Mean percent mitochondrial
p_clustree_percent_mt <- clustree(
  mega_obj_sketch,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_colour = "percent.mt",
  node_colour_aggr = "mean"
)
save_plot(
  plot_object = p_clustree_percent_mt,
  file_stub = "clustree_percent_mt",
  save_dir = plot_dir,
  width = 12,
  height = 8
)

# Mean UMI count
p_clustree_ncount <- clustree(
  mega_obj_sketch,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_colour = "nCount_Spatial.008um",
  node_colour_aggr = "mean"
)
save_plot(
  plot_object = p_clustree_ncount,
  file_stub = "clustree_nCount",
  save_dir = plot_dir,
  width = 12,
  height = 8
)

# Mean detected features
p_clustree_nfeature <- clustree(
  mega_obj_sketch,
  prefix = "cluster_res_",
  layout = "sugiyama",
  node_colour = "nFeature_Spatial.008um",
  node_colour_aggr = "mean"
)
save_plot(
  plot_object = p_clustree_nfeature,
  file_stub = "clustree_nFeature",
  save_dir = plot_dir,
  width = 12,
  height = 8
)




# Promote selected resolution to canonical downstream column


selected_cluster_col <- "cluster_res_0.3" 

if (!selected_cluster_col %in% colnames(mega_obj_sketch[[]])) {
  stop(
    "Selected cluster column not found: ",
    selected_cluster_col
  )
}
mega_obj_sketch[["seurat_cluster.sketched"]] <- mega_obj_sketch[[selected_cluster_col]]
Idents(mega_obj_sketch) <- "seurat_cluster.sketched"







#--------------------------------------------------
# computes silhouette scores for all clusters
#---------------------------------------------------
cluster_cols <- grep(
  "^cluster_res_",
  colnames(mega_obj_sketch[[]]),
  value = TRUE
)

for (cluster_col in cluster_cols) {
  
  message("Computing silhouette for ", cluster_col)
  
  mega_obj_sketch <- add_silhouette(
    object = mega_obj_sketch,
    reduction = embedding_reduction,
    dims = 1:max_dims,
    cluster_column = cluster_col,
    metadata_column = paste0("silhouette_", cluster_col)
  )
  
}


#--------------------------------------------------
# silhouette clustree
#---------------------------------------------------
# silhoutette_clustree <- clustree(
#   mega_obj_sketch,
#   prefix = "cluster_res_",
#   layout = "sugiyama",
#   node_colour = "silhouette_cluster_res_0.5",
#   node_colour_aggr = "mean"
# )
# 
# save_plot(
#   plot_object = silhoutette_clustree,
#   file_stub = "silhoutette_clustree",
#   save_dir = plot_dir,
#   width = 12,
#   height = 8
# )









#creates massive resolution analysis plot
plots <- vector("list", length(cluster_cols))
names(plots) <- cluster_cols

resolution_summary <- data.frame(
  resolution = numeric(length(cluster_cols)),
  mean_silhouette = numeric(length(cluster_cols)),
  median_silhouette = numeric(length(cluster_cols)),
  n_clusters = integer(length(cluster_cols)),
  smallest_cluster = integer(length(cluster_cols)),
  projection_confidence = numeric(length(cluster_cols))
)

for (i in seq_along(cluster_cols)) {
  
  cluster_col <- cluster_cols[i]
  
  silhouette_col <- sub(
    "^cluster_",
    "silhouette_cluster_",
    cluster_col
  )
  
  # ------------------------------------------
  # Subset to cells with silhouette values
  # ------------------------------------------
  
  cells <- colnames(mega_obj_sketch)[
    !is.na(mega_obj_sketch[[silhouette_col]][, 1])
  ]
  
   plot_obj <- subset(
     mega_obj_sketch,
     cells = cells
   )
   
   sil <- plot_obj[[silhouette_col]][,1]
  # ------------------------------------------
  # Cluster ordering
  # ------------------------------------------
  
  cluster_counts <- sort(
    table(plot_obj[[cluster_col]][,1]),
    decreasing = TRUE
  )
  
  cluster_order <- names(cluster_counts)
  
  cluster_percent <-
    100 * cluster_counts / sum(cluster_counts)
  
  cluster_label_map <- setNames(
    paste0(
      cluster_order,
      "\n",
      sprintf("%.1f%%", cluster_percent),
      "\n",
      format(cluster_counts, big.mark = ",")
    ),
    cluster_order
  )
  
  plot_obj$cluster_ordered <- factor(
    plot_obj[[cluster_col]][,1],
    levels = cluster_order
  )
  
  #adding info to resolution summary
  resolution_summary$resolution[i] <-
    as.numeric(sub("cluster_res_", "", cluster_col))
  
  resolution_summary$mean_silhouette[i] <-
    mean(sil, na.rm = TRUE)
  
  resolution_summary$median_silhouette[i] <-
    median(sil, na.rm = TRUE)
  
  resolution_summary$n_clusters[i] <-
    length(cluster_counts)
  
  resolution_summary$smallest_cluster[i] <-
    min(cluster_counts)
  
  
  resolution_summary$percent_negative[i] <-
    mean(sil < 0, na.rm = TRUE)
  
  resolution_summary$percent_good[i] <-
    mean(sil > 0.5, na.rm = TRUE)
  
  # ------------------------------------------
  # Violin plot
  # ------------------------------------------
  
  p_violin <- VlnPlot(
    plot_obj,
    features = silhouette_col,
    group.by = "cluster_ordered",
    pt.size = 0
  ) +
    scale_x_discrete(labels = cluster_label_map) +
    labs(
      title = paste(
        "Resolution",
        sub("cluster_res_", "", cluster_col)
      ),
      x = "Cluster\n(% sketch / count)",
      y = "Approximate silhouette"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 8
      ),
      axis.title.x = element_text(margin = margin(t = 10)),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      )
    )
  
  # ------------------------------------------
  # Histogram of silhouette values
  # ------------------------------------------
  
  hist_df <- data.frame(
    silhouette = sil
  ) %>%
    mutate(
      quality = cut(
        silhouette,
        breaks = c(-Inf, 0, 0.25, 0.5, 0.75, Inf),
        labels = c(
          "< 0",
          "0 - 0.25",
          "0.25 - 0.50",
          "0.50 - 0.75",
          "> 0.75"
        ),
        right = FALSE
      )
    )
  
  p_hist <- ggplot(
    hist_df,
    aes(
      x = silhouette,
      fill = quality
    )
  ) +
    geom_histogram(
      breaks = seq(-1, 1, by = 0.05),
      color = "black",
      linewidth = 0.2
    ) +
    scale_fill_manual(
      values = c(
        "< 0"         = "#d73027",  # red
        "0 - 0.25"    = "#fc8d59",  # orange
        "0.25 - 0.50" = "#fee08b",  # yellow
        "0.50 - 0.75" = "#91cf60",  # light green
        "> 0.75"      = "#1a9850"   # dark green
      ),
      drop = FALSE,
      name = "Silhouette"
    ) +
    geom_vline(
      xintercept = mean(sil, na.rm = TRUE),
      color = "blue",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = 0,
      color = "red4",
      linewidth = 0.6
    ) +
    coord_cartesian(
      xlim = c(-1, 1)
    ) +
    labs(
      x = "Silhouette Score",
      y = "Cells"
    ) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  # ------------------------------------------
  # Summary statistics
  # ------------------------------------------
  

  
  stats_df <- data.frame(
    Metric = c(
      "Clusters",
      "Mean",
      "Median",
      "Min",
      "Max",
      "% < 0",
      "% > 0.50",
      "% > 0.75",
      "Smallest"
    ),
    Value = c(
      length(cluster_counts),
      sprintf("%.3f", mean(sil, na.rm = TRUE)),
      sprintf("%.3f", median(sil, na.rm = TRUE)),
      sprintf("%.3f", min(sil, na.rm = TRUE)),
      sprintf("%.3f", max(sil, na.rm = TRUE)),
      sprintf("%.1f%%", 100 * mean(sil < 0, na.rm = TRUE)),
      sprintf("%.1f%%", 100 * mean(sil > 0.50, na.rm = TRUE)),
      sprintf("%.1f%%", 100 * mean(sil > 0.75, na.rm = TRUE)),
      min(cluster_counts)
    )
  )
  
  table_grob <- tableGrob(
    stats_df,
    rows = NULL,
    theme = ttheme_minimal(
      base_size = 8
    )
  )
  
  p_table <- wrap_elements(table_grob)
  
  # ------------------------------------------
  # Combine
  # ------------------------------------------
  
  plots[[i]] <-
    (p_violin / (p_hist | p_table)) +
    plot_layout(
      heights = c(3, 2),
      widths = c(3, 1)
    )
  
} #loop ends


# Theme for plots of silhouette changes over resolution
summary_theme <-
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 8, 2, 8)
  )

chosen_resolution_line <-
  geom_vline(
    xintercept = selected_cluster_resolution,
    linetype = "dashed",
    color = "red",
    linewidth = 0.6
  )

# Helper for centered panel labels
panel_label <- function(label, y) {
  annotate(
    "text",
    x = mean(range(resolution_summary$resolution)),
    y = y,
    label = label,
    hjust = 0.5,
    vjust = 1,
    fontface = "bold",
    size = 3.5
  )
}

# Mean silhouette
p_mean <-
  ggplot(
    resolution_summary,
    aes(resolution, mean_silhouette)
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = sprintf("%.2f", mean_silhouette)),
    vjust = -0.8,
    size = 3
  ) +
  scale_x_continuous(
    breaks = resolution_summary$resolution
  ) +
  panel_label(
    "Mean Silhouette",
    max(resolution_summary$mean_silhouette, na.rm = TRUE)
  ) +
  chosen_resolution_line +
  summary_theme

# Median silhouette
p_median <-
  ggplot(
    resolution_summary,
    aes(resolution, median_silhouette)
  ) +
  geom_line(linewidth = 0.8, color = "#2C7BB6") +
  geom_point(size = 2.8, color = "#2C7BB6") +
  geom_text(
    aes(label = sprintf("%.2f", median_silhouette)),
    vjust = -0.8,
    size = 3
  ) +
  scale_x_continuous(
    breaks = resolution_summary$resolution
  ) +
  panel_label(
    "Median Silhouette",
    max(resolution_summary$median_silhouette, na.rm = TRUE)
  ) +
  chosen_resolution_line +
  summary_theme

# Percent negative
p_negative <-
  ggplot(
    resolution_summary,
    aes(resolution, percent_negative)
  ) +
  geom_line(linewidth = 0.8, color = "#D73027") +
  geom_point(size = 2.8, color = "#D73027") +
  geom_text(
    aes(label = scales::percent(percent_negative, accuracy = 1)),
    vjust = -0.8,
    size = 3
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  scale_x_continuous(
    breaks = resolution_summary$resolution
  ) +
  panel_label(
    "Negative Silhouettes",
    max(resolution_summary$percent_negative, na.rm = TRUE)
  ) +
  chosen_resolution_line +
  summary_theme

# Number of clusters
p_clusters <-
  ggplot(
    resolution_summary,
    aes(resolution, n_clusters)
  ) +
  geom_line(linewidth = 0.8, color = "#238443") +
  geom_point(size = 2.8, color = "#238443") +
  geom_text(
    aes(label = n_clusters),
    vjust = -0.8,
    size = 3
  ) +
  scale_x_continuous(
    breaks = resolution_summary$resolution
  ) +
  panel_label(
    "Number of Clusters",
    max(resolution_summary$n_clusters, na.rm = TRUE)
  ) +
  chosen_resolution_line +
  summary_theme

# Smallest cluster
p_smallest <-
  ggplot(
    resolution_summary,
    aes(resolution, smallest_cluster)
  ) +
  geom_line(linewidth = 0.8, color = "#7B3294") +
  geom_point(size = 2.8, color = "#7B3294") +
  geom_text(
    aes(label = smallest_cluster),
    vjust = -0.8,
    size = 3
  ) +
  scale_x_continuous(
    breaks = resolution_summary$resolution
  ) +
  panel_label(
    "Smallest Cluster",
    max(resolution_summary$smallest_cluster, na.rm = TRUE)
  ) +
  chosen_resolution_line +
  summary_theme +
  labs(
    x = "Resolution"
  ) +
  theme(
    axis.text.x = element_text(),
    axis.ticks.x = element_line(),
    axis.title.x = element_text()
  )



#Combine the summary dashboard
p_resolution_summary <-
  p_mean /
  p_median /
  p_negative /
  p_clusters /
  p_smallest 


n_rows <- ceiling(length(plots) / 3)



resolution_grid <-
  wrap_plots(
    plots,
    ncol = 2,
    byrow = TRUE
  )

silhouette_summary_plot <-
  wrap_plots(
    p_resolution_summary,
    plot_spacer(),
    resolution_grid,
    ncol = 3,
    widths = c(1.2, 0.2, 8)
  )


save_plot(
  plot_object = silhouette_summary_plot,
  file_stub = "silhouette_summary",
  save_dir = plot_dir,
  width = 26,
  height = 8 * n_rows
)






# TODO maybe make a version of this that uses the projected object? add these umaps
# 
# 
# 
# #----------------------------------------
# # Categorize silhouette scores
# #----------------------------------------
# 
# mega_obj$silhouette_category <- cut(
#   mega_obj$silhouette,
#   breaks = c(-Inf, 0, 0.25, 0.50, 0.75, Inf),
#   labels = c(
#     "< 0",
#     "0 - 0.25",
#     "0.25 - 0.50",
#     "0.50 - 0.75",
#     "> 0.75"
#   ),
#   include.lowest = TRUE,
#   right = FALSE
# )
# 
# #----------------------------------------
# # UMAP colored by silhouette category
# #----------------------------------------
# 
# silhouette_umap_categorical <-
#   DimPlot(
#     mega_obj,
#     reduction = "umap.sketch",
#     group.by = "silhouette_category",
#     pt.size = 0.01,
#     shuffle = TRUE,
#     raster = FALSE
#   ) +
#   scale_color_manual(
#     values = c(
#       "< 0"         = "#d73027",
#       "0 - 0.25"    = "#fc8d59",
#       "0.25 - 0.50" = "#fee08b",
#       "0.50 - 0.75" = "#91cf60",
#       "> 0.75"      = "#1a9850"
#     ),
#     drop = FALSE,
#     name = "Silhouette"
#   ) +
#   ggtitle("Approximate Sketch Silhouette") +
#   theme_classic()
# 
# 
# save_plot(
#   plot_object = silhouette_umap_categorical,
#   file_stub = paste0("silhouette_umap_categorical"),
#   save_dir = plot_dir,
#   width = 12,
#   height = 10,
#   dpi = 1200
# )
# 
# #combine and save
# silhouette_umaps <- umap_silhouette_plot_red_green |
#   silhouette_umap_categorical
# 
# save_plot(
#   plot_object = silhouette_umaps,
#   file_stub = paste0("combined_silhouette_umaps"),
#   save_dir = plot_dir,
#   width = 20,
#   height = 8,
#   dpi = 1200
# )
# 
# 
