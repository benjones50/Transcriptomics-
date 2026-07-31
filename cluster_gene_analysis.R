source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")
plot_dir <- paste0(get_figures_dir(projection_dir),"/cluster_gene_analysis")

#used to look at michalis gene set for initial clustering 

#----------------------------------------
# Define module gene sets
#----------------------------------------

module_gene_sets <- list(
  Urothelium = c("Krt5", "Krt8", "Upk3a"),
  Stroma_fibroblasts = c("Col1a1", "Col1a2", "Dcn"),
  SmoothMuscle = c("Acta2", "Myh11", "Tagln"),
  Endothelium_vessels = c("Pecam1", "Kdr", "Cdh5"),
  Proliferation = c("Mki67", "Top2a"),
  Myeloid_cells = c("Adgre1", "Cd68", "Csf1r"),
  T_cells = c("Cd3e", "Cd8a", "Trac"),
  NK_cytotoxic_lymphocytes = c("Nkg7", "Gzmb", "Prf1")
)

#----------------------------------------
# Compute module scores
#----------------------------------------

for (module_name in names(module_gene_sets)) {
  
  message(paste0(module_name," scored"))
  mega_obj <- AddModuleScore(
    mega_obj,
    features = list(module_gene_sets[[module_name]]),
    name = module_name
  )
  
}

module_features <- paste0(names(module_gene_sets), "1")

#----------------------------------------
# Overall module score UMAP
#----------------------------------------

p_module_scores_umap <- FeaturePlot(
  mega_obj,
  features = module_features,
  combine = TRUE
)

save_plot(
  plot_object = p_module_scores_umap,
  file_stub = "overview_module_scores_umap",
  save_dir = plot_dir,
  width = 14,
  height = 10
)

#----------------------------------------
# Module score + individual gene UMAPs
#----------------------------------------

for (module_name in names(module_gene_sets)) {
  
  genes <- module_gene_sets[[module_name]]
  
  module_feature <- paste0(module_name, "1")
  
  features_to_plot <- c(module_feature, genes)
  
  feature_titles <- c(
    "Module Score",
    genes
  )
  
  p <- FeaturePlot(
    mega_obj,
    features = features_to_plot,
    combine = TRUE,
    ncol = length(features_to_plot)
  ) &
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
  
  # Replace panel titles with nicer labels
  if (inherits(p, "patchwork")) {
    for (i in seq_along(p$patches$plots)) {
      p$patches$plots[[i]] <-
        p$patches$plots[[i]] +
        ggtitle(feature_titles[i])
    }
  }
  
  p <-
    p +
    patchwork::plot_annotation(
      title = paste(module_name, "Module"),
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 16,
          hjust = 0.5
        )
      )
    )
  
  save_plot(
    plot_object = p,
    file_stub = paste0("umap_module_", module_name),
    save_dir = file.path(plot_dir, "module_genes"),
    width = 5 * length(features_to_plot),
    height = 5
  )
  
}

#----------------------------------------------------
# Helper: make one UMAP panel for a single feature
#----------------------------------------------------
make_feature_umap_panel <- function(
    object,
    feature,
    title,
    reduction = "full.umap.sketch",
    pt.size = 0.01,
    legend.position = "right"
) {
  
  p <- FeaturePlot(
    object,
    features = feature,
    reduction = reduction,
    combine = FALSE,
    order = TRUE,
    pt.size = pt.size,
    raster = FALSE
  )[[1]] +
    scale_color_gradientn(
      colours = c(
        "#440154",
        "#3B528B",
        "#21918C",
        "#5DC863",
        "#FDE725"
      )
    ) +
    ggtitle(title) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = legend.position,
      legend.title = element_text(face = "bold")
    )
  
  p
}






tissue_order <- unique(mega_obj[["slide_id"]][, 1])


#----------------------------------------
# Spatial module score + individual gene plots (vertical)
#----------------------------------------
  make_spatial_feature_grid <- function(
    object,
    features,
    row_titles = features,
    image.alpha = 0.05,
    pt.size.factor = 5,
    legend.position = "right",
    groupings = "slide_id" 
  ) {

    stopifnot(length(features) == length(row_titles))

    rows <- vector("list", length(features))



    # Maximum expression across all genes (not module score)
    # used for histograms
    gene_max <- ceiling(
      max(
        FetchData(
          object,
          vars = unlist(module_gene_sets)
        ),
        na.rm = TRUE
      )
    )
    

    
    for (i in seq_along(features)) {

      #----------------------------------------------------
      # Spatial plots
      #----------------------------------------------------
      row_plots <- SpatialFeaturePlot(
        object,
        features = features[i],
        image.alpha = image.alpha,
        pt.size.factor = pt.size.factor,
        combine = FALSE
      )

      row_plots[[1]] <-
        row_plots[[1]] +
        ggtitle(row_titles[i]) +
        theme(
          plot.title = element_text(
            face = "bold",
            hjust = 0.5
          )
        )

      row_plots <- lapply(
        row_plots,
        function(p) {
          p +
            theme(
              legend.position = legend.position,
              legend.title = element_text(face = "bold")
            )
        }
      )
      
      #add legend to last plot
      for (j in seq_along(row_plots)) {
        
        row_plots[[j]] <-
          row_plots[[j]] +
          theme(
            legend.position =
              if (j == length(row_plots))
                "right"
            else
              "none"
          )
        
      }
      spatial_row <- wrap_plots(
        row_plots,
        nrow = 1
      )
      
      
      #----------------------------------------------------
      # Far-left UMAP panel for this feature
      #----------------------------------------------------
      umap_plot <- make_feature_umap_panel(
        object = object,
        feature = features[i],
        title = row_titles[i],
        pt.size = 0.01,
        legend.position = legend.position
      )
      
      
      
      #----------------------------------------------------
      # Matching histogram panel
      #----------------------------------------------------
      plot_df <- FetchData(
        object,
        vars = c(features[i], groupings)
      )

      colnames(plot_df) <- c("expression", "group")
      
      # #----------------------------------------------------
      # # DEBUG: verify expression values for each feature
      # #----------------------------------------------------
      # cat("\n========================================\n")
      # cat("Feature:", features[i], "\n")
      # cat("Class:", class(plot_df$expression), "\n")
      # cat("Default assay:", DefaultAssay(object), "\n")
      # cat("Feature in metadata?:", features[i] %in% colnames(object[[]]), "\n")
      # 
      # cat("\nExpression summary:\n")
      # print(summary(plot_df$expression))
      # 
      # cat("\nUnique values (first 20):\n")
      # print(head(sort(unique(plot_df$expression)), 20))
      # 
      # cat("\nPercent expressing by tissue:\n")
      # print(
      #   plot_df |>
      #     dplyr::group_by(group) |>
      #     dplyr::summarise(
      #       pct = mean(expression > 0) * 100,
      #       n_positive = sum(expression > 0),
      #       n_total = dplyr::n(),
      #       .groups = "drop"
      #     )
      # )
      # 
      # cat("========================================\n")
      # Use a single consistent tissue ordering everywhere
      plot_df$group <- factor(
        plot_df$group,
        levels = tissue_order
      )

      #----------------------------------------------------
      # Summary statistics
      #----------------------------------------------------
      if (i == 1) {

        # Module score: median score per tissue
        summary_df <-
          plot_df |>
          dplyr::group_by(group) |>
          dplyr::summarise(
            value = median(expression, na.rm = TRUE),
            .groups = "drop"
          )

        summary_df$label <- sprintf("%.2f", summary_df$value)

      } else {

        # Gene: percent of bins expressing the gene
        summary_df <-
          plot_df |>
          dplyr::group_by(group) |>
          dplyr::summarise(
            value = mean(expression > 0) * 100,
            .groups = "drop"
          )

        summary_df$label <- sprintf("%.0f%%", summary_df$value)

        # Remove zero-expression bins for histogram only
        plot_df <-
          dplyr::filter(
            plot_df,
            expression > 0
          )
      }

      # Keep ordering identical
      summary_df$group <- factor(
        summary_df$group,
        levels = rev(tissue_order)
      )

      #----------------------------------------------------
      # Histogram
      #----------------------------------------------------
      hist_plot <-
        ggplot(
          plot_df,
          aes(x = expression)
        ) +
        geom_histogram(
          binwidth = 0.25,
          boundary = 0,
          color = "white",
          linewidth = 0.2
        ) +
        facet_wrap(
          ~group,
          ncol = 1,
          scales = "fixed"
        ) +
        theme_classic(base_size = 9) +
        theme(
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          plot.margin = margin(5, 5, 5, 5)
        )

      # Add column header only to first row
      if (i == 1) {
        hist_plot <-
          hist_plot +
          labs(title = "Distribution") +
          theme(
            plot.title = element_text(
              face = "bold",
              hjust = 0.5,
              size = 11
            )
          )
      }

      # Give all genes the same x-axis
      # Leave the module score automatic
      if (i != 1) {
        hist_plot <-
          hist_plot +
          labs(title = "Distribution, 0's removed") +
          theme(
            plot.title = element_text(
              face = "bold",
              hjust = 0.5,
              size = 11
            )
            )+
          coord_cartesian(
            xlim = c(0, gene_max)
          ) +
          scale_x_continuous(
            breaks = pretty(c(0, gene_max), n = 5),
            expand = expansion(mult = c(0, 0.02))
          )

      }

      #----------------------------------------------------
      # Summary panel
      #----------------------------------------------------
      if (i == 1) {

        summary_plot <-
          ggplot(
            summary_df,
            aes(
              x = value,
              y = group,
              fill = value
            )
          ) +
          geom_col(width = 0.7) +
          scale_fill_gradient(
            low = "#DCEAF7",
            high = "#2C7FB8",
            guide = "none"
          )

      } else {

        summary_plot <-
          ggplot(
            summary_df,
            aes(
              x = value,
              y = group
            )
          ) +
          geom_col(
            width = 0.7,
            fill = "grey45"
          )

      }

      # Add common layers
      summary_plot <-
        summary_plot +
        geom_text(
          aes(label = label),
          hjust = -0.15,
          size = 2.5
        ) +
        theme_classic(base_size = 9) +
        theme(
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(5, 5, 5, 5)
        )

      # Column headers
      if (i == 1) {

        summary_plot <-
          summary_plot +
          labs(title = "Median\nModule Score") +
          theme(
            plot.title = element_text(
              face = "bold",
              hjust = 0.5,
              size = 11
            )
          ) +
          coord_cartesian(clip = "off") +
          scale_x_continuous(
            breaks = pretty(summary_df$value, n = 4),
            expand = expansion(mult = c(0, 0.15))
          ) +
          theme(
            plot.margin = margin(5, 20, 5, 5)
          )

      } else {

        summary_plot <-
          summary_plot +
          labs(title = "% Bins\nExpressing") +
          coord_cartesian(clip = "off") +
          scale_x_continuous(
            limits = c(0, 115),
            breaks = c(0, 50, 100),
            expand = expansion(mult = c(0, 0.05))
          ) +
          theme(
            plot.margin = margin(5, 20, 5, 5)
          )

      }


      #----------------------------------------------------
      # Combine UMAP + spatial row + histogram + summary
      #----------------------------------------------------
      rows[[i]] <- patchwork::wrap_plots(
        umap_plot,
        spatial_row,
        hist_plot,
        summary_plot,
        widths = c(2.1, 8, 2, 1.3),
        nrow = 1
      )
    }



    #----------------------------------------------------
    # Bottom QC plot
    #----------------------------------------------------

    tissue_counts <- table(object[[groupings]][, 1])

    tissue_label_map <- setNames(
      paste0(
        names(tissue_counts),
        "\n(n = ",
        format(tissue_counts, big.mark = ","),
        ")"
      ),
      names(tissue_counts)
    )

    qc_df <- FetchData(
      object,
      vars = c("nCount_Spatial.008um", groupings)
    )

    colnames(qc_df) <- c("UMIs", "Tissue")

    qc_df$Tissue <- factor(
      qc_df$Tissue,
      levels = tissue_order
    )


    #vln plot of counts
    counts_plot <-
      ggplot(
        qc_df,
        aes(
          Tissue,
          UMIs,
          fill = Tissue
        )
      ) +
      geom_violin(
        scale = "width",
        trim = FALSE
      ) +
      labs(
        title = "UMI Counts by Tissue",
        x = "Tissue",
        y = "UMIs"
      ) +
      theme_classic(base_size = 11) +
      theme(
        legend.position = "none",
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 8
        )
      )
    
    
    
    #----------------------------------------------------
    # Summary statistics table
    #----------------------------------------------------
    
    table_df <- FetchData(
      object,
      vars = c(
        module_feature,
        genes,
        "nCount_Spatial.008um",
        groupings
      )
    )
    
    colnames(table_df) <-
      c(
        "ModuleScore",
        genes,
        "UMIs",
        "Tissue"
      )
    #consistent ordering 
    table_df$Tissue <- factor(
      table_df$Tissue,
      levels = tissue_order
    )
    
    table_df <-
      table_df |>
      dplyr::mutate(
        genes_detected =
          rowSums(
            dplyr::across(all_of(genes)) > 0
          ),
        frac_detected =
          genes_detected / length(genes),
        triple_positive =
          genes_detected == length(genes)
      ) |>
      dplyr::group_by(Tissue) |>
      dplyr::summarise(
        `Median Module Score` = median(ModuleScore),
        `Mean Genes Detected` = mean(genes_detected),
        `Mean UMIs` = mean(UMIs),
        `Median UMIs` = median(UMIs),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        across(where(is.numeric), ~round(.x, 2))
      )
    table_df <-
      table_df |>
      dplyr::arrange(Tissue)
    
    
    table_plot <-
      patchwork::wrap_elements(
        gridExtra::tableGrob(
          table_df,
          rows = NULL,
          theme = gridExtra::ttheme_minimal(
            base_size = 9
          )
        )
      )
    
    

    #----------------------------------------------------
    # Bottom row aligned to the 4-column layout
    #----------------------------------------------------
    bottom_row <-
      patchwork::plot_spacer() |
      counts_plot |
      patchwork::plot_spacer() |
      table_plot +
      plot_layout(
        widths = c(2.1, 8, 2, 1.3)
      )

    wrap_plots(
      c(rows, list(bottom_row)),
      ncol = 1,
      heights = c(rep(6, length(rows)), 3)
    )

  }


# used for testing
#----------------------------------------
# Create and save ONLY the Smooth Muscle figure
#----------------------------------------

module_name <- "SmoothMuscle"

genes <- module_gene_sets[[module_name]]
module_feature <- paste0(module_name, "1")

module_title <- dplyr::recode(
  module_feature,
  "Urothelium1" = "Urothelium",
  "Stroma_fibroblasts1" = "Fibroblasts",
  "SmoothMuscle1" = "Smooth Muscle",
  "Endothelium_vessels1" = "Endothelium",
  "Proliferation1" = "Proliferation",
  "Myeloid_cells1" = "Myeloid Cells",
  "T_cells1" = "T Cells",
  "NK_cytotoxic_lymphocytes1" = "NK Cells",
  "Bach11" = "Bach1",
  "Hmox11" = "Hmox1"
)

p <- make_spatial_feature_grid(
  object = mega_obj,
  features = c(module_feature, genes),
  row_titles = c("Module Score", genes),
  image.alpha = 0.1,
  pt.size.factor = 5
) +
  plot_annotation(
    title = paste(module_title, "Module"),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      )
    )
  )

save_plot(
  plot_object = p,
  file_stub = paste0("spatial_module_", module_name),
  save_dir = file.path(plot_dir, "spatial_module_genes"),
  width = 24,
  height = 4 * (length(genes) + 1) + 3
)














#TODO remove, used for testing
print("stop")
print("stop")
print("stop")
print("stop")
print("stop")
print("stop")

#----------------------------------------
# Create and save module + gene figures
#----------------------------------------

for (module_name in names(module_gene_sets)) {
  
  genes <- module_gene_sets[[module_name]]
  module_feature <- paste0(module_name, "1")
  
  module_title <- dplyr::recode(
    module_feature,
    "Urothelium1" = "Urothelium",
    "Stroma_fibroblasts1" = "Fibroblasts",
    "SmoothMuscle1" = "Smooth Muscle",
    "Endothelium_vessels1" = "Endothelium",
    "Proliferation1" = "Proliferation",
    "Myeloid_cells1" = "Myeloid Cells",
    "T_cells1" = "T Cells",
    "NK_cytotoxic_lymphocytes1" = "NK Cells",
    "Bach11" = "Bach1",
    "Hmox11" = "Hmox1"
  )
  
  p <- make_spatial_feature_grid(
    object = mega_obj,
    features = c(module_feature, genes),
    row_titles = c("Module Score", genes),
    image.alpha = 0.05,
    pt.size.factor = 3
  ) +
    plot_annotation(
      title = paste(module_title, "Module"),
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 16,
          hjust = 0.5
        )
      )
    )
  
  save_plot(
    plot_object = p,
    file_stub = paste0("spatial_module_", module_name),
    save_dir = file.path(plot_dir, "spatial_module_genes"),
    width = 24,
    height = 4 * (length(genes) + 1) + 3
  )
}




#----------------------------------------
# Violin by cluster
#----------------------------------------

p_module_scores_violin <- VlnPlot(
  mega_obj,
  features = module_features,
  group.by = "seurat_cluster.projected",
  pt.size = 0
)

save_plot(
  plot_object = p_module_scores_violin,
  file_stub = "module_scores_violin_by_cluster",
  save_dir = plot_dir,
  width = 16,
  height = 12
)








#------------------------------------------------------------
# Heatmap of Relative Module Enrichment by Projected Cluster
#
# Mean module scores are computed for each projected cluster
# and then standardized (z-score) across clusters within each
# module.
#
# Interpretation:
#   Red   = cluster is relatively enriched for the module
#   White = average enrichment
#   Blue  = cluster is relatively depleted for the module
#------------------------------------------------------------
module_features <- paste0(
  setdiff(
    names(module_gene_sets),
    c("Bach1", "Hmox1")
  ),
  "1"
)

cluster_col <- "seurat_cluster.projected"

heatmap_df <- mega_obj[[]] |>
  dplyr::select(
    all_of(cluster_col),
    all_of(module_features)
  ) |>
  dplyr::rename(cluster = all_of(cluster_col)) |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    dplyr::across(
      all_of(module_features),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = -cluster,
    names_to = "Module",
    values_to = "MeanScore"
  ) |>
  dplyr::group_by(Module) |>
  dplyr::mutate(
    ZScore = as.numeric(scale(MeanScore))
  ) |>
  dplyr::ungroup()

#------------------------------------------------------------
# Preserve ordering and use cleaner display names
#------------------------------------------------------------

heatmap_df$Module <- factor(
  heatmap_df$Module,
  levels = module_features,
  labels = c(
    "Urothelium",
    "Fibroblasts",
    "Smooth Muscle",
    "Endothelium",
    "Proliferation",
    "Myeloid",
    "T Cells",
    "NK Cells"
  )
)

heatmap_df$cluster <- factor(
  heatmap_df$cluster,
  levels = sort(unique(as.numeric(as.character(heatmap_df$cluster))))
)

#------------------------------------------------------------
# Plot heatmap
#------------------------------------------------------------

p_module_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = cluster,
    y = Module,
    fill = ZScore
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "Relative\nEnrichment\n(Z-score)"
  ) +
  labs(
    title = "Biological Module Enrichment Across Projected Clusters",
    subtitle = "Mean module scores were standardized (z-score) across projected clusters.",
    x = "Projected Cluster",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    axis.title.x = element_text(
      face = "bold"
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )

save_plot(
  plot_object = p_module_heatmap,
  file_stub = "projected_cluster_module_heatmap_zscore",
  save_dir = plot_dir,
  width = 10,
  height = 6
)














