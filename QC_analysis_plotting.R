source(file.path(code_dir, "config.R"))


#this file plots different representations of a single samples qc metrics


#loads / creates objects

# Load QC-annotated object
object <- load_visium_object(
  sample_output_dir = sample_output_dir,
  sample_name = sample_name,
  sample_tissue = sample_tissue,
  raw_data_dir = raw_data_dir,
  bin_size = bin_size,
  version = "qc"
)

# Load filtered object 
filtered_object <- load_visium_object(
  sample_output_dir = sample_output_dir,
  sample_name = sample_name,
  sample_tissue = sample_tissue,
  raw_data_dir = raw_data_dir,
  bin_size = bin_size,
  version = "filtered",
  min_counts = min_counts,
  max_counts = max_counts,
  min_features = min_features,
  max_features = max_features,
  max_percent_mt = max_percent_mt
)





# ----------------------------------------------------------
# Determine QC column names
# ----------------------------------------------------------

count_col <- paste0("nCount_",DefaultAssay(object))
feature_col <- paste0("nFeature_",DefaultAssay(object))



# Visualize QC metrics


#creates spatial feature plots for counts, featues, and mt %
count_spatial <- SpatialFeaturePlot(
  object,
  features = count_col,
  pt.size = pt_size,
  image.alpha = 0.05
) + ggtitle("nCount (Raw UMI Counts)")

feature_spatial <- SpatialFeaturePlot(
  object,
  features = feature_col,
  pt.size = pt_size,
  image.alpha = 0.05
) + ggtitle("nFeatures (# Genes)")

mt_spatial <- SpatialFeaturePlot(
  object,
  features = "percent.mt",
  pt.size = pt_size,
  image.alpha = 0.05
) + ggtitle("% Mitochondrial mRNA")

#row of spatial plots
spatial_row <- count_spatial |
  feature_spatial |
  mt_spatial



# Visualizes the same QC metrics as 3 violin plots
qc_vln <- VlnPlot( 
  object, 
  features = c(count_col, feature_col, "percent.mt"), 
  pt.size = 0, 
  ncol = 3 ) & 
  theme(axis.text = element_text(size = 10)) & 
  NoLegend() 

#will potentially get this warning msg: 
#Warning message: # Removed x rows containing non-finite outside the scale range
#due to some bins being empty (mt % is dividing by 0, NaN results)

#combos vlns with their spatial counterpart
vln_spatial_qc_plot <- qc_vln /
  spatial_row +
  plot_annotation(
    title = paste0(
      sample_name,
      " QC Summary (",
      bin_size,
      " um bins)"
    )
  )

#vln_spatial_qc_plot


#saves combined qc plot
save_plot(
  QC_overview_plot,
  paste0("QC_summary_pdf_before", sample_name),
  sample_output_dir,
  width = 16,
  height = 10,
  dpi = 600
)





# =========================
# QC Summary Statistics
# =========================

#loads a summary of overall before/stats of bins / umi's with filtering for sample

#uses the QC_summary_maker file
qc_summary <- get_qc_summary(
  object = object,
  filtered_object = filtered_object,
  count_col = count_col,
  sample_output_dir = sample_output_dir,
  sample_name = sample_name
)

qc_summary
  


# =========================
# Visualize QC Impact
# =========================

# Label each bin according to whether it passed QC
object$QC <- "Removed"
object$QC[colnames(filtered_object)] <- "Kept"


# Compact QC summary text
qc_text <- paste0(
  "QC Thresholds\n\n",
  "Counts: ", min_counts, " - ", format(max_counts, big.mark = ","), "\n",
  "Features: ", min_features, " - ", format(max_features, big.mark = ","), "\n",
  "MT%: < ", max_percent_mt, "\n\n",
  "Retained bins: ", qc_summary$pct_bins_retained, "%\n",
  "Retained UMIs: ", qc_summary$pct_umis_retained, "%"
)


# Spatial plot of kept vs removed
spatial_plot <- SpatialDimPlot(
  object,
  group.by = "QC",
  cols = c(
    "Removed" = "blue",
    "Kept" = "limegreen"
  ),
  pt.size.factor = pt_size,
  image.alpha = 0
) +
  ggtitle(
    paste0(
      sample_name,
      " | ",
      bin_size,
      " um bins | QC Filtering Results"
    )
  ) &
  theme(
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# QC information panel
qc_panel <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 1,
    label = qc_text,
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void()

# Combine plots into one overview, per sample
QC_overview_plot <- spatial_plot | qc_panel +
  plot_layout(widths = c(6, 1))


#saves QC_overview plot
save_plot(
  QC_overview_plot,
  "QC_overview_plot",
  sample_output_dir,
  width = 10,
  height = 6,
  dpi = 600
)




# Histogram of UMI count distribution, showing red line for minimum

count_df <- data.frame(
  counts = object[[count_col]][,1]
)

hist_plot <- ggplot(
  count_df,
  aes(x = counts)
) +
  geom_histogram(
    bins = 200
  ) +
  geom_vline(
    xintercept = min_counts,
    linetype = "dashed",
    linewidth = 0.75,
    color = "red"
  ) +
  annotate(
    "text",
    x = min_counts,
    y = Inf,
    label = paste0("Min = ", min_counts),
    angle = 90,
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = paste0(
      sample_name,
      " UMI Count Distribution"
    ),
    subtitle = paste0(
      "Minimum count threshold = ",
      min_counts
    ),
    x = "UMI Counts per Bin",
    y = "Number of Bins"
  ) +
  theme_bw(base_size = 14)


#saves plot
save_plot(
  hist_plot,
  paste0("UMI_count_histogram_", sample_name),
  sample_output_dir,
  width = 10,
  height = 6,
  dpi = 600
)



#feature scatter plots of qc

plot1 <- FeatureScatter(
  object,
  feature1 = count_col,
  feature2 = "percent.mt"
) +
  coord_cartesian(ylim = c(0, 40)) +
  ggtitle(
    paste0(sample_name, " - Counts vs % Mitochondrial")
  )

plot2 <- FeatureScatter(
  object,
  feature1 = count_col,
  feature2 = feature_col
) +
  ggtitle(
    paste0(sample_name, " - Counts vs Features")
  )

scatter_plot <- plot1 + plot2


#saving scatters 

#uses save plot function
save_plot(
  scatter_plot,
  paste0("QC_scatter_plot_", sample_name),
  sample_output_dir,
  width = 10,
  height = 5
)

                          