# ============================================================
# SEGMENTED CELL OBJECT TESTING SCRIPT
# ============================================================

source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")

# ============================================================
# Load objects
# ============================================================

object <- load_visium_object(
  sample_output_dir = sample_output_dir,
  sample_name = sample_name,
  analysis_mode = analysis_mode,
  sample_tissue = sample_tissue,
  raw_data_dir = raw_data_dir,
  bin_size = bin_size,
  version = "qc",
  force_rebuild = FALSE
)




# ============================================================
# Basic object information
# ============================================================

cat("\n==============================\n")
cat("OBJECT SUMMARY\n")
cat("==============================\n")

print(object)

cat("\nAssays:\n")
print(Assays(object))

cat("\nDefault assay:\n")
print(DefaultAssay(object))

cat("\nImages:\n")
print(Images(object))

cat("\nDimensions:\n")
print(dim(object))

cat("\nMetadata columns:\n")
print(colnames(object@meta.data))

# ============================================================
# Verify image attachment
# ============================================================

cat("\n==============================\n")
cat("IMAGE CHECK\n")
cat("==============================\n")

img <- object@images[[1]]

cat("\nImage class:\n")
print(class(img))

cat("\nBoundary names:\n")
print(names(img@boundaries))

cat("\nCells in image:\n")
print(length(Cells(img)))

cat("\nCells in assay:\n")
print(ncol(object))

cat("\nIdentical ordering:\n")
print(
  identical(
    Cells(img),
    colnames(object)
  )
)

# ============================================================
# QC metric checks
# ============================================================

cat("\n==============================\n")
cat("QC CHECK\n")
cat("==============================\n")

count_col <- paste0(
  "nCount_",
  DefaultAssay(object)
)

feature_col <- paste0(
  "nFeature_",
  DefaultAssay(object)
)

summary(
  object[[count_col]][,1]
)

summary(
  object[[feature_col]][,1]
)

if ("percent.mt" %in% colnames(object@meta.data)) {
  
  cat("\npercent.mt summary:\n")
  
  print(
    summary(
      object$percent.mt
    )
  )
}

# ============================================================
# Look at gene names
# ============================================================

cat("\n==============================\n")
cat("GENE CHECK\n")
cat("==============================\n")

cat("\nFirst 20 genes:\n")

print(
  head(
    rownames(object),
    20
  )
)

cat("\nCommon epithelial markers:\n") #from gpt, just for testing. 

print(
  intersect(
    c(
      "Epcam",
      "Krt8",
      "Krt18",
      "Krt19",
      "Krt20"
    ),
    rownames(object)
  )
)

# ============================================================
# Violin QC plots
# ============================================================

qc_vln <- VlnPlot(
  object,
  features = c(
    count_col,
    feature_col,
    "percent.mt"
  ),
  pt.size = 0
)
qc_vln


# ============================================================
# Segmentation-only plot, outline of tissue
# ============================================================

seg_plot <- ImageDimPlot(
  object,
  border.color = NA,
  axes = FALSE
)

seg_plot

# ============================================================
# Plot a marker gene that actually exists
# ============================================================

candidate_genes <- c(
  "Epcam",
  "Krt8",
  "Krt18",
  "Krt19",
  "Krt20",
  "Pecam1",
  "Ptprc",
  "Col1a1"
)

gene_to_plot <- candidate_genes[
  candidate_genes %in% rownames(object)
][3]

cat("\nGene selected for plotting:\n")
print(gene_to_plot)

if (!is.na(gene_to_plot)) {
  
  gene_plot <- ImageFeaturePlot(
    object,
    features = gene_to_plot,
    fov = Images(object)[1],
    boundaries = "segmentations",
    border.color = NA,
    crop = TRUE
  )
  
  print(gene_plot)
}

# ============================================================
# Save diagnostic plots
# ============================================================

save_plot(
  qc_vln,
  "test_segmented_QC_violin",
  sample_output_dir,
  width = 12,
  height = 4
)

save_plot(
  seg_plot,
  "test_segmented_cell_boundaries",
  sample_output_dir,
  width = 8,
  height = 8
)


save_plot(
  gene_plot,
  paste0("test_gene_", gene_to_plot),
  sample_output_dir,
  width = 8,
  height = 8
)

qc_vln
seg_plot
gene_plot



#qc plots for segmented cells

# ============================================================
# Determine QC column names
# ============================================================

count_col <- paste0(
  "nCount_",
  DefaultAssay(object)
)

feature_col <- paste0(
  "nFeature_",
  DefaultAssay(object)
)

# ============================================================
# Spatial QC plots
# ============================================================

count_spatial <- ImageFeaturePlot(
  object,
  features = count_col,
  fov = Images(object)[1],
  boundaries = "segmentations",
  border.color = NA,
  crop = TRUE
) +
  ggtitle("nCount (Raw UMI Counts)")

feature_spatial <- ImageFeaturePlot(
  object,
  features = feature_col,
  fov = Images(object)[1],
  boundaries = "segmentations",
  border.color = NA,
  crop = TRUE
) +
  ggtitle("nFeatures (# Genes)")

mt_spatial <- ImageFeaturePlot(
  object,
  features = "percent.mt",
  fov = Images(object)[1],
  boundaries = "segmentations",
  border.color = NA,
  crop = TRUE
) +
  ggtitle("% Mitochondrial mRNA")

# ============================================================
# Combine spatial plots
# ============================================================

spatial_row <-
  count_spatial |
  feature_spatial |
  mt_spatial

# ============================================================
# Violin plots
# ============================================================

qc_vln <- VlnPlot(
  object,
  features = c(
    count_col,
    feature_col,
    "percent.mt"
  ),
  pt.size = 0,
  ncol = 3
) &
  theme(
    axis.text = element_text(size = 10)
  ) &
  NoLegend()

# ============================================================
# Combined QC figure
# ============================================================

QC_overview_plot <-
  qc_vln /
  spatial_row +
  plot_annotation(
    title = paste0(
      sample_name,
      "_QC_Summary"
    )
  )

save_plot(
  QC_overview_plot,
  "test_seg_QC_overview_plot",
  sample_output_dir,
  width = 8,
  height = 8
)

QC_overview_plot
