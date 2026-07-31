#config file
# library(SeuratObject)
# library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(BPCells)
library(future)
library(harmony)
library(ggrepel)

#directory where data is being pulled from 

output_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/out"
data_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/data"
code_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code"
raw_data_dir <- "/projects/nagy_lab_projects/projects_benjones/bladder_cancer_spatial_expression/data/spaceranger_count_manual_align_outs/"

#used for holding polygons to seperate tissues
polygon_file <- "/projects/nagy_lab_projects/projects_benjones/share/data/tissue_separation_polygon_coords.csv"

#loads helper function to make saving plots ezpz 
source(file.path(code_dir, "plot_saver.R"))
source(file.path(code_dir, "seurat_object_loader.R"))
source(file.path(code_dir,"QC_summary_maker.R"))
source(file.path(code_dir,"cell_segmentation_custom.R"))
source(file.path(code_dir, "QC_metrics_and_filtering.R"))

source(file.path(code_dir,"object_normalizer.R"))
source(file.path(code_dir,"seurat_spatial_fixes.R"))

source(file.path(code_dir,"new_tissue_seperator_loader.R"))
source(file.path(code_dir,"mega_umap_plotter.R"))

source(file.path(code_dir,"convert_to_bpcells.R")) #X

source(file.path(code_dir,"file_system.R"))

source(file.path(code_dir,"pca_analysis.R"))
source(file.path(code_dir,"sketch_pca.R"))

source(file.path(code_dir,"clustering.R"))

source(file.path(code_dir,"generic_tissue_renamer.R"))

source(file.path(code_dir,"silhouette.R"))

source(file.path(code_dir,"add_experimental_metadata.R"))
source(file.path(code_dir,"paneled_umap_plotter2.R"))


source(file.path(code_dir,"resolve_normalization_method.R"))

source(file.path(code_dir,"harmonizer.R"))
source(file.path(code_dir,"cluster_tree.R"))
source(file.path(code_dir,"extract_color_palette.R"))


source(file.path(code_dir,"new_diff_expression_clusters_volc.R"))


files_list <- list.files(raw_data_dir) #creates a list of the samples files
files_list 





# # ============================================================
# # Pipeline configuration
# #
# # These are the only options that should normally need changing.
# # The rest of the pipeline automatically adapts based on these.
# # ============================================================

# Analysis mode

# "binned" or "segmented_cells"
analysis_mode <- "binned"

# "lognorm" "sct"
normalization <- "lognorm"
# "pca" or "harmony"
embedding     <- "harmony"       


# Supported:
#   "umap"
#   "tsne"
visualization <- "umap"

# ============================================================
# Sample selection
# ============================================================

#must set to NULL when using segmented cells
bin_size <- 8

sample_tissue_number <- 5 # num of sample in the list

sample_tissue <- files_list[sample_tissue_number]


#configurables for qc
#should add functionality to print this in graphs #TODO
min_counts <- 15
max_counts <- 1500

min_features <- 15
max_features <- 1250

max_percent_mt <- 20

# ============================================================
# Sample naming #used for file output
# ============================================================

if (analysis_mode == "binned") {
  
  sample_name <- paste0(
    bin_size,
    "um_",
    sample_tissue
  )
  
  mode_output_dir <- file.path(
    output_dir,
    paste0(bin_size, "um")
  )
  
} else {
  
  sample_name <- paste0("seg_",sample_tissue)
  
  mode_output_dir <- file.path(
    output_dir,
    "segmented_cells"
  )
  
}

# ============================================================
# Creates analysis-mode folder (2um, 8um, 16um, segmented)
# ============================================================

dir.create(
  mode_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Creates folders on a per sample basis
# ============================================================

sample_output_dir <- file.path(
  mode_output_dir,
  sample_name
)

dir.create(
  sample_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

#directory for grouped outputs
mega_dir <- build_mega_object_dir(
  analysis_mode = analysis_mode,
  normalization = normalization,
  bin_size = bin_size,
  min_counts = min_counts,
  max_counts = max_counts,
  min_features = min_features,
  max_features = max_features,
  max_percent_mt = max_percent_mt
) 







#changes pt size for spatial graphs, makes bin sizes look better
#maybe in the future should use this for all spatial graphs
#may need adjusting

if (analysis_mode == "segmented_cells") {
  
  pt_size <- 6
  
} else {
  
  pt_size <- switch(
    as.character(bin_size),
    "2"  = 0.1,
    "8"  = 0.5,
    "16" = 2
  )
  
}


normalization_method <- switch(
  normalization,
  lognorm = "LogNormalize",
  sct     = "SCT"
)

embedding_reduction <- switch(
  embedding,
  pca     = "pca.sketch",
  harmony = "harmony"
)
projected_reduction <- switch(
  embedding,
  pca     = "projected.pca",
  harmony = "projected.harmony"   # if you decide to name it this
)



visualization_reduction <- switch(
  visualization,
  umap = "umap.sketch",
  tsne = "tsne.sketch"
)

visualization_reduction <- switch(
  visualization,
  umap = "umap.sketch",
  tsne = "tsne.sketch"
)

projected_visualization_reduction <- switch(
  visualization,
  umap = "full.umap.sketch",
  tsne = "full.tsne.sketch"
)



# ============================================================
# Configuration Summary
# ============================================================

# ============================================================
# Configuration Summary
# ============================================================

cat(
  "\n",
  "============================================================\n",
  "           Spatial Analysis Configuration\n",
  "============================================================\n",
  
  "Pipeline\n",
  sprintf("  Analysis mode      : %s\n", analysis_mode),
  sprintf("  Normalization      : %s (%s)\n",
          normalization,
          normalization_method),
  sprintf("  Embedding          : %s (%s)\n",
          embedding,
          embedding_reduction),
  
  "\n",
  
  "Sample\n",
  sprintf("  Sample number      : %d / %d\n",
          sample_tissue_number,
          length(files_list)),
  sprintf("  Sample             : %s\n", sample_tissue),
  sprintf("  Sample name        : %s\n", sample_name),
  sprintf("  Bin size           : %s\n",
          if (is.null(bin_size))
            "Segmented cells"
          else
            paste0(bin_size, " um")),
  
  "\n",
  
  "Quality Control\n",
  sprintf("  Counts             : %d - %d\n",
          min_counts,
          max_counts),
  sprintf("  Features           : %d - %d\n",
          min_features,
          max_features),
  sprintf("  Max %% MT           : %.1f%%\n",
          max_percent_mt),
  
  "\n",
  
  "Directories\n",
  sprintf("  Output             : %s\n", sample_output_dir),
  "\n",
  sprintf("  Mega             : %s\n", mega_dir),
  
  
  "============================================================\n\n",
  sep = ""
)






