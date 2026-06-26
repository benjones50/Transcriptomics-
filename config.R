#config file


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)



#directory where data is being pulled from 

output_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/out"
data_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/data"
code_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code"
raw_data_dir <- "/projects/nagy_lab_projects/projects_benjones/bladder_cancer_spatial_expression/data/spaceranger_count_manual_align_outs/"


#loads helper function to make saving plots ezpz 
source(file.path(code_dir, "plot_saver.R"))
source(file.path(code_dir, "seurat_object_loader.R"))
source(file.path(code_dir,"QC_summary_maker.R"))
source(file.path(code_dir,"cell_segmentation_custom.R"))
source(file.path(code_dir, "QC_metrics_and_filtering.R"))

source(file.path(code_dir,"object_normalizer.R"))
source(file.path(code_dir,"seurat_spatial_fixes.R"))



files_list <- list.files(raw_data_dir) #creates a list of the samples files
files_list 



# ============================================================
# Analysis mode
# ============================================================
# "binned" or "segmented_cells"
analysis_mode <- "segmented_cells"

# ============================================================
# Sample selection
# ============================================================

#must set to NULL when using segmented cells
bin_size <- NULL

sample_tissue_number <- 1

sample_tissue <- files_list[sample_tissue_number]

# ============================================================
# Sample naming
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



#configurables for qc
#should add functionality to print this in graphs #TODO
min_counts <- 15
max_counts <- 1500

min_features <- 15
max_features <- 1250

max_percent_mt <- 20




#changes pt size for spatial graphs, makes bin sizes look better
#maybe in the future should use this for all spatial graphs
#may need adjusting

if (analysis_mode == "segmented_cells") {
  
  pt_size <- 6
  
} else {
  
  pt_size <- switch(
    as.character(bin_size),
    "2"  = 0.5,
    "8"  = 1,
    "16" = 2
  )
  
}
