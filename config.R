#config file


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gatepoints)



#directory where data is being pulled from 

output_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/out"
data_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/data"
code_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code"
raw_data_dir <- "/projects/nagy_lab_projects/projects_benjones/bladder_cancer_spatial_expression/data/spaceranger_count_manual_align_outs/"


#loads helper function to make saving plots ezpz 
source(file.path(code_dir, "plot_saver.R"))
source(file.path(code_dir, "seurat_object_loader.R"))
source(file.path(code_dir,"QC_summary_maker.R"))

files_list <- list.files(raw_data_dir) #creates a list of the samples files
files_list 


#selects bin size to analyze
bin_size <- 8

#which sample from the files_list is currently being analyzed (if not all)
sample_tissue_number <- 1



#name of current file being analyzed, will allow for looping through all of the data / bins
sample_name <- paste0(bin_size,"um_",files_list[sample_tissue_number])
sample_tissue <- files_list[sample_tissue_number]



# specific output path for bin sizes
bin_output_dir <- file.path(
  output_dir,
  paste0(bin_size, "um")
)

#creates folder for bin size output
dir.create(bin_output_dir, recursive = TRUE, showWarnings = FALSE)




# Create sample-specific output path
sample_output_dir <- file.path(
  bin_output_dir,
  sample_name
)

#creates folder for sample output
dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)





#loads seurat file
#source((paste0(code_dir,"/seurat_object_loader.R")))



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
pt_size <- switch(
  as.character(bin_size),
  "2"  = .5,
  "8"  = 5,
  "16" = 7,
)

