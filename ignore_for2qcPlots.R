list.files(output_dir)

temp_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/out - june17th_and_before/8um_NL_D_prime/"

list.files(temp_dir)[1]


object <- readRDS(paste0(temp_dir,list.files(temp_dir)[1]))


count_col <- paste0("nCount_", DefaultAssay(object))
feature_col <- paste0("nFeature_", DefaultAssay(object))

VlnPlot(
  object,
  features = c(count_col, feature_col),
  pt.size = 0
)






