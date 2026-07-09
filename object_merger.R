source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")

# #tells seurat to use v5 class by default, idk if necessary
# options(Seurat.object.assay.version = "v5") 


# Configure the future framework to allow Seurat functions (e.g. SCTransform)
# to use up to 6 parallel R worker processes when they support parallel execution.
#plan(multisession, workers = 6)
#NormalizeData, ScaleData, JackStraw, FindMarkers, FindIntegrationAnchors, FindClusters 
#plan(sequential)

options(future.globals.maxSize = 1000 * 1024^2)

force_rebuild <- FALSE


# Put your 6 slide/sample folder names here
slide_ids <- files_list[1:6]

load_one_slide <- function(sample_tissue) {
  #current sample
  print(sample_tissue)
  
  obj <- load_normalized_object(
    sample_tissue   = sample_tissue,
    raw_data_dir    = raw_data_dir,
    analysis_mode   = analysis_mode,
    bin_size        = bin_size,
    normalization   = normalization,
    min_counts      = min_counts,
    max_counts      = max_counts,
    min_features    = min_features,
    max_features    = max_features,
    max_percent_mt  = max_percent_mt,
    force_rebuild   = force_rebuild
  )
  
  
  # Always store easy-to-merge metadata
  obj$slide_id <- sample_tissue
  obj$analysis_mode <- analysis_mode
  obj$bin_size <- bin_size
  
  
  # adds metadata columns specifying tissue origin
  tissue_seperated_results <- assign_tissues_by_polygon( #from new_tissue_seperator_loader
    obj,
    sample_tissue,
    polygon_file, #from config
    output_dir = mode_output_dir,
  )
  
  obj <- tissue_seperated_results$object
  
  
  
  
  
  
  
  #builds a path to the objects rds file, used to store bp cells file
  # object_file <- build_object_path(
  #   sample_tissue = sample_tissue,
  #   stage = normalization,
  #   analysis_mode = analysis_mode,
  #   bin_size = bin_size,
  #   min_counts = min_counts,
  #   max_counts = max_counts,
  #   min_features = min_features,
  #   max_features = max_features,
  #   max_percent_mt = max_percent_mt
  # )
  # 
  
  #makes individual object use bpcells to store on file
  # obj <- convert_to_bpcells(
  #   object = obj,
  #   object_file = object_file
  # )
  
  return(obj)
  
}



obj_list <- setNames(lapply(slide_ids, load_one_slide), slide_ids)


# class(obj_list)
# 
# class(obj_list[[1]])
# 
# str(obj_list[[1]], max.level = 1)

#debug individual slides, irgnore
#obj_list <- setNames(lapply(slide_ids[5], load_one_slide), slide_ids[5])

# slide_ids[6]
# 
# obj_list <- setNames(
#   lapply(slide_ids[6], load_one_slide),
#   slide_ids[6]
# )



# Merge into one master object

mega_obj <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list),
  project = "6_objects"
)

rm(obj_list)

#some info to make sure it worked
colnames(mega_obj[[]])

table(mega_obj$slide_id)

validObject(mega_obj)

Reductions(mega_obj)

Images(mega_obj)

Layers(mega_obj[[DefaultAssay(mega_obj)]])

class(mega_obj)




# At this point you have one Seurat object with all tissues tagged in metadata.
# If you want a unified assay after merging and are NOT doing integration:

assay <- DefaultAssay(mega_obj)

if (inherits(mega_obj[[assay]], "Assay5")) { #checks if sct vs lognorm basically
  mega_obj[[assay]] <- JoinLayers(mega_obj[[assay]])
}
Layers(mega_obj[[assay]])





# If you want batch correction / co-embedding across slides:
# mega_obj[["RNA"]] <- split(mega_obj[["RNA"]], f = mega_obj$slide_id)
# mega_obj <- FindVariableFeatures(mega_obj)
# mega_obj <- RunPCA(mega_obj)
# mega_obj <- IntegrateLayers(
#   object = mega_obj,
#   method = RPCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "integrated.rpca",
#   verbose = FALSE
# )
# mega_obj[["RNA"]] <- JoinLayers(mega_obj[["RNA"]])

# Save it


save_mega_object(
  object = mega_obj,
  directory = mega_dir
)





