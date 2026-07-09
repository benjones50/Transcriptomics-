##
# Read Visium HD spatial raw data (segmented cell outputs)
# Bladder cancer in conditional KO mice (Dr. Michalis Sarris)
#
# By Oscar Ospina
# Created: Mar 26, 2026
# Modified: Jun 16, 2026
#

# Large data set handling
options(future.globals.maxSize=1e10)

library('Seurat')
library('BPCells')
library('tidyverse')

# Get sample names
snames = list.dirs('./data/spaceranger_count_manual_align_outs/', recursive=FALSE)
snames = stringr::str_extract(snames, '[_0-9A-Za-z]+$')

# List Visium HD Space Ranger output directories
vhd_fps = sapply(snames, function(i){
  fp_tmp = list.files(paste0('./data/spaceranger_count_manual_align_outs/', i, '/outs/'), 
                      recursive=TRUE, pattern='filtered_feature_cell_matrix.h5', full.names=TRUE)
  return(fp_tmp)
})

# # Load data and create on-disk matrices
# lapply(snames, function(i){
#   fp_tmp = vhd_fps[[i]]
#   vhd_h5_tmp = open_matrix_10x_hdf5(path=fp_tmp)
#   outfp_tmp = paste0('./data/ondisk_visium_hd_mtx_segmented/', i)
#   write_matrix_dir(mat=vhd_h5_tmp, dir=outfp_tmp, overwrite=TRUE)
# })

# Read matrices
vhd_ls = lapply(snames, function(i){
  mtx_tmp = paste0('./data/ondisk_visium_hd_mtx_segmented/', i)
  vhd_mtx = open_matrix_dir(dir=mtx_tmp)
  return(vhd_mtx)
})
names(vhd_ls) = snames

# Create Seurat objects
seurat_ls = lapply(snames, function(i){
  mat_tmp = vhd_ls[[i]]
  # Read gene symbols
  fp_tmp = paste0('./data/spaceranger_count_manual_align_outs/', i, '/outs/segmented_outputs/filtered_feature_cell_matrix/features.tsv.gz')
  feat_tmp = readr::read_delim(fp_tmp, delim="\t", col_names=FALSE) %>% 
    dplyr::filter(X3 == "Gene Expression") %>%
    dplyr::select(ensembl=X1, symbol=X2)
  # Match the matrix rownames (Ensembl IDs) to the features table
  idx = match(rownames(mat_tmp), feat_tmp$ensembl)
  gene_symbol = feat_tmp$symbol[idx]
  # If symbol is missing keep the Ensembl ID
  gene_symbol[is.na(gene_symbol) | gene_symbol == ""] = rownames(mat_tmp)[is.na(gene_symbol) | gene_symbol == ""]
  # Collapse duplicated gene symbols by summing counts
  ## Build a sparse mapping matrix: gene_symbol x ensembl_id
  f = factor(gene_symbol, levels=unique(gene_symbol))
  map_mat = Matrix::sparseMatrix(i=as.integer(f), j=seq_along(f), x=1,
                                 dims=c(nlevels(f), length(f)),
                                 dimnames=list(levels(f), rownames(mat_tmp)))
  ## Sum rows with the same gene symbol
  mat_gene = map_mat %*% mat_tmp
  seu_tmp = CreateSeuratObject(counts=mat_gene, assay="Spatial", project=i)
  # Add cell polygons
  out_fp = paste0('./data/spaceranger_count_manual_align_outs/', i, '/outs')
  sp_fp = paste0(out_fp, '/segmented_outputs/spatial/')
  pol_obj = Read10X_Segmentations(image.dir=sp_fp, 
                                  data.dir=out_fp, 
                                  slice=paste0("slice1.", i),
                                  assay='Spatial')
  # Remove zero-count cells
  common_tmp = intersect(Cells(pol_obj), colnames(seu_tmp))
  pol_obj = subset(pol_obj, cells=common_tmp)
  # Make sure cells are in the same order... Otherwise throw error
  if(!identical(Cells(pol_obj), colnames(seu_tmp))){
    stop("Cell IDs in Seurat object and polygon object are not the same or not in the same order")
  } else{
    seu_tmp[[paste0("slice1.", i)]] = pol_obj 
  }
  return(seu_tmp)
})
names(seurat_ls) = snames

# Save Seurat object list
saveRDS(seurat_ls, file='./data/bladder_cancer_seurat_segmented_outs.RDS')

lapply(seurat_ls, dim)
############ Data dimensions:
# $NL_A_prime
# [1] 19065 57689
# 
# $NL_A1_prime
# [1] 19065 39236
# 
# $NL_B
# [1] 19065 67895
# 
# $NL_C
# [1] 19065 47609
# 
# $NL_D_prime
# [1]  19065 137198
# 
# $NL_E
# [1] 19065 58988