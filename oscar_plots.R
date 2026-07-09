##
# Quality control of Visium HD data (pancreas diabetes type I and healthy samples)
# Dr. Matthew Poy's project
#
# By Oscar Ospina
# Created: Mar 26, 2026
# Modified: Mar 26, 2026
#

# Large data set handling
options(future.globals.maxSize=1e10)

library('Seurat')
library('BPCells')

dir.create('./results/quality_ctrl/', recursive=TRUE, showWarnings=FALSE)

# Read Seurat object
seurat_obj = readRDS(file='./data/pancreas_diabetes_tpI_visiumhd.RDS')

dim(seurat_obj)
############ Data dimensions: 18,132 genes x 335,111 pixels

# Calculate other MT-DNA percentage stats
seurat_obj[['percent_100_mt']] = PercentageFeatureSet(seurat_obj, pattern='^MT\\-')

# Plot distribution of counts and genes per pixel
vp1 = VlnPlot(seurat_obj, feature='nCount_RNA', group.by='orig.ident', alpha=0.1, pt.size=0.1) +
  labs(title='RNA counts per cell', x='', y='Counts')
vp2 = VlnPlot(seurat_obj, feature='nFeature_RNA', group.by='orig.ident', alpha=0.1, pt.size=0.1) +  
  labs(title='Detected genes per cell', x='', y='Genes')
vp3 = VlnPlot(seurat_obj, feature='percent_100_mt', group.by='orig.ident', alpha=0.1, pt.size=0.1) +  
  labs(title='Percent of MT RNA counts per cell', x='', y='Counts')

graphics.off()
pdf('./results/quality_ctrl/count_and_gene_per_pixel_distributions.pdf', width=14, height=6)
print(ggpubr::ggarrange(vp1, vp2, vp3, ncol=3))
dev.off()

rm(vp1, vp2, vp3) # Clean env

# Plot spatial distributions of counts and detected genes
sp1 = ggplot(seurat_obj@meta.data) + 
  geom_point(aes(y=pxl_row_in_fullres, x=pxl_col_in_fullres, color=nCount_RNA), size=0.1) + 
  khroma::scale_color_YlOrBr(midpoint=100) +
  scale_y_reverse() + coord_equal() + 
  theme_void() +
  facet_wrap(~orig.ident)

sp2 = ggplot(seurat_obj@meta.data) + 
  geom_point(aes(y=pxl_row_in_fullres, x=pxl_col_in_fullres, color=nFeature_RNA), size=0.1) + 
  khroma::scale_color_YlOrBr(midpoint=100) +
  scale_y_reverse() + coord_equal() + 
  theme_void() +
  facet_wrap(~orig.ident)

sp3 = ggplot(seurat_obj@meta.data) + 
  geom_point(aes(y=pxl_row_in_fullres, x=pxl_col_in_fullres, color=percent_100_mt), size=0.1) + 
  khroma::scale_color_YlOrBr(midpoint=20) +
  scale_y_reverse() + coord_equal() + 
  theme_void() +
  facet_wrap(~orig.ident)

graphics.off()
png('./results/quality_ctrl/count_and_gene_per_pixel_spatial.png', width=30, height=6, res=300, units='in')
print(ggpubr::ggarrange(sp1, sp2, sp3, ncol=3))
dev.off()

rm(sp1, sp2, sp3) # Clean env

# Filter using hard thresholds
seurat_obj = subset(seurat_obj, nCount_RNA >= 200 &
                      nFeature_RNA >= 200 #&
                    #percent_100_mt <= 20 # Keep for now, since very low MT-DNA counts overall
)

# Plot distribution of counts and genes per pixel
vp1 = VlnPlot(seurat_obj, feature='nCount_RNA', group.by='orig.ident', alpha=0.1, pt.size=0.1, ) +
  labs(title='RNA counts per cell', x='', y='Counts')
vp2 = VlnPlot(seurat_obj, feature='nFeature_RNA', group.by='orig.ident', alpha=0.1, pt.size=0.1) +  
  labs(title='Detected genes per cell', x='', y='Genes')
vp3 = VlnPlot(seurat_obj, feature='percent_100_mt', group.by='orig.ident', alpha=0.1, pt.size=0.1) +  
  labs(title='Percent of MT RNA counts per cell', x='', y='Counts')

graphics.off()
pdf('./results/quality_ctrl/count_and_gene_per_pixel_distributions_post_hard_filter.pdf', width=14, height=6)
print(ggpubr::ggarrange(vp1, vp2, vp3, ncol=3))
dev.off()

# Plot spatial distributions of counts and detected genes
sp1 = ggplot(seurat_obj@meta.data) + 
  geom_point(aes(y=pxl_row_in_fullres, x=pxl_col_in_fullres, color=nCount_RNA), size=0.1) + 
  khroma::scale_color_YlOrBr(midpoint=100) +
  scale_y_reverse() + coord_equal() + 
  theme_void() +
  facet_wrap(~orig.ident)

sp2 = ggplot(seurat_obj@meta.data) + 
  geom_point(aes(y=pxl_row_in_fullres, x=pxl_col_in_fullres, color=nFeature_RNA), size=0.1) + 
  khroma::scale_color_YlOrBr(midpoint=100) +
  scale_y_reverse() + coord_equal() + 
  theme_void() +
  facet_wrap(~orig.ident)

sp3 = ggplot(seurat_obj@meta.data) + 
  geom_point(aes(y=pxl_row_in_fullres, x=pxl_col_in_fullres, color=percent_100_mt), size=0.1) + 
  khroma::scale_color_YlOrBr(midpoint=20) +
  scale_y_reverse() + coord_equal() + 
  theme_void() +
  facet_wrap(~orig.ident)

graphics.off()
png('./results/quality_ctrl/count_and_gene_per_pixel_spatial_post_hard_filter.png', width=30, height=6, res=300, units='in')
print(ggpubr::ggarrange(sp1, sp2, sp3, ncol=3))
dev.off()

dim(seurat_obj)
############ Data dimensions: 18,132 genes x 264,482 pixels

# Remove mitochondrial genes
seurat_obj = seurat_obj[grep('^MT\\-', rownames(seurat_obj), value=TRUE, invert=TRUE), ]

# Remove zero-count genes
genes_keep = names(which(Matrix::rowSums(GetAssayData(JoinLayers(seurat_obj), layer='counts')) > 0))
seurat_obj = seurat_obj[genes_keep, ]

dim(seurat_obj)
###### Data set at this point: 18,100 genes x 264,482 cells ###### 

# Normalize with log-transformation
seurat_obj = NormalizeData(seurat_obj, 
                           scale.factor=median(seurat_obj$nCount_RNA),
                           assay='RNA')

# Save normalized data set
saveRDS(seurat_obj, './data/pancreas_diabetes_tpI_visiumhd_lognorm.RDS')

