all is WIP, especially this readme

* = outdated or not currently functional

main files:

config.R                     - settings file, bin sizes, directories, etc, loads helper functions



used for making plots: 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*plots_using_all_tissues      - for making plots using all 6 tissues
*100mt_investigation          - used when I was looking at the bins that had 100% mt counts
*QC_analysis_plotting         - used for plotting vlns and spatial plots of each tissue, before filtering, and to show what is removed during filtering. all of individual

segmented_cell_tester.R       - testing segmented cell outputs, and spatialdimplot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



helpful functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seurat_spatial_fixes.R       - if Read10X_Segmentations() fails, attempts to build object manually, "coords_x_orientation" problem, same for if Load10X_Spatial() fails, uses                                           build_binned_visium_object instead

*seurat_spatial_plot_fixer.R - no coords_x_orientation was treated as evidence that the object was outdated, but our custom objects are ignoring this area, instead If the slot                                        exists, require it to be "horizontal", otherwise assume compatibility and continue. allows compatible objects without the slot to be plotted.




object_normalizer.R          - used to apply log normalization and SCT transform to objects


seurat_object_loader.R       - used to load in different seurat objects - uses QC_metrics_and_filtering to generate plots with info for QC or filters if not already saved as rds files

  used by seurat_object_loader:
    QC_metrics_and_filtering.R   - adds the QC metrics (just mt counts, easy to add other stuff), also has a function to apply filters. 
    QC_summary_maker             - creates a df with info about how many bins/umis before vs after filtering, used for some plots
    cell_segmentation_custom.R   - helper function for seurat_object_loader to load segmented cell data



plot_saver.R                 - saving the plots got repetative, this makes it so I dont have to write the code to save to a pdf and png, this function will do both


gene_marker_qc_analysis.R    - not really used yet, will allow us to look at specific umis related to cell type markers that are useful but potentaill being filtered out.


new_tissue_seperator.R       - allows us to more easily differentiate between the 2 tissues on the slides  




