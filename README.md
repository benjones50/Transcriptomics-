all is WIP, especially this readme

main files:

config.R                     - settings file, bin sizes, directories, etc, loads helper functions


used for making plots: 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plots_using_all_tissues      - for making plots using all 6 tissues
100mt_investigation          - used when I was looking at the bins that had 100% mt counts
QC_analysis_plotting         - used for plotting vlns and spatial plots of each tissue, before filtering, and to show what is removed during filtering. all of individual tissue slides
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


helpful functions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seurat_object_loader.R       - used to load in different seurat objects - uses QC_metrics_and_filtering to generate plots with QC info or filters if not already saved as rds files
  both used by seurat_object_loader:
    QC_metrics_and_filtering.R   - adds the QC metrics (just mt counts, easy to add other stuff), also has a function to apply filters. 
    QC_summary_maker             - creates a df with info about how many bins/umis before vs after filtering, used for some plots


plot_saver.R                 - saving the plots got repetative, this makes it so I dont have to write the code to save to a pdf and png, this function will do both


gene_marker_qc_analysis.R    - not really used yet, will allow us to look at specific umis related to cell type markers that are useful but potentaill being filtered out.


new_tissue_seperator.R       - allows us to more easily differentiate between the 2 tissues on the slides  




