#filtering while evaluating markers


#directory where data is being pulled from 

#gets directories
source(file.path(code_dir, "config.R"))



# filtered object must be created in initial_filtering first? maybe should source that file?
# Load filtered Seurat object
filtered_object <- readRDS(
  file.path(
    sample_output_dir,
    paste0("visium_seurat_filtered_", sample_name, ".rds")
  )
)





# Check it loaded correctly
filtered_object


# Load mouse TAM reference database
tam_table <- read.csv(
  paste0(data_dir,"/Mouse_TAM_Reference_Database_LONG_MouseOnly.csv"),
  stringsAsFactors = FALSE
)



# Keep only genes present in the dataset
tam_genes <- unique(tam_table$gene)

#finds TAM genes that were in the table
tam_genes_present <- intersect(
  tam_genes,
  rownames(object)
)


#finds the genes not found in the data at all for that found in visium trial
tam_genes_missing <- setdiff(
  tam_genes,
  rownames(object)
)

cat("Total TAM genes:", length(tam_genes), "\n")
cat("Present:", length(tam_genes_present), "\n")
cat("Missing:", length(tam_genes_missing), "\n")





#begin checking before vs after qc 

# num of tams genes overall detected before qc overall
tam_genes_detected_before <- sum(
  tam_genes_present %in%
    rownames(object)[
      Matrix::rowSums(
        GetAssayData(object, layer = "counts") > 0
      ) > 0
    ]
)

# num of tams genes overall detected after qc overall
tam_genes_detected_after <- sum(
  tam_genes_present %in%
    rownames(filtered_object)[
      Matrix::rowSums(
        GetAssayData(filtered_object, layer = "counts") > 0
      ) > 0
    ]
)

# num tam umis detected before vs after overall
tam_umi_before <- sum(
  GetAssayData(
    object,
    layer = "counts"
  )[tam_genes_present, ]
)

tam_umi_after <- sum(
  GetAssayData(
    filtered_object,
    layer = "counts"
  )[tam_genes_present, ]
)


# percentages of umis and genes overall
pct_tam_genes_retained_overall<- round(
  100 *
    tam_genes_detected_after /
    tam_genes_detected_before,
  2
)
pct_tam_umi_retained_overall <- round(
  100 * tam_umi_after / tam_umi_before,
  2
)


#puts general overview statistics of all tam genes into a table
tam_general_qc_summary <- data.frame(
  metric = c(
    "TAM marker genes",
    "TAM marker UMIs"
  ),
  
  before = c(
    tam_genes_detected_before,
    tam_umi_before
  ),
  
  after = c(
    tam_genes_detected_after,
    tam_umi_after
  ),
  
  pct_retained = c(
    pct_tam_genes_retained_overall,
    pct_tam_umi_retained_overall
  )
)

tam_general_qc_summary

#TODO save this


#specific by gene qc dataframe creation

by_gene_marker_qc_summary <- data.frame()

for(gene in tam_genes_present){
  
  # nums umi before vs after
  umi_before <- sum(
    GetAssayData(
      object,
      layer = "counts"
    )[gene, ]
  )
  
  umi_after <- sum(
    GetAssayData(
      filtered_object,
      layer = "counts"
    )[gene, ]
  )
  # num bins before and after
  bins_before <- sum(
    GetAssayData(
      object,
      layer = "counts"
    )[gene, ] > 0
  )
  
  bins_after <- sum(
    GetAssayData(
      filtered_object,
      layer = "counts"
    )[gene, ] > 0
  )
  
  #puts num bins before / after in a data table
  by_gene_marker_qc_summary <- rbind(
    by_gene_marker_qc_summary,
    data.frame(
      gene = gene,
      
      #umi / genes before vs after, #lost, and %lost
      umi_before = umi_before,
      umi_after = umi_after,
      umi_lost = umi_before - umi_after,
      pct_umi_retained =
        round(
          100 * umi_after / umi_before,
          2
        ),
      
      bins_before = bins_before,
      bins_after = bins_after,
      bins_lost = bins_before - bins_after,
      pct_bins_retained =
        round(
          100 * bins_after / bins_before,
          2
        )
    )
  )
}

#TODO save these
by_gene_marker_qc_summary



#sorts the gene qc summary by biggest losers for easier analysis

#biggest losers by umi #
top10_absolute <- by_gene_marker_qc_summary[
  order(
    by_gene_marker_qc_summary$umi_lost,
    decreasing = TRUE
  ),
]

top10_absolute <- head(top10_absolute, 10)

#biggest losers by umi %
top10_percent <- by_gene_marker_qc_summary[
  order(
    by_gene_marker_qc_summary$pct_umi_retained
  ),
]

top10_percent <- head(top10_percent, 10)

#TODO save these^^^









