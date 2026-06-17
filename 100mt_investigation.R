# 100% mt counts investigation


#load directories and settings
source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")

#loads seurat object
source(file.path(code_dir, "seurat_object_loader.R"))






#code to evualate the mt only bins

qc_all_raw <- readRDS(
  file.path(
    bin_output_dir,
    paste0("qc_all_samples_raw_", bin_size, "um.rds")
  )
)

cat("\n================ QC SUMMARY ================\n")

cat(
  "\nTotal bins: ",
  format(nrow(qc_all_raw), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "Empty bins (nCount = 0): ",
  format(sum(qc_all_raw$nCount == 0), big.mark = ","),
  " (",
  round(
    100 * mean(qc_all_raw$nCount == 0),
    2
  ),
  "%)\n",
  sep = ""
)

cat(
  "Non-empty bins: ",
  format(sum(qc_all_raw$nCount > 0), big.mark = ","),
  "\n",
  sep = ""
)

mt100 <- qc_all_raw %>%
  filter(percent.mt == 100)

cat(
  "\n100% mitochondrial bins: ",
  format(nrow(mt100), big.mark = ","),
  " (",
  round(
    100 * nrow(mt100) / sum(qc_all_raw$nCount > 0),
    3
  ),
  "% of non-empty bins)\n",
  sep = ""
)

cat("\nDistribution of UMI counts among 100% mitochondrial bins:\n")

print(
  mt100 %>%
    count(nCount, sort = TRUE) %>%
    mutate(
      percent = round(
        100 * n / sum(n),
        2
      )
    )
)


#per tissue evaluation of mc counts
qc_all_raw %>%
  filter(sample == sample_name) %>%
  summarise(
    n_rows = n(),
    n_na = sum(is.na(percent.mt)),
    n_mt100 = sum(percent.mt == 100, na.rm = TRUE)
  )



qc_all_raw %>%
  filter(percent.mt == 100) %>%
  count(sample, sort = TRUE)


object[["percent.mt"]] <- PercentageFeatureSet(
  object,
  pattern = "^mt-"
)

# SpatialFeaturePlot(
#   object,
#   features = "percent.mt"
# )



#makes a plot of just the 100% mt values
mt100_obj <- subset(
  object,
  subset = percent.mt == 100
)

n_mt100 <- ncol(mt100_obj)
n_total <- ncol(object)

mt100_plot <- ImageDimPlot(mt100_obj) +
  labs(
    title = paste0("8um_",sample_tissue, " - 100% Mitochondrial Bins"),
    subtitle = paste0(
      "n = ", format(n_mt100, big.mark = ","),
      " (",
      round(100 * n_mt100 / n_total, 3),
      "% of all bins)"
    )
  )

mt100_plot

#saves mt plots

save_plot(
  mt100_plot, #plot
  paste0("mt100_bins_", sample_name), #name of plot
  save_dir = sample_output_dir, #where saving
  width = 8,
  height = 8,
  dpi = 600
)


