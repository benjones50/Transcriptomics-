library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ggbeeswarm)


# ------------------------------------------------------------
# Load project configuration
#
#   config.R defines:
#   directories
#   files_list
#   bin_size
# ------------------------------------------------------------
source("/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/code/config.R")



# ------------------------------------------------------------
# Create an empty list to store per-sample QC tables
#
# Each element will be a dataframe for one sample.
# Later we combine them into one big table.
# ------------------------------------------------------------

qc_long_file <- file.path(
  bin_output_dir,
  paste0("qc_long_all_samples_", bin_size, "um.rds")
)


if (file.exists(qc_long_file)) {
  
  message("Loading existing QC table...")
  qc_long <- readRDS(qc_long_file)
  
} else {
  
  message("Building QC table from Seurat objects...")
  
  
qc_list <- list()

# ------------------------------------------------------------
# Loop through all samples
# ------------------------------------------------------------
for (i in seq_along(files_list)) {
  # ----------------------------------------------------------
  # Identify the current sample
  # ----------------------------------------------------------
  
  sample_tissue <- files_list[i]
  sample_name <- paste0(bin_size, "um_", sample_tissue)
  
  
  # ----------------------------------------------------------
  # Load QC-annotated object
  #
  # The loader will:
  #   - load an existing QC object if present
  #   - otherwise build it from the raw object
  #
  # The returned object already contains:
  #   - nCount_
  #   - nFeature_
  #   - percent.mt
  # ----------------------------------------------------------
  
  object <- load_visium_object(
    sample_output_dir = sample_output_dir,
    sample_name = sample_name,
    sample_tissue = sample_tissue,
    raw_data_dir = raw_data_dir,
    bin_size = bin_size,
    version = "qc"
  )
  
  # ----------------------------------------------------------
  # Determine QC column names
  # ----------------------------------------------------------
  
  count_col <- paste0("nCount_",DefaultAssay(object))
  feature_col <- paste0("nFeature_",DefaultAssay(object))
  
  # ----------------------------------------------------------
  # Extract QC metadata
  # ----------------------------------------------------------
  
  qc_df <- object@meta.data %>%
    dplyr::select(
      all_of(
        c(
          count_col,
          feature_col,
          "percent.mt"
        )
      )
    ) %>%
    rename(
      nCount = all_of(count_col),
      nFeature = all_of(feature_col)
    ) %>%
    mutate(
      sample = sample_name,
      sample_tissue = sample_tissue
    )
  
  # ----------------------------------------------------------
  # Store QC table
  # ----------------------------------------------------------
  
  qc_list[[sample_name]] <- qc_df
}






# ------------------------------------------------------------
# Combine all sample QC tables into one dataframe
#
# Each row is one bin/cell/spot.
# ------------------------------------------------------------
qc_all_raw <- bind_rows(qc_list)

# ------------------------------------------------------------
# Keep a filtered version for plotting
#
# Empty bins have:
#   nCount = 0
#   percent.mt = NaN
#
# Those bins are useful to preserve in the raw table, but they
# make QC plots harder to interpret. therefore removed
# ------------------------------------------------------------
qc_all <- qc_all_raw %>%
  filter(nCount > 0)

# ------------------------------------------------------------
# Save both versions
# ------------------------------------------------------------
saveRDS(
  qc_all_raw,
  file.path(bin_output_dir, paste0("qc_all_samples_raw_", bin_size, "um.rds"))
)

saveRDS(
  qc_all,
  file.path(bin_output_dir, paste0("qc_all_samples_filtered_", bin_size, "um.rds"))
)

# ------------------------------------------------------------
# Convert from wide format to long format
#
# Before:
#   sample | nCount | nFeature | percent.mt
#
# After:
#   sample | metric | value
# ------------------------------------------------------------
qc_long <- qc_all %>%
  pivot_longer(
    cols = c(nCount, nFeature, percent.mt),
    names_to = "metric",
    values_to = "value"
  )

# ------------------------------------------------------------
# Save the long-format QC table
# ------------------------------------------------------------
saveRDS(
  qc_long,
  file.path(bin_output_dir, paste0("qc_long_all_samples_", bin_size, "um.rds"))
)

}




# ============================================================
# 1) VIOLIN PLOTS
# ============================================================

# ------------------------------------------------------------
# Violin plot of all QC metrics across samples
#
# Each facet is one metric.
# Each violin is one sample.
# ------------------------------------------------------------
qc_vln_all <- ggplot(
  qc_long,
  aes(x = sample, y = value)
) +
  geom_violin(
    trim = TRUE,
    scale = "width"
  ) +
  facet_wrap(
    ~metric,
    scales = "free_y",
    nrow = 1
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = paste0("QC Metrics Across Samples (", bin_size, " um bins)"),
    x = "Sample",
    y = "Value"
  )


save_plot(qc_vln_all, "QC_vln_all_samples", bin_output_dir, width = 16, height = 6)


# ============================================================
# 2) HISTOGRAMS
# ============================================================

# ------------------------------------------------------------
# Histogram version of the QC distributions
#
# Histograms are very helpful for seeing:
#   - peaks
#   - tails
#   - threshold cutoffs
#   - multimodal distributions
#
# Using alpha helps when many samples overlap.
# ------------------------------------------------------------
qc_hist_all <- ggplot(
  qc_long,
  aes(x = value, fill = sample)
) +
  geom_histogram(
    bins = 80,
    alpha = 0.3,
    position = "identity"
  ) +
  facet_wrap(
    ~metric,
    scales = "free_x",
    nrow = 1
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("QC Histograms Across Samples (", bin_size, " um bins)"),
    x = "Value",
    y = "Count"
  )

save_plot(qc_hist_all, "QC_hist_all_samples",bin_output_dir, width = 16, height = 6)



# ============================================================
# 3) RIDGE PLOTS
# ============================================================
# ------------------------------------------------------------
# Ridge plots show the distribution for each sample stacked
# vertically.
#
# The full distribution is displayed.
# ------------------------------------------------------------

qc_ridge_all <- ggplot(
  qc_long,
  aes(
    x = value,
    y = sample,
    fill = sample
  )
) +
  ggridges::geom_density_ridges(
    alpha = 0.7,
    scale = 1.0,
    rel_min_height = 0.01,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~metric,
    scales = "free_x",
    ncol = 1
  ) +
  theme_bw() +
  labs(
    title = paste0(
      "QC Ridge Plots Across Samples (",
      bin_size,
      " um bins)"
    ),
    x = "Value",
    y = "Sample"
  )

qc_ridge_all

save_plot(qc_ridge_all,"QC_ridge_all_samples",bin_output_dir,width = 10,height = 12)



# ============================================================
# PLOT-ONLY TRIM AT THE 99TH PERCENTILE
# ============================================================

# ------------------------------------------------------------
# Compute a per-metric cutoff at the 99th percentile
# This removes the top 1% of values from the plots only.
# Your saved qc_long table stays unchanged.
# ------------------------------------------------------------
qc_caps <- qc_long %>%
  group_by(metric) %>%
  summarise(
    x_cap = quantile(value, probs = 0.99, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Keep only values at or below the cutoff for plotting
# ------------------------------------------------------------
qc_long_plot <- qc_long %>%
  left_join(qc_caps, by = "metric") %>%
  filter(value <= x_cap) %>%
  select(-x_cap)

# ============================================================
# 1) VIOLIN PLOTS
# ============================================================

qc_vln_all <- ggplot(
  qc_long_plot,
  aes(x = sample, y = value)
) +
  geom_violin(
    trim = TRUE,
    scale = "width"
  ) +
  facet_wrap(
    ~metric,
    scales = "free_y",
    nrow = 1
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = paste0("QC Metrics Across Samples (", bin_size, " um bins; top 1% trimmed)"),
    x = "Sample",
    y = "Value"
  )

#qc_vln_all
save_plot(qc_vln_all, "QC_vln_all_samples_trimmed",  bin_output_dir, width = 16, height = 6)


# ============================================================
# 2) HISTOGRAMS
# ============================================================

qc_hist_all <- ggplot(
  qc_long_plot,
  aes(x = value, fill = sample)
) +
  geom_histogram(
    bins = 80,
    alpha = 0.25,
    position = "identity"
  ) +
  facet_wrap(
    ~metric,
    scales = "free_x",
    nrow = 1
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("QC Histograms Across Samples (", bin_size, " um bins; top 1% trimmed)"),
    x = "Value",
    y = "Count"
  )

#qc_hist_all
save_plot(qc_hist_all, "QC_hist_all_samples_trimmed", bin_output_dir, width = 16, height = 6)


# ============================================================
# 3) RIDGE PLOTS
# ============================================================

qc_ridge_all <- ggplot(
  qc_long_plot,
  aes(
    x = value,
    y = sample,
    fill = sample
  )
) +
  ggridges::geom_density_ridges(
    alpha = 0.7,
    scale = 1.0,
    rel_min_height = 0.01,
    show.legend = FALSE,
    trim = TRUE
  ) +
  facet_wrap(
    ~metric,
    scales = "free_x",
    ncol = 1
  ) +
  theme_bw() +
  labs(
    title = paste0("QC Ridge Plots Across Samples (", bin_size, " um bins; top 1% trimmed)"),
    x = "Value",
    y = "Sample"
  )

#qc_ridge_all
save_plot(qc_ridge_all, "QC_ridge_all_samples_trimmed", bin_output_dir, width = 10, height = 12)







