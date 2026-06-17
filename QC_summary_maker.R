# ------------------------------------------------------------
# Create or load QC summary statistics
# ------------------------------------------------------------
#
# This function:
#   1) checks whether the saved QC summary already exists
#   2) loads it if it does
#   3) otherwise calculates it from object and filtered_object
#   4) saves it 
#   5) returns the summary table
#


get_qc_summary <- function(
    object,
    filtered_object,
    count_col,
    sample_output_dir,
    sample_name,
    force_rebuild = FALSE
) {
  
  qc_summary_file <- file.path(
    sample_output_dir,
    paste0("qc_summary_table_", sample_name, ".rds")
  )
  
  if (file.exists(qc_summary_file) && !force_rebuild) {
    message("Loading existing QC summary...")
    return(readRDS(qc_summary_file))
  }
  
  #calculates stats comparing before vs after filtering
  before_bins <- ncol(object)
  after_bins  <- ncol(filtered_object)
  
  before_umis <- sum(object[[count_col]][, 1])
  after_umis  <- sum(filtered_object[[count_col]][, 1])
  
  pct_bins_retained <- round(100 * after_bins / before_bins, 2)
  pct_umis_retained <- round(100 * after_umis / before_umis, 2)
  
  
  #combos stats into df
  qc_summary <- data.frame(
    before_bins       = before_bins,
    after_bins        = after_bins,
    removed_bins      = before_bins - after_bins,
    pct_bins_retained = pct_bins_retained,
    pct_bins_removed  = round(100 - pct_bins_retained, 2),
    pct_umis_retained = pct_umis_retained
  )
  #saves rds file
  saveRDS(qc_summary, qc_summary_file)
  qc_summary
}