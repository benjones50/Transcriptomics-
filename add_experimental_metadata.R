
# ------------------------------------------------------------
# Add experimental metadata
#
# Creates standardized experimental metadata columns from the
# original tissue names. This includes:
#   - experimental_group
#   - tissue_id
#   - experimental_group_tissue
#   - genotype
#   - treatment
#   - timepoint
#   - phenotype
#
# Also converts the relevant metadata columns into ordered
# factors for consistent plotting.
#
# Returns the modified Seurat object.
# ------------------------------------------------------------

add_experimental_metadata <- function(
    object,
    force_rebuild = FALSE
) {
  
  
  #checking if data already exists
  required_columns <- c(
    "experimental_group",
    "tissue_id",
    "experimental_group_tissue",
    "genotype",
    "treatment",
    "timepoint",
    "phenotype"
  )
  
  if (
    !force_rebuild &&
    all(required_columns %in% colnames(object[[]]))
  ) {
    message("Experimental metadata already exists. Skipping.")
    return(object)
  }
  
  
  
  
  
  
  lookup <- c(
    "NL_A_prime_left"   = "WT_Control__Untreated__12_20wks__Normal_Homeostasis",
    "NL_A_prime_right"  = "cKO_BACH1__Untreated__12_20wks__Normal_Homeostasis",
    
    "NL_A1_prime_left"  = "cKO_BACH1__BBN__2wks__Dysplasia_Inflammation",
    "NL_A1_prime_right" = "BACH1_Control__BBN__2wks__Dysplasia_Inflammation",
    
    "NL_B_left"         = "BACH1_Control__BBN__15wks__Cis_Early_Invasion",
    "NL_B_right"        = "cKO_BACH1__BBN__15wks__Cis_Early_Invasion",
    
    "NL_C_left"         = "BACH1_Control__BBN__15wks__Cis_Early_Invasion",
    "NL_C_right"        = "cKO_BACH1__BBN__15wks__Cis_Early_Invasion",
    
    "NL_D_prime_left"   = "Nrf2_Control__BBN__19wks__Invasion",
    "NL_D_prime_right"  = "BACH1_Control__BBN__19wks__Invasion",
    
    "NL_E_left"         = "cKO_BACH1__BBN__19wks__Invasion",
    "NL_E_right"        = "cKO_BACH1__BBN__19wks__Invasion"
  )
  
  # ----------------------------------------------------------
  # Create experimental group from tissue names
  # ----------------------------------------------------------
  
  object <- rename_metadata_values( #generic_tissue_renamer
    object = object,
    column = "tissue",
    lookup = lookup,
    new_column = "experimental_group"
  )
  
  # ----------------------------------------------------------
  # Create compact tissue identifier
  # ----------------------------------------------------------
  
  object$tissue_id <- object$tissue
  object$tissue_id <- gsub("^NL_", "", object$tissue_id)
  object$tissue_id <- gsub("_prime", "", object$tissue_id)
  object$tissue_id <- gsub("_left$", "_L", object$tissue_id)
  object$tissue_id <- gsub("_right$", "_R", object$tissue_id)
  object$tissue_id <- gsub("_", ".", object$tissue_id)
  
  # ----------------------------------------------------------
  # Create unique sample identifier
  # ----------------------------------------------------------
  
  object$experimental_group_tissue <- paste(
    object$experimental_group,
    gsub(" ", "_", object$tissue_id),
    sep = "__"
  )
  
  # ----------------------------------------------------------
  # Split experimental group into individual metadata columns
  # ----------------------------------------------------------
  
  object <- split_metadata_column(
    object = object,
    column = "experimental_group",
    into = c(
      "genotype",
      "treatment",
      "timepoint",
      "phenotype"
    )
  )
  
  # ----------------------------------------------------------
  # Order treatments
  # ----------------------------------------------------------
  
  object$treatment <- factor(
    object$treatment,
    levels = c(
      "Untreated",
      "BBN"
    )
  )
  # ----------------------------------------------------------
  # Order timepoints chronologically
  # ----------------------------------------------------------
  
  object$timepoint <- factor(
    object$timepoint,
    levels = c(
      "2wks",
      "12_20wks",
      "15wks",
      "19wks"
    )
  )
  
  # ----------------------------------------------------------
  # Order experimental groups by timepoint
  # ----------------------------------------------------------
  
  experimental_group_levels <- unique(
    object@meta.data[
      order(
        object$treatment,
        object$timepoint,
        object$experimental_group
      ),
      "experimental_group"
    ]
  )
  
  object$experimental_group <- factor(
    object$experimental_group,
    levels = experimental_group_levels
  )
  
  # ----------------------------------------------------------
  # Order experimental_group_tissue by timepoint then tissue
  # ----------------------------------------------------------
  
  group_tissue_levels <- unique(
    object@meta.data[
      order(
        object$treatment,
        object$timepoint,
        object$experimental_group,
        object$tissue_id
      ),
      "experimental_group_tissue"
    ]
  )
  
  object$experimental_group_tissue <- factor(
    object$experimental_group_tissue,
    levels = group_tissue_levels
  )
  

  return(object)
}