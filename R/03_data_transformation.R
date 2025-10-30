###################################################################
# Data Transformation for Microbiome of Miscanthus
#
# Transforming to Relative abundance and calculating
# diversity indices. Output is dataframes with relative abundance
# of AVS and observed, Shannon, Simpson's and Inv. Simpson
# diversity indices.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-28
###################################################################

source("R/utils/00_setup.R")

# Relative Abundace Transformation -----------------------------
transform_nested_to_relab <- function(nested_list) {
  purrr::map(nested_list, function(project_list) {
    purrr::map(project_list, function(physeq_obj) {
      transform_sample_counts(physeq_obj, function(x) x / sum(x))
    })
  })
}

main_mxg_relab_psq_list <- transform_nested_to_relab(main_mxg_physeq_list)

# Verify transformations
purrr::map(main_mxg_relab_psq_list, function(project_list) {
  purrr::map(project_list, function(physeq_obj) {
    sample_sums(physeq_obj)
  })
})

# Phyloseq to DF
# main_physeq_relabdf_list <- purrr::map(
#   main_physeq_relab_list,
#   function(project_list) {
#     purrr::map(project_list, function(physeq_obj) {
#       physeq2df(physeq_obj)
#     })
#   }
# )

# Calculate Diversity indices ---------------------------------

main_mxg_relab_df_list <- compute_diversity_nested(
  main_mxg_relab_psq_list,
  drop = "ASV_",
  first_asv_col = "ASV_1"
)

# Save
save(
  main_mxg_physeq_list,
  main_mxg_relab_psq_list,
  main_mxg_relab_df_list,
  file = "data/output/processed/rdata/main_mxg_psq_df_list.rda"
)
