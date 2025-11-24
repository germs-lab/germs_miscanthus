###################################################################
# Data Transformation for Microbiome of Miscanthus
#
# Transforming to Relative abundance.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-28
###################################################################

source("R/utils/00_setup.R")

# Relative Abundace Transformation -----------------------------
transform_to_relab <- function(nested_list) {
  purrr::map(nested_list, function(project_list) {
    purrr::map(project_list, function(physeq_obj) {
      transform_sample_counts(physeq_obj, function(x) x / sum(x))
    })
  })
}

transform_to_hellinger <- function(nested_list) {
  purrr::map(nested_list, function(project_list) {
    purrr::map(project_list, function(physeq_obj) {
      otu_mat <- otu_table(physeq_obj) %>%
        as.matrix()

      # Apply Hellinger transformation
      # decostand expects rows as samples, so transpose if needed
      if (taxa_are_rows(physeq_obj)) {
        otu_mat <- t(otu_mat)
      }

      # Apply Hellinger transformation
      otu_hellinger <- vegan::decostand(
        otu_mat,
        method = "hellinger"
      ) %>%
        t() # Transpose back to taxa_are_rows = TRUE

      otu_table(physeq_obj) <- otu_table(otu_hellinger, taxa_are_rows = TRUE)

      return(physeq_obj)
    })
  })
}

# Ratalive abundance
main_mxg_relab_psq_list <- transform_to_relab(main_mxg_physeq_list)

# Hellinger transformation
main_mxg_hellgr_psq_list <- transform_to_hellinger(main_mxg_physeq_list)
main_hellgr_physeq_list <- transform_to_hellinger(main_physeq_list)


# Verify transformations
purrr::map(main_mxg_relab_psq_list, function(project_list) {
  purrr::map(project_list, function(physeq_obj) {
    sample_sums(physeq_obj)
  })
})

main_hellgr_physeq_list$mxg_ef$ef_16S_physeq %>%
  otu_table() %>%
  as.matrix() %>%
  head() # Should see decimal values between 0 and 1


# Phyloseq to DF
# main_physeq_relabdf_list <- purrr::map(
#   main_physeq_relab_list,
#   function(project_list) {
#     purrr::map(project_list, function(physeq_obj) {
#       physeq2df(physeq_obj)
#     })
#   }
# )

# Save
save(
  main_mxg_physeq_list,
  main_mxg_relab_psq_list,
  main_mxg_hellgr_psq_list,
  file = "data/output/rdata/main_mxg_psq_list.rda"
)

save(
  main_hellgr_physeq_list,
  file = "data/output/rdata/main_hellgr_physeq_list.rda"
)
