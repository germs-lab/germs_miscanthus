################################################################################
# Rebuild Phyloseq object and Transformation for Microbiome of Miscanthus
#
# This script focuses on reassigning to the 16S_DNA phyloseq objects the newly generated OTU (02_main_subset_to_fasta.R) and taxonomy tables (04_redo_dada2.R). The new objects contain the concatenated OTU table of all projects and it's corresponding taxonomy table. Objects are then subset to the specific samples for each project, thus retaining only the OTUs present in each. OTU1 in project EF should be the same as OTU1 in LAMPS 2018. They were assigned and renamed all together.

# We also perform Relative abundance and Hellinger transformations to "main_physeq_list" and "main_mxg_physe_list", the former contains all project data and metadata, the latter contains all data and metadata focused on MXG crop.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-28
################################################################################

source("R/utils/00_setup.R")

# Rebuild phyloseq objects ----
## Clean up to taxonomy table ----
new_taxa_16S_DNA <- read.csv("data/output/taxonomy/new_16S_DNA_tax.csv") |>
  rename(sequence = X) |>
  janitor::clean_names()

rownames(new_taxa_16S_DNA) <- new_taxa_16S_DNA$sequence

new_taxa_16S_DNA <- tax_table(as.matrix(new_taxa_16S_DNA))

## New main_16S_DNA_physeq_list ----
main_16S_physeq_list <- purrr::map(
  main_16S_physeq_list,
  function(project_names) {
    tax_table(project_names) <- tax_table(new_taxa_16S_DNA)
    project_names <- add_refseq(project_names, tag = NA, seq_to_name = TRUE)

    return(project_names)
  }
)

## Inspect ----
purrr::map(main_16S_physeq_list, function(project_names) {
  dim(tax_table(project_names))
})
# All should have 8 columns (now it includes species column)

# Relative Abundance Transformation ----
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

## Relative abundance ----
main_mxg_relab_psq_list <- transform_to_relab(main_mxg_physeq_list)

## Hellinger transformation ----
main_mxg_hellgr_psq_list <- transform_to_hellinger(main_mxg_physeq_list)
main_hellgr_physeq_list <- transform_to_hellinger(main_physeq_list)

## Verify transformations ----
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

# Save objects ----
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

save(main_16S_physeq_list, file = "data/output/rdata/main_16S_physeq_list.rda")
