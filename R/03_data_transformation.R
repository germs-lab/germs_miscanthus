###################################################################
# Data Transformation for Microbiome of Miscanthus
#
# Transforming to Relative abundance and calculating
# diversity indices. Output is dataframes with relative abundance
# of AVS and observed, Shannon, Simpson's and Inv. Simpson
# diversity indices.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-28#
###################################################################

source("R/utils/00_setup.R")


# Summary statistics for datasets
main_physeq_list <- list(
  ef_physeq_list = ef_physeq_list,
  lamps_2018_physeq_list = lamps_2018_physeq_list,
  lamps_2022_physeq_list = lamps_2022_physeq_list
)

# Energy Farm Collab
# Examine the structure of your phyloseq list
str(main_physeq_list, max.level = 2)
names(main_physeq_list)
length(main_physeq_list)


# Explore phyloseq
purrr::iwalk(main_physeq_list, function(project_list, project_name) {
  cat("### PROJECT:", project_name, "###\n")
  purrr::iwalk(project_list, function(physeq_obj, region_name) {
    full_name <- paste(project_name, region_name, sep = "_")
    explore_phyloseq_list(physeq_obj, full_name)
  })
})

# Summary table
physeq_summary <- purrr::imap_dfr(
  main_physeq_list,
  function(project_list, project_name) {
    purrr::imap_dfr(
      project_list,
      function(physeq_obj, region_name) {
        data.frame(
          #physeq_list = project_name,
          region = gsub(".*_([^_]+)_physeq$", "\\1", region_name),
          n_taxa = ntaxa(physeq_obj),
          n_samples = nsamples(physeq_obj),
          total_reads = sum(sample_sums(physeq_obj)),
          min_reads_per_sample = min(sample_sums(physeq_obj)),
          max_reads_per_sample = max(sample_sums(physeq_obj)),
          mean_reads_per_sample = mean(sample_sums(physeq_obj)),
          median_reads_per_sample = median(sample_sums(physeq_obj))
        )
      },
      .id = "project_name"
    )
  },
  .id = "project_list"
)


physeq_summary


# Relative Abundace Transformation -----------------------------
transform_nested_to_relab <- function(nested_list) {
  purrr::map(nested_list, function(project_list) {
    purrr::map(project_list, function(physeq_obj) {
      transform_sample_counts(physeq_obj, function(x) x / sum(x))
    })
  })
}

main_relab_physeq_list <- transform_nested_to_relab(main_physeq_list)

# Verify transformations
purrr::map(main_physeq_relab_list, function(project_list) {
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

main_relab_df_list <- compute_diversity_nested(
  main_physeq_relab_list,
  drop = "ASV_",
  first_asv_col = "ASV_1"
)

# Save
save(
  main_relab_df_list,
  file = "data/output/processed/rdata/main_relab_df_list.rda"
)
