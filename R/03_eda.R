###################################################################
# EDA for Microbiome of MIscanthus
#
#
#
#Author: Bolívar Aponte Rolón
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
purrr::iwalk(main_physeq_list, function(study_list, study_name) {
  cat("### STUDY:", study_name, "###\n")
  purrr::iwalk(study_list, function(physeq_obj, region_name) {
    full_name <- paste(study_name, region_name, sep = "_")
    explore_phyloseq(physeq_obj, full_name)
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
  purrr::map(nested_list, function(study_list) {
    purrr::map(study_list, function(physeq_obj) {
      transform_sample_counts(physeq_obj, function(x) x / sum(x))
    })
  })
}

main_physeq_relab_list <- transform_nested_to_relab(main_physeq_list)

# Verify transformations
purrr::map(main_physeq_relab_list, function(study_list) {
  purrr::map(study_list, function(physeq_obj) {
    sample_sums(physeq_obj)
  })
})

# Phyloseq to DF

main_physeq_relabdf_list <- purrr::map(
  main_physeq_relab_list,
  function(study_list) {
    purrr::map(study_list, function(physeq_obj) {
      physeq2df(physeq_obj)
    })
  }
)

# Calculate Diversity indices ---------------------------------
#TODO
test <- compute_diversity(
  main_physeq_relabdf_list$ef_physeq_list$ef_16S_physeq,
  drop = "ASV_",
  first_asv_col = "ASV_1"
)


# Function to apply and flatten results
compute_diversity_nested <- function(
  nested_list,
  drop,
  first_asv_col = NULL
) {
  results <- list()

  # Level 1: main categories (ef_physeq_list, lamps_2018_physeq_list, etc.)
  purrr::iwalk(nested_list, function(level1_list, level1_name) {
    # Level 2: study groups within each main category
    purrr::iwalk(level1_list, function(level2_list, level2_name) {
      # Level 3: individual phyloseq objects
      purrr::iwalk(level2_list, function(dataset, level3_name) {
        # Create unique name for flattened structure
        result_name <- paste(level1_name, level2_name, level3_name, sep = "_")
        result_name <- gsub("_physeq", "", result_name) # clean up naming

        # results <- physeq2df(physeq_obj)
        # Apply diversity computation (assuming dataset is already a data frame)
        results[[result_name]] <<- compute_diversity(
          dataset = dataset,
          drop = drop,
          first_asv_col = first_asv_col
        )
      })
    })
  })

  return(results)
}

# Usage:
flat_diversity_results <- compute_diversity_nested(
  main_physeq_relabdf_list,
  drop = "ASV_",
  first_asv_col = "ASV_1"
)
test_single <- compute_diversity(
  main_physeq_relabdf_list$ef_physeq_list$ef_16S_physeq,
  drop = "ASV_",
  first_asv_col = "ASV_1"
)

# Add some debug prints to see where it fails
compute_diversity_nested_debug <- function(
  nested_list,
  drop,
  first_asv_col = NULL
) {
  results <- list()

  purrr::iwalk(nested_list, function(level1_list, level1_name) {
    cat("Processing level1:", level1_name, "\n")
    purrr::iwalk(level1_list, function(level2_list, level2_name) {
      cat("  Processing level2:", level2_name, "\n")
      purrr::iwalk(level2_list, function(dataset, level3_name) {
        cat(
          "    Processing level3:",
          level3_name,
          "- class:",
          class(dataset)[1],
          "\n"
        )
        # Just return the class instead of processing for now
        result_name <- paste(level1_name, level2_name, level3_name, sep = "_")
        results[[result_name]] <<- class(dataset)[1]
      })
    })
  })

  return(results)
}

# Run the debug version
debug_results <- compute_diversity_nested_debug(
  main_physeq_relabdf_list,
  drop = "ASV_",
  first_asv_col = "ASV_1"
)

# main_dataframes <- imap(
#   datasets,
#   ~ compute_diversity(.x, drop = unwanted_cols, first_asv_col = "ASV_1")
# )
