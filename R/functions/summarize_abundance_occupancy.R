summarize_abundance_occupancy <- function(physeq_obj, crop_column = "crop") {
  results <- list()

  # Extract OTU table (ASVs as rows, samples as columns)
  otu_mat <- as.data.frame(otu_table(physeq_obj))

  overall_df <- otu_mat %>%
    reframe(
      abundance = rowSums(.),
      occupancy = 1 * rowSums(. > 0),
      occupancy_prop = occupancy / ncol(.), #aka Prevalence
      mean_abundance = abundance / occupancy,
      mean_rel_abundance = apply(
        vegan::decostand(., method = "total", MARGIN = 2),
        1,
        mean
      )
    )
  rownames(overall_df) <- rownames(
    otu_mat
  )

  # Selecting crops
  crops <- sample_data(physeq_obj) %>%
    as.matrix() %>%
    as.data.frame() %>%
    pull(crop_column) %>%
    unique()

  by_crop_list <- purrr::map(crops, function(crop_name) {
    # Index
    sample_idx <- which(sample_data(physeq_obj)[[crop_column]] == crop_name)

    if (length(sample_idx) == 0) {
      return(NULL) # Skip if no samples
    }

    # Subset and clean
    otu_mat_by_crop <- prune_samples(
      sample_names(physeq_obj)[sample_idx],
      physeq_obj
    ) %>%
      phyloseq::prune_taxa(taxa_sums(.) > 0, .) %>%
      otu_table() %>%
      as.data.frame()

    # Summary by crop
    crop_metrics <- otu_mat_by_crop %>%
      reframe(
        abundance = rowSums(.),
        occupancy = 1 * rowSums(. > 0),
        occupancy_prop = occupancy / ncol(.), #aka Prevalence
        mean_abundance = abundance / occupancy,
        mean_rel_abundance = apply(
          vegan::decostand(., method = "total", MARGIN = 2),
          1,
          mean
        )
      )

    rownames(crop_metrics) <- rownames(
      otu_mat_by_crop
    )

    return(crop_metrics)
  })

  names(by_crop_list) <- crops

  # Results
  results <- list(
    overall = overall_df,
    by_crop = by_crop_list
  )

  return(results)
}
