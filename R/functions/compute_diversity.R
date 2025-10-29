compute_diversity <- function(
  dataset,
  drop,
  first_asv_col = NULL
) {
  # Abundance block (ASVs) as numeric matrix
  abund <- dataset %>% select(starts_with(drop))
  abund_mat <- as.matrix(abund)

  out <- dataset %>%
    mutate(
      observed = rowSums(abund_mat > 0, na.rm = TRUE),
      shannon = vegan::diversity(abund_mat, index = "shannon"),
      simpson = vegan::diversity(abund_mat, index = "simpson"),
      inv_simpson = vegan::diversity(abund_mat, index = "invsimpson")
    )

  # Figure out where to relocate the new columns
  if (is.null(first_asv_col)) {
    # put before the first ASV column after the dropped block
    first_asv_col <- colnames(dataset)[min(setdiff(
      seq_len(ncol(dataset)),
      drop
    ))]
  }

  out <- out %>%
    relocate(
      any_of(c("observed", "shannon", "simpson", "inv_simpson")),
      .before = all_of(first_asv_col)
    )

  return(out)
}


compute_diversity_nested <- function(
  nested_list,
  drop,
  first_asv_col = NULL
) {
  results <- list()

  # Level 1: project groups within each main category
  purrr::iwalk(nested_list, function(level1_list, level1_name) {
    # Level 2: individual phyloseq objects
    purrr::iwalk(level1_list, function(dataset, level2_name) {
      # Create unique name for flattened structure
      result_name <- paste(level2_name, "df", sep = "_")
      result_name <- gsub("_physeq", "", result_name) # clean up naming
      # Phyloseq to DF
      physeq_df <- physeq2df(dataset)
      # Apply diversity computation
      results[[result_name]] <<- compute_diversity(
        dataset = physeq_df,
        drop = drop,
        first_asv_col = first_asv_col
      )
    })
  })

  return(results)
}
