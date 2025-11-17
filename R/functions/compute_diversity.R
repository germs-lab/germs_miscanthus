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
  results <- purrr::imap(nested_list, function(level1_list, level1_name) {
    # Level 1: project groups within each main category

    purrr::imap(level1_list, function(dataset, level2_name) {
      # Level 2: individual phyloseq objects

      # Phyloseq to DF
      physeq_df <- physeq2df(dataset)

      compute_diversity(
        dataset = physeq_df,
        drop = drop,
        first_asv_col = first_asv_col
      )
    })
  })

  # Flatten the nested structure and set names
  results <- purrr::list_flatten(results, name_spec = "{inner}_df") %>%
    purrr::set_names(gsub("_physeq", "", names(.)))

  return(results)
}
