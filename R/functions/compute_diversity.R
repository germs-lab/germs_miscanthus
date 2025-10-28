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
