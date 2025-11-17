physeq2df <- function(physeq) {
  conditional_df <- if ("sample_id" %in% names(phyloseq::sample_data(physeq))) {
    physeq %>% phyloseq::sample_data() %>% data.frame()
  } else {
    physeq %>%
      phyloseq::sample_data() %>%
      data.frame() %>%
      rownames_to_column(var = "sample_id")
  }

  physeq_df <- physeq %>%
    otu_table() %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(., var = "iter_id") %>%
    dplyr::left_join(
      .,
      conditional_df,
      by = c("iter_id" = "sample_id")
    ) %>%
    mutate(sample_id = iter_id) %>%
    relocate(sample_id, .after = iter_id) %>%
    relocate(., !starts_with("ASV_"), .after = sample_id) %>%
    column_to_rownames(., var = "iter_id")

  return(physeq_df)
}
