export_otu_table <- function(physeq_obj, taxa_as_rows = FALSE) {
  if (taxa_as_rows) {
    result <- physeq_obj %>%
      otu_table() %>%
      as.data.frame() %>%
      rownames_to_column(., var = "taxa_id")
  } else {
    result <- physeq_obj %>%
      otu_table() %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(., var = "sample_id")
  }

  return(result)
}
