export_otu_table <- function(physeq_obj) {
  result <- physeq_obj %>%
    otu_table() %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(., var = "sample_id")

  return(result)
}
