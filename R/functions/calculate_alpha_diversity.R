# Calculate alpha diversity metrics for nested phyloseq objects

alpha_workflow <- function(physeq_obj) {
  alpha_div <- estimate_richness(
    physeq_obj,
    measures = c("Observed", "Shannon", "Simpson", "InvSimpson")
  ) %>%
    janitor::clean_names() %>%
    rownames_to_column(var = "sample_id") %>%
    mutate(
      # estimate_richness() messes up sample names, now we need to clean them
      sample_id = str_replace_all(sample_id, c("^X" = "", "\\." = "-"))
    ) %>%
    column_to_rownames(var = "sample_id")

  # Add sample data
  sample_df <- data.frame(sample_data(physeq_obj))
  alpha_div$sample_id <- rownames(alpha_div)

  # Merge with sample metadata
  alpha_div_full <- merge(
    alpha_div,
    sample_df,
    by.x = "sample_id",
    by.y = "row.names"
  )
  return(alpha_div_full)
}
