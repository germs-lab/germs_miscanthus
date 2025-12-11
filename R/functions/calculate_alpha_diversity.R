# Calculate alpha diversity metrics for nested phyloseq objects

calculate_alpha_diversity <- function(project_list) {
  purrr::imap(project_list, function(physeq_obj, physeq_name) {
    # Calculate alpha diversity indices
    alpha_div_full <- alpha_workflow(physeq_obj = physeq_obj)

    # Add project and region information
    alpha_div_full <- alpha_div_full %>%
      mutate(
        project = gsub("([^_]+).*", "\\1", physeq_name),
        physeq_name = physeq_name,
        sample_date = if ("sample_date" %in% names(.)) {
          as.character(sample_date)
        } else {
          NA_character_
        }
      ) %>%
      relocate(c(observed:inv_simpson), .before = project)

    return(alpha_div_full)
  })
}


calculate_alpha_diversity_nested <- function(nested_list) {
  purrr::imap(nested_list, function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      # Calculate alpha diversity indices
      # physeq_obj <- main_mxg_physeq_list$mxg_ef$ef_16S_physeq
      alpha_div_full <- alpha_workflow(physeq_obj = physeq_obj)

      # Add project and region information
      alpha_div_full <- alpha_div_full %>%
        mutate(
          project = project_name,
          region = gsub(".*_([^_]+)_physeq$", "\\1", physeq_name),
          physeq_name = physeq_name,
          sample_date = if ("sample_date" %in% names(.)) {
            as.character(sample_date)
          } else {
            NA_character_
          }
        ) %>%
        relocate(c(observed:inv_simpson), .before = project)

      return(alpha_div_full)
    })
  })
}


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
