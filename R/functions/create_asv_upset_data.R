create_asv_upset_data <- function(physeq_list) {
  # Create base presence list
  asv_lists <- imap(physeq_list, function(psq, project_name) {
    # Get unique crops in this project
    crops <- sample_data(psq)$crop |> unique()

    # For each crop, get the ASVs present (non-zero reads)
    crop_asvs <- map(crops, function(crop_name) {
      # Subset phyloseq to this crop
      psq_subset <- prune_samples(sample_data(psq)$crop == crop_name, psq)

      # Extract ASVs with non-zero reads
      asvs_present <- taxa_names(psq_subset)[taxa_sums(psq_subset) > 0]

      asvs_present
    })

    if (grepl("^lamps", project_name, ignore.case = TRUE)) {
      # LAMPS:  extract first two parts (lamps_2018 or lamps_2022)
      crop_asvs <- set_names(
        crop_asvs,
        paste0(
          gsub("^([^_]+_[^_]+).*", "\\1", project_name),
          "_",
          crops
        )
      )
    } else {
      # (ef, mxg, etc. )
      crop_asvs <- set_names(
        crop_asvs,
        paste0(gsub("^([^_]+).*", "\\1", project_name), "_", crops)
      )
    }

    crop_asvs
  })
  # Extract taxonomic information from the first phyloseq object
  first_psq <- physeq_list[[1]]
  tax_table_df <- tax_table(first_psq) |>
    as.data.frame() |>
    rownames_to_column("asv")

  # Create a detailed attributes dataframe using imap + list_rbind
  asv_attributes <- imap(asv_lists, function(project_crops, project_name) {
    imap(project_crops, function(asvs, crop_name) {
      tibble(
        asv = asvs,
        project = project_name,
        crop = crop_name
      ) |>
        left_join(tax_table_df, by = "asv")
    }) |>
      list_rbind(names_to = "crop_combo")
  }) |>
    list_rbind(names_to = "project_id")

  list(
    presence_list = asv_lists,
    attributes = asv_attributes
  )
}
