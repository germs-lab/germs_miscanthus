create_binary_df_from_flat <- function(flat_list) {
  all_items <- unique(unlist(flat_list))
  set_names <- names(flat_list)
  binary_df <- map_dfc(set_names, function(set_name) {
    tibble(!!set_name := as.integer(all_items %in% flat_list[[set_name]]))
  })
  binary_df$item <- all_items
  binary_df |> relocate(item)
}

list_to_binary_df <- function(presence_list) {
  flattened <- purrr::flatten(presence_list)
  create_binary_df_from_flat(flattened)
}

create_upset_df <- function(presence_list, attributes_df, tax_level = "asv") {
  binary_df <- list_to_binary_df(presence_list)

  tax_info <- attributes_df |>
    select(asv, phylum, class, genus) |>
    distinct()

  if (tax_level == "asv") {
    upset_df <- binary_df |>
      rename(asv = item) |>
      left_join(tax_info, by = "asv")
  } else {
    upset_df <- binary_df |>
      rename(!!tax_level := item)
  }

  upset_df
}

create_tax_upset_df <- function(presence_list) {
  binary_df <- list_to_binary_df(presence_list) |>
    rename(taxon = item)
  binary_df
}

create_tax_presence_list <- function(physeq_list, tax_level = "phylum") {
  imap(physeq_list, function(psq, project_name) {
    crops <- sample_data(psq)$crop |> unique()

    crop_taxa <- map(crops, function(crop_name) {
      psq_subset <- prune_samples(sample_data(psq)$crop == crop_name, psq)
      tax_df <- tax_table(psq_subset) |>
        as.data.frame() |>
        rownames_to_column("asv")

      present_asvs <- taxa_names(psq_subset)[taxa_sums(psq_subset) > 0]
      tax_df |>
        filter(asv %in% present_asvs) |>
        pull(!!sym(tax_level)) |>
        unique() |>
        na.omit() |>
        as.character()
    })

    if (grepl("^lamps", project_name, ignore.case = TRUE)) {
      prefix <- gsub("^([^_]+_[^_]+).*", "\\1", project_name)
    } else {
      prefix <- gsub("^([^_]+).*", "\\1", project_name)
    }
    set_names(crop_taxa, paste0(prefix, "_", crops))
  })
}

get_mxg_sets <- function(presence_list) {
  flattened <- purrr::flatten(presence_list)
  mxg_patterns <- c("MXG", "Miscanthus", "_M$")
  mxg_idx <- grep(paste(mxg_patterns, collapse = "|"), names(flattened))
  flattened[mxg_idx]
}
