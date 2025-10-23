build_updated_metadata <- function(
  csv,
  xlsx,
  .distinct = TRUE,
  target_region,
  rna = FALSE,
  .project
) {
  # Read and clean CSV metadata
  meta_csv <- read_csv(
    csv,
    col_types = cols(.default = col_character()),
    show_col_types = TRUE
  ) %>%
    janitor::clean_names() %>%
    select(sample_id, sample_date, soil_type) %>%
    as_tibble()

  # Read and clean Excel metadata
  meta_xls <- read_xlsx(
    xlsx,
    sheet = "useful_metadata",
    .name_repair = "unique",
    col_types = "text"
  ) %>%
    janitor::clean_names() %>%
    rename(sample_id = id) %>%
    mutate(
      # Standardize plant names to tokens
      plant = case_when(
        plant == "Corn" ~ "C",
        plant == "Miscanthus" ~ "M",
        TRUE ~ plant
      ),
      # Add prefixes to plot and nitrogen concentration
      plot = paste0("P", plot),
      nitrogen_conc = paste0("N", nitrogen_conc)
    ) %>%
    select(!samples)

  # Join metadata and create sequence IDs
  meta <- meta_xls %>%
    left_join(., meta_csv, by = "sample_id") %>%
    relocate(c(year:nitrogen_conc), .after = plant) %>%
    relocate(sample_date, .after = replicate) %>%
    mutate(
      # Normalize re-extraction marker
      sample_id = if_else(
        str_starts(sample_id, "Reextracted-"),
        str_replace(sample_id, "^Reextracted-", "re-"),
        sample_id
      ),
      # For re- samples, copy sample_date from the original sample
      sample_date = case_when(
        str_starts(sample_id, "re-") & is.na(sample_date) ~
          {
            # Extract the number after "re-"
            original_id <- str_extract(sample_id, "(?<=re-)\\d+")
            # Find the sample_date from the matching original sample
            sample_date[match(original_id, sample_id)]
          },
        TRUE ~ sample_date
      ),
      target_region = target_region,
      project = .project
    ) %>%
    relocate(target_region, .after = nucleotide) %>%
    tidyr::unite(
      "sequence_id",
      c(sample_id, plant, year, plot, nitrogen_conc, replicate, sample_date),
      sep = "_",
      remove = FALSE,
      na.rm = TRUE
    ) %>%
    relocate(sequence_id, .before = sample_id)

  # RNA-specific processing
  if (rna) {
    meta <- meta %>%
      mutate(
        # Fix potential typo in sequence file naming
        sample_date = case_when(
          str_equal(sample_date, "20180429") ~ "20180430",
          TRUE ~ sample_date
        )
      ) %>%
      tidyr::unite(
        "sequence_id",
        c(
          sample_id,
          plant,
          year,
          plot,
          nitrogen_conc,
          replicate,
          sample_date
        ),
        sep = "_",
        remove = FALSE,
        na.rm = TRUE
      ) %>%
      relocate(sequence_id, .before = sample_id) %>%
      mutate(
        # Add RNA suffix to sequence_id
        sequence_id = paste(sequence_id, "RNA", sep = "_")
      )
  }

  # Remove duplicates if requested
  if (.distinct) {
    meta <- meta %>%
      distinct(sequence_id, .keep_all = TRUE)
  }

  # Data validation checks
  stopifnot(!any(is.na(meta$sample_id)))

  if (anyDuplicated(meta$sample_id)) {
    warning("Duplicate sample_id detected.")
  }

  if (any(is.na(meta$plant))) {
    warning("Missing plant for some rows after join.")
  }

  return(meta)
}
