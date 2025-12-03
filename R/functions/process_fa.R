process_fa <- function(
  region,
  path = "data/output/sequences/",
  combined_suffix = "_combined_asv_renamed.fa",
  prefix = NULL,
  target_suffix = NULL,
  .all = NULL
) {
  # Determine file selection logic
  if (is.null(.all)) {
    # Process specific region
    fa_files <- fs::dir_ls(path) %>%
      str_subset(paste0(region, target_suffix))

    output_path <- paste0(
      path,
      prefix,
      region,
      combined_suffix
    )
  } else {
    # Process all regions combined
    fa_files <- fs::dir_ls(path) %>%
      str_subset(target_suffix)

    output_path <- paste0(
      path,
      prefix,
      .all,
      combined_suffix
    )
  }

  # Read and combine sequences using read.fasta
  all_sequences_df <- purrr::map_dfr(fa_files, function(file) {
    phylotools::read.fasta(file) %>%
      dplyr::mutate(source_file = fs::path_file(file))
  })

  # Deduplicate sequences
  unique_df <- all_sequences_df[!duplicated(all_sequences_df$seq.text), ]

  # Rename ASVs
  unique_df$seq.name <- paste0("ASV_", seq_len(nrow(unique_df)))

  # Save to file using dat2fasta
  phylotools::dat2fasta(unique_df, outfile = output_path)

  # Return summary info
  list(
    region = if (is.null(.all)) region else "all",
    n_files = length(fa_files),
    n_sequences = nrow(unique_df),
    output_path = output_path
  )
}
