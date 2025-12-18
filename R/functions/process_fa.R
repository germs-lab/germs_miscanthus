#' Process and Combine FASTA Files
#'
#' Reads existing FASTA files, combines them, removes duplicates, and exports
#' with explicit naming indicating the combination parameters.
#'
#' @param region Character vector of target regions (e.g., "16S", "ITS", "AMF").
#' @param path Directory containing input FASTA files.
#' @param vec_files Optional. Character string to filter files by name pattern.
#' @param crop_subset Character string indicating crop subset in source files
#'   (e.g., "mxg_only", "all_crops"). Used to filter input files.
#' @param output_prefix Character prefix for output file (e.g., "combined_mxg_").
#' @param vec_files Optional. Vector of targete files to concatenate.
#' @param nucleotide Optional. Nucleotide type to filter by (e.g., "DNA", "RNA").
#' @param exclude_str Optional. Pattern to exclude files (used with combine_all = TRUE).
#' @param combine_all Logical. If TRUE, combines all specified regions into one file.
#'   Default is FALSE (process each region separately).
#' @param rename_headers Logical. If TRUE, renames sequence headers to ASV_1, ASV_2, etc.
#'   Default is FALSE.
#'
#' @return A list containing:
#'   - region: The region(s) processed
#'   - n_files: Number of input files processed
#'   - n_sequences: Number of unique sequences in output
#'   - output_path: Path to the output file
#'
#' @details
#' Output naming convention:
#' - Single region: {output_prefix}{region}_{nucleotide}_combined.fa
#' - All regions: {output_prefix}all_regions_{nucleotide}_combined.fa
#'
#' Examples:
#' - combined_mxg_16S_DNA_combined.fa (Miscanthus 16S DNA combined)
#' - combined_all_crops_AMF_DNA_combined.fa (All crops AMF DNA combined)
#' - combined_mxg_all_regions_DNA_combined.fa (Miscanthus all regions DNA combined)
#'
#' @export
process_fa <- function(
  region,
  path = "data/output/sequences/",
  vec_files = NULL,
  crop_subset = NULL,
  output_prefix = "combined_",
  extra_id = NULL,
  nucleotide = NULL,
  exclude_str = NULL,
  combine_all = FALSE,
  rename_headers = FALSE
) {
  # Find matching files
  fa_files <- fs::dir_ls(path, regexp = "\\.fa$")

  # Filter by region - handle vector of regions with grepl pattern
  if (length(region) > 1) {
    # Create regex pattern for multiple regions: (16S|ITS|AMF)
    region_pattern <- paste0("(", paste(region, collapse = "|"), ")")
    fa_files <- fa_files[grepl(region_pattern, basename(fa_files))]
  } else {
    # Single region
    fa_files <- fa_files[grepl(region, basename(fa_files))]
  }
  # Exclude pattern applies when combine_all = TRUE (to avoid mixing RNA with DNA)
  if (combine_all && !is.null(exclude_str)) {
    fa_files <- fa_files[!grepl(exclude_str, basename(fa_files))]
  }

  # Filter by crop subset if provided
  if (!is.null(crop_subset)) {
    fa_files <- fa_files[grepl(crop_subset, basename(fa_files))]
  }

  # Filter by nucleotide if provided
  if (!is.null(nucleotide)) {
    fa_files <- fa_files[grepl(nucleotide, basename(fa_files))]
  }

  if (!is.null(vec_files)) {
    vec_pattern <- paste0("(", paste(vec_files, collapse = "|"), ")")
    fa_files <- fa_files[grepl(vec_pattern, basename(fa_files))]
  }

  if (length(fa_files) == 0) {
    warning(
      "No FASTA files found matching criteria:\n",
      "  region: ",
      paste(region, collapse = ", "),
      "\n",
      "  crop_subset: ",
      if (is.null(crop_subset)) "any" else crop_subset,
      "\n",
      "  nucleotide: ",
      if (is.null(nucleotide)) "any" else nucleotide
    )
    return(NULL)
  }

  if (length(fa_files) == 1) {
    message(
      "Only 1 file found. Proceeding with single file (no deduplication needed)."
    )
  }

  # Print check of target files
  cat("Target files to combine:\n")
  cat(paste0("  - ", basename(fa_files), "\n"))
  cat("Total files:", length(fa_files), "\n\n")

  # Build output filename
  region_label <- if (combine_all) "all_regions" else region
  extra_id <- if (!is.null(extra_id)) paste0("_", extra_id) else ""
  nuc_label <- if (!is.null(nucleotide)) paste0("_", nucleotide) else ""
  crop_label <- if (!is.null(crop_subset)) paste0("_", crop_subset) else ""

  if (length(fa_files) == 1) {
    output_filename <- basename(fa_files)
    output_path <- fa_files[[1]]
    message("Overwriting file of the same name.")
  } else {
    output_filename <- paste0(
      output_prefix,
      region_label,
      extra_id,
      nuc_label,
      crop_label,
      "_combined.fa"
    )
    output_path <- file.path(path, output_filename)
  }
  # Read and combine sequences
  all_sequences_df <- purrr::map_dfr(fa_files, function(file) {
    phylotools::read.fasta(file) %>%
      dplyr::mutate(source_file = fs::path_file(file))
  })

  # Deduplicate sequences (keep first occurrence)
  unique_df <- all_sequences_df[!duplicated(all_sequences_df$seq.text), ]

  # Rename headers if requested
  if (rename_headers) {
    unique_df$seq.name <- paste0("ASV_", seq_len(nrow(unique_df)))
  }

  # Export to file
  phylotools::dat2fasta(unique_df, outfile = output_path)

  message(
    "Combined ",
    length(fa_files),
    " files -> ",
    output_filename,
    " (",
    nrow(unique_df),
    " unique sequences from ",
    nrow(all_sequences_df),
    " total)"
  )

  # Return summary
  list(
    region = region_label,
    n_files = length(fa_files),
    n_sequences = nrow(unique_df),
    n_total_sequences = nrow(all_sequences_df),
    input_files = basename(fa_files),
    output_path = output_path
  )
}
