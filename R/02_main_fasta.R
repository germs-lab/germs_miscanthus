###################################################################
# From Phyloseq to FASTA
# Author: Bolívar Aponte Rolón
# Date: 2025-10-23
###################################################################

source("R/utils/00_setup.R")


# Subset phyloseq by plant type
# To subset the actual phyloseq object
mxg_lamps_2018 <- purrr::map(
  lamps_2018_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, plant == "M")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
)

mxg_lamps_2022 <- purrr::map(
  lamps_2022_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, treatment == "Miscanthus")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
)

mxg_ef <- purrr::map(
  ef_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "MXG")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
)

# Phyloseq to FASTA
# create output dir

# LAMPS
refseq2fasta(
  mxg_lamps_2018,
  extra_id = "_mxg",
  out_dir = "data/output/processed/sequences"
)
refseq2fasta(
  mxg_lamps_2022,
  extra_id = "_mxg",
  out_dir = "data/output/processed/sequences"
)


# Energy Farm Collab
refseq2fasta(
  mxg_ef,
  extra_id = "_mxg",
  out_dir = "data/output/processed/sequences"
)


# Concatenate and export per target region
target_regions <- c("16S", "ITS", "AMF")

process_fa <- function(region, .all = NULL) {
  # Determine file selection logic
  if (is.null(.all)) {
    # Process specific region
    fa_files <- fs::dir_ls("data/output/processed/sequences/") %>%
      str_subset(paste0(region, "_mxg.fa"))

    output_path <- paste0(
      "data/output/processed/sequences/mxg_",
      region,
      "_combined_asv_renamed.fa"
    )
  } else {
    # Process all regions combined
    fa_files <- fs::dir_ls("data/output/processed/sequences/") %>%
      str_subset("_mxg.fa")

    output_path <- paste0(
      "data/output/processed/sequences/mxg_",
      "all",
      "_combined_asv_renamed.fa"
    )
  }

  # Read and combine sequences
  fa_list <- purrr::map(fa_files, readr::read_lines) %>%
    set_names(fs::path_file(fa_files))

  all_sequences <- unlist(fa_list)

  # Deduplicate sequences if processing all regions
  if (!is.null(.all)) {
    all_sequences <- deduplicate_sequences(all_sequences)
  }

  # Rename ASVs
  header_indices <- which(stringr::str_starts(all_sequences, ">"))
  new_asv_names <- paste0("ASV_", seq_along(header_indices))
  all_sequences[header_indices] <- paste0(">", new_asv_names)

  # Save to file
  readr::write_lines(all_sequences, output_path)

  # Return summary info
  list(
    region = if (is.null(.all)) region else "all",
    n_files = length(fa_files),
    n_sequences = length(header_indices),
    output_path = output_path
  )
}

# Helper function for deduplication
deduplicate_sequences <- function(all_sequences) {
  header_indices <- which(stringr::str_starts(all_sequences, ">"))

  # Create data frame directly with sequence pairs
  sequence_df <- data.frame(
    header = character(length(header_indices)),
    sequence = character(length(header_indices)),
    stringsAsFactors = FALSE
  )

  sequence_pairs <- purrr::map_dfr(seq_along(header_indices), function(i) {
    start_pos <- header_indices[i]
    end_pos <- if (i < length(header_indices)) {
      header_indices[i + 1] - 1
    } else {
      length(all_sequences)
    }
    sequence_df$header[i] <- all_sequences[start_pos]
    sequence_df$sequence[i] <- paste(
      all_sequences[(start_pos + 1):end_pos],
      collapse = ""
    )
  })

  unique_pairs <- sequence_pairs %>%
    dplyr::distinct(sequence, .keep_all = TRUE)

  as.vector(rbind(unique_pairs$header, unique_pairs$sequence))
}

# Export
results <- purrr::map(target_regions, process_fa) %>%
  set_names(target_regions)

result2 <- list(all_regions = process_fa(region = target_regions, .all = "all"))
purrr::map_dfr(result2, ~ data.frame(.x), .id = "region")

# View summary
purrr::map_dfr(result, ~ data.frame(.x), .id = "region")
