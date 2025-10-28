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

  # Read and combine sequences using read.fasta
  all_sequences_df <- purrr::map_dfr(fa_files, function(file) {
    phylotools::read.fasta(file) %>%
      dplyr::mutate(source_file = fs::path_file(file))
  })

  # Deduplicate sequences if processing all regions
  if (!is.null(.all)) {
    all_sequences_df <- all_sequences_df %>%
      dplyr::distinct(seq.text, .keep_all = TRUE)
  }

  # Rename ASVs
  all_sequences_df$seq.name <- paste0("ASV_", seq_len(nrow(all_sequences_df)))

  # Save to file using dat2fasta
  phylotools::dat2fasta(all_sequences_df, outfile = output_path)

  tryCatch(
    Biostrings::writeXStringSet(
      seqs,
      filepath = outfile,
      format = "fasta",
      append = FALSE,
      compress = FALSE
    ),
    error = function(e) {
      message("Failed to write ", outfile, ": ", conditionMessage(e))
    }
  )

  # Return summary info
  list(
    region = if (is.null(.all)) region else "all",
    n_files = length(fa_files),
    n_sequences = nrow(all_sequences_df),
    output_path = output_path
  )
}


# Export
results <- purrr::map(target_regions, process_fa) %>%
  set_names(target_regions)

result2 <- list(all_regions = process_fa(region = target_regions, .all = "all"))
purrr::map_dfr(result2, ~ data.frame(.x), .id = "region")

# View summary
purrr::map_dfr(result, ~ data.frame(.x), .id = "region")
