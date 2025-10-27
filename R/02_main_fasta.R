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


# Concatenate and export as single file
fa_mxg_files <- fs::dir_ls("data/output/processed/sequences/") %>%
  str_subset("_mxg.fa")

main_fa_list <- purrr::map(fa_mxg_files, readr::read_lines) %>%
  set_names(fs::path_file(fa_mxg_files))

all_sequences <- unlist(main_fa_list)

# Find header lines (start with ">")
header_indices <- which(stringr::str_starts(all_sequences, ">"))

new_asv_names <- paste0("ASV_", seq_along(header_indices))

# Replace headers
all_sequences[header_indices] <- paste0(">", new_asv_names)

# Save to file
readr::write_lines(
  all_sequences,
  "data/output/processed/sequences/combined_mxg_seqs_renamed.fa"
)
