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
    ps_subset <- subset_samples(.x, crop == "M")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
)

mxg_lamps_2022 <- purrr::map(
  lamps_2022_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "Miscanthus")
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

# Main list for downstream analyses
# Main list for downstream analyses
main_physeq_list <- list(
  ef_physeq_list = ef_physeq_list,
  lamps_2018_physeq_list = lamps_2018_physeq_list,
  lamps_2022_physeq_list = lamps_2022_physeq_list
)

save(
  main_physeq_list,
  file = "data/output/processed/rdata/main_physeq_list.rda"
)

main_mxg_physeq_list <- list(
  mxg_ef = mxg_ef,
  mxg_lamps_2018 = mxg_lamps_2018,
  mxg_lamps_2022 = mxg_lamps_2022
)

save(
  main_mxg_physeq_list,
  file = "data/output/processed/rdata/main_mxg_physeq_list.rda"
)

# Phyloseq to FASTA
# create output dir

# LAMPS
# Full
refseq2fasta(
  lamps_2018_physeq_list,
  out_dir = "data/output/processed/sequences"
)

refseq2fasta(
  lamps_2022_physeq_list,
  out_dir = "data/output/processed/sequences"
)

# MXG
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
# Full
refseq2fasta(
  ef_physeq_list,
  out_dir = "data/output/processed/sequences"
)

# MXG
refseq2fasta(
  mxg_ef,
  extra_id = "_mxg",
  out_dir = "data/output/processed/sequences"
)


# Concatenate and export per target region
target_regions <- c("16S", "ITS", "AMF")

# Export and summary
results <- purrr::map(target_regions, process_fa) %>%
  set_names(target_regions)
purrr::map_dfr(results, ~ data.frame(.x), .id = "region")

# All together
result <- list(all_regions = process_fa(region = target_regions, .all = "all"))
purrr::map_dfr(result, ~ data.frame(.x), .id = "region")
