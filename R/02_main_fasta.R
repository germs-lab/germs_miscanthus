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

# rownames(taxa_16S) <- paste0("ASV_", 1:nrow(taxa_16S))

# asv_fasta <- seqtab2fasta(seqtab_nochim_16S)

# seqtab_nochim_16S <- t(seqtab_nochim_16S) # Retaining sequences and asigning shorthand ASV names

# row.names(seqtab_nochim_16S) <- sub(">", "", asv_fasta$asv_headers)

# asv_16S <- otu_table(seqtab_nochim_16S, taxa_are_rows = TRUE)

# ## Save FASTA file
# dir.create("data/output/processed/sequences/", recursive = TRUE)
# write(
#   asv_fasta$asv_fasta,
#   file.path("data/output/processed/sequences/energy_farm_collab_16S.fa")
# )
