###################################################################
# From Phyloseq to FASTA
# Author: Bolívar Aponte Rolón
# Date: 2025-10-23
###################################################################

source("R/utils/00_setup.R")

# Phyloseq to FASTA
# create output dir

# LAMPS 2022
refseq2fasta(
  lamps_2022_physeq_list,
  out_dir = "data/output/processed/sequences"
)


# Energy Farm Collab
refseq2fasta(ef_physeq_list, out_dir = "data/output/processed/sequences")

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
