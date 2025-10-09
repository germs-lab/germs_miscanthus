###################################################################
# Import of R object from DADA2 pipeline for phyloseq generation
#
#
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-01
###################################################################

source("R/utils/00_setup.R")

#===============================================
# Files from Energy Farm Collab
#===============================================

# Metadata ------------------------
# For 16S and AMF
# "data/input/energy_farm_collab/files_for_phyloseq_16S/ef2024_sampledata.csv" and
# "data/input/energy_farm_collab/files_for_phyloseq_AMF/ef2024_sampledata.csv" are the same

ef_metadata <- read.csv("data/input/energy_farm_collab/files_for_phyloseq_16S/ef2024_sampledata.csv") %>%
  janitor::clean_names(.) %>%
  column_to_rownames(var = "label_id")


ef_metadata <- sample_data(ef_metadata)

#----------------------------------
# 16S
#----------------------------------

# Taxonomy ------------------------
taxa_16S <- readRDS(
  "data/input/energy_farm_collab/files_for_phyloseq_16S/taxa.rds"
) %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sequence") %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv

rownames(taxa_16S) <- paste0("ASV_", 1:nrow(taxa_16S))

taxa_16S <- tax_table(as.matrix(taxa_16S))

# ASVs ------------------------
seqtab_nochim_16S <- readRDS(
  "data/input/energy_farm_collab/files_for_phyloseq_16S/seqtab.nochim.rds"
)
# Create ASV table
## Cleaning ASV names for FASTA file

asv_fasta <- seqtab2fasta(seqtab_nochim_16S)

seqtab_nochim_16S <- t(seqtab_nochim_16S) # Retaining sequences and asigning shorthand ASV names

row.names(seqtab_nochim_16S) <- sub(">", "", asv_fasta$asv_headers)

asv_16S <- otu_table(seqtab_nochim_16S, taxa_are_rows = TRUE)


## Save FASTA file
dir.create("data/output/processed/sequences/", recursive = TRUE)
write(
  asv_fasta$asv_fasta,
  file.path("data/output/processed/sequences/energy_farm_collab_16S.fa")
)

# 16S Phyloseq object ------------------------
# Checking that metadata and asvs have the same number of samples
samples_missing_metadata <- base::setdiff(
  colnames(seqtab_nochim_16S),
  rownames(ef_metadata)
)
# We have 10 missing samples

ef_16S_physeq <- phyloseq(asv_16S, taxa_16S, ef_metadata)

otu_table(ef_16S_physeq)
tax_table(ef_16S_physeq)
sample_data(ef_16S_physeq)

# Any other missing or dropped sample?
post_physeq_missing <- base::setdiff(
  sample_names(ef_metadata),
  sample_names(ef_16S_physeq)
)
missing_sample <- post_physeq_missing %in% colnames(seqtab_nochim_16S)


#----------------------------------
# AMF (only forward reads)
#----------------------------------

# Taxonomy ------------------------
taxa_AMF <- readRDS(
  "data/input/energy_farm_collab/files_for_phyloseq_AMF/taxa.rds"
) %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sequence") %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv

rownames(taxa_AMF) <- paste0("ASV_", 1:nrow(taxa_AMF))

taxa_AMF <- tax_table(as.matrix(taxa_AMF))

# ASVs ------------------------
seqtab_nochim_AMF <- readRDS("data/input/energy_farm_collab/files_for_phyloseq_AMF/seqtabf.nochim.rds") # Forward reads only

# seqtab_nochim_AMF needs a sample name clean up to match metadata
rownames(seqtab_nochim_AMF) <-gsub("_S.*", "", rownames(seqtab_nochim_AMF))

# Remove "method" samples 
seqtab_nochim_AMF <- seqtab_nochim_AMF[!grepl("^method([1-9]|10)$", rownames(seqtab_nochim_AMF)), ]

# Create ASV table
## Cleaning ASV names for FASTA file

asv_fasta <- seqtab2fasta(seqtab_nochim_AMF)

seqtab_nochim_AMF <- t(seqtab_nochim_AMF) # Retaining sequences and asigning shorthand ASV names

row.names(seqtab_nochim_AMF) <- sub(">", "", asv_fasta$asv_headers)

asv_AMF <- otu_table(seqtab_nochim_AMF, taxa_are_rows = TRUE)


## Save FASTA file
dir.create("data/output/processed/sequences/", recursive = TRUE)
write(
  asv_fasta$asv_fasta,
  file.path("data/output/processed/sequences/energy_farm_collab_AMF.fa")
)

# AMF Phyloseq object ------------------------
# Checking that metadata and asvs have the same number of samples

samples_missing_metadata <- base::setdiff(
  colnames(seqtab_nochim_AMF),
  rownames(ef_metadata)
)
# No missing samples

ef_AMF_physeq <- phyloseq(asv_AMF, taxa_AMF, ef_metadata)

otu_table(ef_AMF_physeq)
tax_table(ef_AMF_physeq)
sample_data(ef_AMF_physeq)

# Any other missing or dropped sample?
post_physeq_missing <- base::setdiff(
  sample_names(ef_metadata),
  sample_names(ef_AMF_physeq)
)
missing_sample <- post_physeq_missing %in% colnames(seqtab_nochim_AMF)
# No missing samples

#--------------------------------------------------------
# Save results
#--------------------------------------------------------

# # Save phyloseq object as RDS file
# Name change
ef_physeq_list <- list(
  ef_16S_physeq = ef_16S_physeq,
  ef_AMF_physeq = ef_AMF_physeq
)

#dir.create("data/output/processed/rdata/phyloseq/", recursive = TRUE)
save(
  ef_physeq_list,
  file = "data/output/processed/rdata/phyloseq/ef_physeq_list.rda"
)

# Save phyloseqs independently
purrr::iwalk(
  ef_physeq_list,
  ~ {
    assign(.y, .x)
    save(list = .y, file = file.path("data/output/processed/rdata/phyloseq/", paste0(.y, ".rda")))
  }
)
