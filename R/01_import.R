###################################################################
# Import of R object from DADA2 pipeline for phyloseq generation
#
#
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-01
###################################################################
###################################################################

source("R/utils/00_setup.R")

#-----------------------------------------------
# Files from Energy Farm Collab
#-----------------------------------------------

# Metadata ------------------------
# For 16S and AMF
# "data/input/energy_farm_collab/files_for_phyloseq_16S/ef2024_sampledata.csv" and
# "data/input/energy_farm_collab/files_for_phyloseq_AMF/ef2024_sampledata.csv" are the same

ef_metadata <- read.csv(
  "data/input/energy_farm_collab/files_for_phyloseq_16S/ef2024_sampledata.csv"
) %>%
  janitor::clean_names(.) %>%
  column_to_rownames(var = "label_id")

metadata_16S <- sample_data(ef_metadata)

#================
# 16S
#================
# Taxonomy ------------------------
taxa <- readRDS(
  "data/input/energy_farm_collab/files_for_phyloseq_16S/taxa.rds"
) %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sequence") %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv

rownames(taxa) <- paste0("ASV_", 1:nrow(taxa))

tax_16S <- tax_table(as.matrix(taxa))

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
  rownames(ef_16S_metadata)
)

ef_16S_physeq <- phyloseq(asv_16S, tax_16S, metadata_16S)

# Any other missing or dropped sample?
post_physeq_missing <- base::setdiff(
  sample_names(metadata_16S),
  sample_names(ef_16S_physeq)
)
missing_sample <- post_physeq_missing %in% colnames(seqtab_nochim_16S)


# Phyloseq

#================
# AMF (only forward reads)
#================
# Metadata
# data/input/energy_farm_collab/files_for_phyloseq_16S/ef2024_sampledata.csv
#
# Taxonomy
#
# ASVs
#
# Phyloseq

#--------------------------------------------------------
# Step 12: Create phyloseq object
#--------------------------------------------------------

#--------------------------------------------------------
# Step 13: Add sample metadata
#--------------------------------------------------------

# Load metadata from CSV file
metadata_1 <- readxl::read_xlsx(
  "data/input/metadata_CABBI_SABR2023_DNA.xlsx"
) |>
  janitor::clean_names() |>
  mutate(
    id = gsub(
      pattern = "^(.*?)_R\\d+.*$",
      replacement = "\\1",
      x = id
    ),
    across(where(is.numeric), as.factor)
  ) |>
  distinct(id, .keep_all = TRUE) |>
  rename(sample_id = id)
# column_to_rownames(var = "id")

metadata_2 <- readxl::read_xlsx(
  "data/input/overall_data.xlsx",
  sheet = "OVERALL",
) |>
  janitor::clean_names() |>
  mutate(
    initial_order = as.character(initial_order),
    sampling = as.character(sampling),
    plot = as.character(plot),
    column = as.character(column),
    plate_number = as.character(plate_number),
    sample = sprintf('SABR_%s', sample),
    across(dna_conc_ng_ul:gwc_g_g, as.numeric)
  ) |>
  rename(sample_id = sample) |>
  select(c(sample_id, dna_conc_ng_ul:gwc_g_g))

metadata_join <- dplyr::left_join(metadata_1, metadata_2, by = "sample_id") |>
  column_to_rownames(var = "sample_id")

sabr_2023_metadata_clean <- metadata_join # Name change for saving purposes
save(
  sabr_2023_metadata_clean,
  file = "data/output/processed/metadata/sabr_2023_metadata_clean.rda"
)

# Check for sample name consistency between phyloseq and metadata
head(sample_names(physeq))
head(rownames(metadata_join))
base::setdiff(sample_names(physeq), rownames(metadata_join)) # In phyloseq but not in metadata
base::setdiff(rownames(metadata_join), sample_names(physeq)) # In metadata but not in phyloseq

# Convert metadata to phyloseq sample_data format
sampledata <- sample_data(metadata_join)
sampledata

# Add metadata to phyloseq object
physeq <- merge_phyloseq(physeq, sampledata)
physeq

# Verify phyloseq components
otu_table(physeq) # ASV table
sample_data(physeq) # Sample metadata
tax_table(physeq) # Taxonomy table

# #--------------------------------------------------------
# # Step 14: Save results
# #--------------------------------------------------------

# # Save phyloseq object as RDS file
# Name change
sabr_2023_physeq <- physeq # Name change for saving purposes
remove(physeq)
save(
  sabr_2023_physeq,
  file = "data/output/processed/rdata/phyloseq/sabr_2023_physeq.rda"
)
