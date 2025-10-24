###################################################################
# Import of R objects, .csv and .xlsx for Miscanthus overview
#
# Import sequence, phyloseqs and metadata files from:
# - Energy Farm Collab project
# - LAMPS: 2018
#   - 16S
#     - DNA
#     - RNA
# For LAMPS metadata: Clean and standardize metadata to match sequence file  names. See README in ~/data/input/LAMPS/README.md

# Author: Bolívar Aponte Rolón
# Date: 2025-10-01
###################################################################

source("R/utils/00_setup.R")

# =============================================================================
# METADATA PROCESSING FOR LAMPS SEQUENCING DATA
# =============================================================================
#----------------------------------------------------------------------------
# LAMPS: 2018
#----------------------------------------------------------------------------
# 16S -----------------------------------------------------------------------

# DNA amplicon processing ---------------------------------------------------

# Define file paths
csv_path <- "data/input/LAMPS/DNA_amplicon_sequencing/bacterial-metadata.csv"
xlsx_path <- "data/input/LAMPS/DNA_amplicon_sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx"

# Build DNA metadata
dna_updated_metadata <- build_updated_metadata(
  csv = csv_path,
  xlsx = xlsx_path,
  .distinct = FALSE,
  target_region = "16S",
  .project = "LAMPS_2018"
)

# Get actual sequence file names
dna_seq_names <- readr::read_tsv(
  "data/input/LAMPS/DNA_amplicon_sequencing/dna_seq_filenames.tsv",
  col_types = readr::cols(.default = readr::col_character())
) %>%
  dplyr::pull(1)


# Clean sequence file names and correct inverted fields
dna_seq_names_clean <- dna_seq_names %>%
  str_remove("_CABBI\\.R[12]\\.fastq\\.gz$") %>%
  # Correct BULK_DNA position: move from after date to before date
  ifelse(
    str_detect(., "_\\d{8}_BULK_DNA$"), # ends with date_BULK_DNA
    str_replace(., "^(.+)_(\\d{8})_BULK_DNA$", "\\1_BULK_DNA_\\2"),
    . # keep as-is if already correct
  )

# Validate matches between cleaned file names and metadata
matches <- dna_seq_names_clean %in%
  {
    dna_updated_metadata %>% pull(sequence_id)
  }

# Check for mismatches
dna_mismatches <- setdiff(dna_seq_names_clean, {
  dna_updated_metadata %>% pull(sequence_id)
})

# No mismatches expected, all sequence files should have corresponding metadata
if (length(dna_mismatches) > 0) {
  warning(
    "Mismatched sequence files found: ",
    paste(dna_mismatches, collapse = ", ")
  )
}

# Save DNA metadata
# write_csv(
#   dna_updated_metadata %>%
#     distinct(sequence_id, .keep_all = TRUE),
#   "data/input/LAMPS/DNA_amplicon_sequencing/dna_updated_metadata.csv",
#   na = ""
# )

# RNA amplicon processing ---------------------------------------------------

# Define file paths (note: uses same CSV as DNA)
csv_path <- "data/input/LAMPS/DNA_amplicon_sequencing/bacterial-metadata.csv"
xlsx_path <- "data/input/LAMPS/RNA_amplicon_sequencing/0.metadata_CABBI_LAMPS_RNA.xlsx"

# Build RNA metadata
rna_updated_metadata <- build_updated_metadata(
  csv = csv_path,
  xlsx = xlsx_path,
  .distinct = FALSE,
  target_region = "16S",
  rna = TRUE,
  .project = "LAMPS_2018"
)

# Duplicate metadata rows to match R1/R2 paired-end files
rna_duplicated_df <- rna_updated_metadata %>%
  slice(rep(1:n(), each = 2))

# Get actual RNA sequence file names
rna_seq_names <- readr::read_tsv(
  "data/input/LAMPS/RNA_amplicon_sequencing/rna_seq_filenames.tsv",
  col_types = readr::cols(.default = readr::col_character())
) %>%
  dplyr::pull(1)


# Clean RNA sequence file names and correct inverted fields
rna_seq_names_clean <- rna_seq_names %>%
  str_remove("_CABBI\\_L001\\_R[12]\\_001\\.fastq\\.gz$") %>%
  # Correct N/P field order: swap nitrogen (N) and plot (P) positions
  ifelse(
    str_detect(., "_(N\\d+)_(P\\d+)_"), # detects N0, N200, N400, etc. before P
    str_replace(., "^(.+)_(N\\d+)_(P\\d+)_(.+)$", "\\1_\\3_\\2_\\4"),
    . # keep as-is if already correct (P before N)
  )


# Validate matches between cleaned file names and metadata
rna_matches <- rna_seq_names_clean %in%
  {
    rna_duplicated_df %>% pull(sequence_id)
  }

# Check for mismatches
rna_mismatches <- setdiff(rna_seq_names_clean, {
  rna_duplicated_df %>% pull(sequence_id)
})


if (length(rna_mismatches) > 0) {
  warning(
    "Mismatched RNA sequence files found: ",
    paste(rna_mismatches, collapse = ", ")
  )
}
# Note: Missing sequence sample replicates cause irregular A, B, C replicate ordering, thus mismatches. This was resolved by updated field entry in "replicate" column in Excel.

# Save RNA metadata
# write_csv(
#   rna_updated_metadata,
#   "data/input/LAMPS/RNA_amplicon_sequencing/rna_updated_metadata.csv",
#   na = ""
# )

# --------------------------------------------------------------------------
# Metadata update to LAMPS 2018 phyloseq object ----------------------------
# --------------------------------------------------------------------------

ps_16S_LAMPS <- readRDS(
  "data/input/LAMPS/16S_rRNA_phyloseq_data/ps_16S_LAMP.rds"
)

#head(sample_data(ps_16S_LAMPS))

# One last little clean up
dna_updated_metadata <- dna_updated_metadata %>%
  distinct(sequence_id, .keep_all = TRUE) %>%
  mutate(
    sample_id = paste(sample_id, nucleotide, sep = "_"),
    rownames = sample_id
  )

rna_updated_metadata <- rna_updated_metadata %>%
  mutate(
    sample_id = paste(sample_id, nucleotide, sep = "_"),
    rownames = sample_id
  )

complete_16S_metadata <- rbind(dna_updated_metadata, rna_updated_metadata)

# Clean up to taxonomy table
lamps_2018_tax <- tax_table(ps_16S_LAMPS) %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sequence") %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv
rownames(lamps_2018_tax) <- lamps_2018_tax %>%
  pull(sequence)


base::setdiff(
  rownames(otu_table(ps_16S_LAMPS)),
  rownames(lamps_2018_tax)
)

# Updated phyloseq
# Check class for replacement
class(tax_table(ps_16S_LAMPS))
class(tax_table(as.matrix(lamps_2018_tax)))

tax_table(ps_16S_LAMPS) <- tax_table(as.matrix(lamps_2018_tax))

lamps_2018_16S_physeq <- regen_physeq(
  ps_16S_LAMPS,
  complete_16S_metadata,
  rownames = "rownames"
)

# Add refseq
lamps_2018_16S_physeq <- add_refseq(lamps_2018_16S_physeq, tag = "ASV_")

save(
  lamps_2018_16S_physeq,
  file = "data/output/processed/rdata/phyloseq/lamps_2018_16S_physeq.rda"
)

# ITS -----------------------------------------------------------------------

# Import metadata

its_metadata <- read_xlsx(
  "data/input/LAMPS/ITS_sequencing/0.metadata.xlsx",
  sheet = "Sample_name",
  col_types = "text"
) %>%
  janitor::clean_names() %>%
  rename(sample_id = name_of_sequence) %>%
  mutate(
    sequence_id = sample_id,
    sample_id = str_remove(sample_id, "_(.+)$"),
    plant = case_when(
      plant == "Corn" ~ "C",
      plant == "Miscanthus" ~ "M",
      TRUE ~ plant
    ),
    nitrogen_conc = paste0("N", nitrogen_conc),
    target_region = "ITS",
    project = "LAMPS_2018"
  ) %>%
  select(!samples) %>%
  relocate(sequence_id, .before = sample_id) %>%
  relocate(plant, .before = year) %>%
  relocate(plot, .before = nitrogen_conc) %>%
  relocate(target_region, .after = nucleotide)


# Sequence file names
its_seq_names <- readr::read_tsv(
  "data/input/LAMPS/ITS_sequencing/its_seq_filenames.tsv",
  col_types = readr::cols(.default = readr::col_character())
) %>%
  dplyr::pull(1)

# Validate matches between cleaned file names and metadata
its_seq_names_clean <- its_seq_names %>%
  #basename() %>%
  str_remove("_S(.+)$")


its_mismatches <- setdiff(its_seq_names_clean, {
  its_metadata %>% pull(sequence_id)
})

# No mismatches

write_csv(
  its_metadata,
  "data/input/LAMPS/ITS_sequencing/its_updated_metadata.csv",
  na = ""
)
#TODO
load("data/input/LAMPS/ITS_sequencing/manuscript_all_data.rda")

intersect(sample_names(ps.b), {
  its_seq_names_clean %>% str_remove("_(.+)$")
})


# I think that Millican subseted the 16S data to only the samples that matched the ITS sampling effort.

name_test <- its_metadata %>%
  pull(sample_id) %>%
  str_remove("_RNA")

test <- prune_samples(sample_names(lamps_2018_16S_physeq), sample_names(ps.b))
#----------------------------------------------------------------------------
# LAMPS: 2022
#----------------------------------------------------------------------------

# Phyloseq objects
# 16S
lamps_2022_16S <- readRDS("data/input/LAMPS/2022/data/LAMPS_EPS_16S.rds")

# AMF
lamps_2022_AMF <- readRDS("data/input/LAMPS/2022/data/LAMPS_EPS_AMF.rds")

# Metadata
lamps_metadata_2022 <- read_xlsx(
  "data/input/LAMPS/2022/data/LAMPS_EPS_metadata.xlsx",
  sheet = "sample_data"
) %>%
  janitor::clean_names() %>%
  select(sample_id:replicate) %>%
  mutate(
    replicate = as.character(replicate),
    year = "2022",
    sequence_id = sample_id
  ) %>% # Sampling year
  relocate(sequence_id, .before = sample_id)


# Fix taxa names to be readable
taxa_names(lamps_2022_16S) <- paste0("ASV_", seq(ntaxa(lamps_2022_16S)))
taxa_names(lamps_2022_AMF) <- paste0("ASV_", seq(ntaxa(lamps_2022_AMF)))

# Let's reconstruct the phyloseq object with simpler metadata

nphyseq_lamps_2022_16S <- regen_physeq(
  lamps_2022_16S,
  sample_metadata = lamps_metadata_2022,
)

sample_names(nphyseq_lamps_2022_16S)


nphyseq_lamps_2022_AMF <- regen_physeq(
  lamps_2022_AMF,
  sample_metadata = lamps_metadata_2022
)

sample_names(nphyseq_lamps_2022_AMF)


lamps_2022_physeq_list <- list(
  lamps_2022_16S_physeq = nphyseq_lamps_2022_16S,
  lamps_2022_AMF_physeq = nphyseq_lamps_2022_AMF
)


save(
  lamps_2022_physeq_list,
  file = "data/output/processed/rdata/phyloseq/lamps_2022_physeq_list.rda"
)


# Save phyloseqs independently
purrr::iwalk(
  lamps_2022_physeq_list,
  ~ {
    assign(.y, .x)
    save(
      list = .y,
      file = file.path(
        "data/output/processed/rdata/phyloseq/",
        paste0(.y, ".rda")
      )
    )
  }
)


# =============================================================================
# Files from Energy Farm Collab
# =============================================================================

# Metadata ---------------------------------------------------
# For 16S and AMF
# "data/input/energy_farm_collab/files_for_phyloseq_16S/ef2024_sampledata.csv" and
# "data/input/energy_farm_collab/files_for_phyloseq_AMF/ef2024_sampledata.csv" are the same

ef_metadata <- read.csv(
  "data/input/energy_farm_collab/files_for_phyloseq_16S/ef2024_sampledata.csv"
) %>%
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
rownames(taxa_16S) <- taxa_16S %>%
  pull(sequence)

taxa_16S <- tax_table(as.matrix(taxa_16S))

# ASVs ------------------------
seqtab_nochim_16S <- readRDS(
  "data/input/energy_farm_collab/files_for_phyloseq_16S/seqtab.nochim.rds"
)
# Create ASV table

seqtab_nochim_16S <- t(seqtab_nochim_16S) # Retaining sequences and asigning shorthand ASV names

asv_16S <- otu_table(seqtab_nochim_16S, taxa_are_rows = TRUE)

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

# Add refseq
ef_16S_physeq <- add_refseq(ef_16S_physeq, tag = "ASV_")

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
  "data/input/energy_farm_collab/files_for_phyloseq_AMF/taxaf.rds"
) %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sequence") %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv

rownames(taxa_AMF) <- taxa_AMF %>%
  pull(sequence)

taxa_AMF <- tax_table(as.matrix(taxa_AMF))

# ASVs ------------------------
seqtab_nochim_AMF <- readRDS(
  "data/input/energy_farm_collab/files_for_phyloseq_AMF/seqtabf.nochim.rds"
) # Forward reads only

# seqtab_nochim_AMF needs a sample name clean up to match metadata
rownames(seqtab_nochim_AMF) <- gsub("_S.*", "", rownames(seqtab_nochim_AMF))

# Remove "method" samples
seqtab_nochim_AMF <- seqtab_nochim_AMF[
  !grepl("^method([1-9]|10)$", rownames(seqtab_nochim_AMF)),
]

# Create ASV table
## Cleaning ASV names for FASTA file
seqtab_nochim_AMF <- t(seqtab_nochim_AMF) # Retaining sequences and asigning shorthand ASV names

asv_AMF <- otu_table(seqtab_nochim_AMF, taxa_are_rows = TRUE)


# AMF Phyloseq object ------------------------
# Checking that metadata and asvs have the same number of samples

samples_missing_metadata <- base::setdiff(
  colnames(seqtab_nochim_AMF),
  rownames(ef_metadata)
)

non_matching_asvs <- base::setdiff(
  rownames(asv_AMF),
  rownames(taxa_AMF)
)
# No missing samples

ef_AMF_physeq <- phyloseq(asv_AMF, taxa_AMF, ef_metadata)

otu_table(ef_AMF_physeq)
tax_table(ef_AMF_physeq)
sample_data(ef_AMF_physeq)

# Add refseq
ef_AMF_physeq <- add_refseq(ef_AMF_physeq, tag = "ASV_")


# Any other missing or dropped sample?
post_physeq_missing <- base::setdiff(
  sample_names(ef_metadata),
  sample_names(ef_AMF_physeq)
)
missing_sample <- post_physeq_missing %in% colnames(seqtab_nochim_AMF)
# No missing samples

# Save results ---------------------------------------------------

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
    save(
      list = .y,
      file = file.path(
        "data/output/processed/rdata/phyloseq/",
        paste0(.y, ".rda")
      )
    )
  }
)
