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


# =============================================================================
# METADATA PROCESSING FOR LAMPS SEQUENCING DATA
# =============================================================================

# Updated metadata worlflow
build_updated_metadata <- function(
  csv,
  xlsx,
  .distinct = TRUE,
  target_region,
  rna = FALSE
) {
  # Read and clean CSV metadata
  meta_csv <- read_csv(
    csv,
    col_types = cols(.default = col_character()),
    show_col_types = TRUE
  ) %>%
    janitor::clean_names() %>%
    select(sample_id, sample_date, soil_type) %>%
    as_tibble()

  # Read and clean Excel metadata
  meta_xls <- read_xlsx(
    xlsx,
    sheet = "useful_metadata",
    .name_repair = "unique",
    col_types = "text"
  ) %>%
    janitor::clean_names() %>%
    rename(sample_id = id) %>%
    mutate(
      # Standardize plant names to tokens
      plant = case_when(
        plant == "Corn" ~ "C",
        plant == "Miscanthus" ~ "M",
        TRUE ~ plant
      ),
      # Add prefixes to plot and nitrogen concentration
      plot = paste0("P", plot),
      nitrogen_conc = paste0("N", nitrogen_conc)
    ) %>%
    select(!samples)

  # Join metadata and create sequence IDs
  meta <- meta_xls %>%
    left_join(., meta_csv, by = "sample_id") %>%
    relocate(c(year:nitrogen_conc), .after = plant) %>%
    relocate(sample_date, .after = replicate) %>%
    mutate(
      # Normalize re-extraction marker
      sample_id = if_else(
        str_starts(sample_id, "Reextracted-"),
        str_replace(sample_id, "^Reextracted-", "re-"),
        sample_id
      ),
      # For re- samples, copy sample_date from the original sample
      sample_date = case_when(
        str_starts(sample_id, "re-") & is.na(sample_date) ~
          {
            # Extract the number after "re-"
            original_id <- str_extract(sample_id, "(?<=re-)\\d+")
            # Find the sample_date from the matching original sample
            sample_date[match(original_id, sample_id)]
          },
        TRUE ~ sample_date
      ),
      target_region = target_region
    ) %>%
    relocate(target_region, .after = nucleotide) %>%
    tidyr::unite(
      "sequence_id",
      c(sample_id, plant, year, plot, nitrogen_conc, replicate, sample_date),
      sep = "_",
      remove = FALSE,
      na.rm = TRUE
    ) %>%
    relocate(sequence_id, .before = sample_id)

  # RNA-specific processing
  if (rna) {
    meta <- meta %>%
      mutate(
        # Fix potential typo in sequence file naming
        sample_date = case_when(
          str_equal(sample_date, "20180429") ~ "20180430",
          TRUE ~ sample_date
        )
      ) %>%
      tidyr::unite(
        "sequence_id",
        c(
          sample_id,
          plant,
          year,
          plot,
          nitrogen_conc,
          replicate,
          sample_date
        ),
        sep = "_",
        remove = FALSE,
        na.rm = TRUE
      ) %>%
      relocate(sequence_id, .before = sample_id) %>%
      mutate(
        # Add RNA suffix to sequence_id
        sequence_id = paste(sequence_id, "RNA", sep = "_")
      )
  }

  # Remove duplicates if requested
  if (.distinct) {
    meta <- meta %>%
      distinct(sequence_id, .keep_all = TRUE)
  }

  # Data validation checks
  stopifnot(!any(is.na(meta$sample_id)))

  if (anyDuplicated(meta$sample_id)) {
    warning("Duplicate sample_id detected.")
  }

  if (any(is.na(meta$plant))) {
    warning("Missing plant for some rows after join.")
  }

  return(meta)
}


regen_physeq <- function(physeq, sample_metadata, rownames = "sample_id") {
  phyloseq(
    otu_table(physeq),
    tax_table(physeq),
    sample_data(sample_metadata %>% column_to_rownames(., var = rownames))
  )
}


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
  target_region = "16S"
)

# Get actual sequence file names
dna_seq_names <- fs::dir_ls(
  "data/input/LAMPS/DNA_amplicon_sequencing/raw_sequences"
)

# Clean sequence file names and correct inverted fields
dna_seq_names_clean <- dna_seq_names %>%
  basename() %>%
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
write_csv(
  dna_updated_metadata %>%
    distinct(sequence_id, .keep_all = TRUE),
  "data/input/LAMPS/DNA_amplicon_sequencing/dna_updated_metadata.csv",
  na = ""
)

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
  rna = TRUE
)

# Duplicate metadata rows to match R1/R2 paired-end files
rna_duplicated_df <- rna_updated_metadata %>%
  slice(rep(1:n(), each = 2))

# Get actual RNA sequence file names
rna_seq_names <- fs::dir_ls(
  "data/input/LAMPS/RNA_amplicon_sequencing/raw_sequences"
)

# Clean RNA sequence file names and correct inverted fields
rna_seq_names_clean <- rna_seq_names %>%
  basename() %>%
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
write_csv(
  rna_updated_metadata,
  "data/input/LAMPS/RNA_amplicon_sequencing/rna_updated_metadata.csv",
  na = ""
)

# --------------------------------------------------------------------------
# Metadata update to LAMPS 2018 phyloseq object ----------------------------
# --------------------------------------------------------------------------

ps_16S_LAMPS <- readRDS("data/input/LAMPS/ps_16S_LAMP.rds")

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


# Updated phyloseq

lamps_2018_16S_physeq <- regen_physeq(
  ps_16S_LAMPS,
  complete_16S_metadata,
  rownames = "rownames"
)

#head(sample_data(lamps_2018_physeq))

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
    target_region = "ITS"
  ) %>%
  select(!samples) %>%
  relocate(sequence_id, .before = sample_id) %>%
  relocate(plant, .before = year) %>%
  relocate(plot, .before = nitrogen_conc) %>%
  relocate(target_region, .after = nucleotide)


# Sequence file names
its_seq_names <- fs::dir_ls(
  "data/input/LAMPS/ITS_sequencing/raw_sequences"
)

# Validate matches between cleaned file names and metadata
its_seq_names_clean <- its_seq_names %>%
  basename() %>%
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

#----------------------------------------------------------------------------
# LAMPS: 2023
#----------------------------------------------------------------------------

# Phyloseq objects
# 16S
lamps_16S_2022 <- readRDS("data/input/LAMPS/2022/data/LAMPS_EPS_16S.rds")

# AMF
lamps_AMF_2022 <- readRDS("data/input/LAMPS/2022/data/LAMPS_EPS_AMF.rds")

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


# Let's reconstruct the phyloseq object with simpler metadata

nphyseq_lamps_16S_2022 <- regen_physeq(
  lamps_16S_2022,
  sample_metadata = lamps_metadata_2022,
)

sample_names(nphyseq_lamps_16S_2022)


nphyseq_lamps_AMF_2022 <- regen_physeq(
  lamps_AMF_2022,
  sample_metadata = lamps_metadata_2022
)

sample_names(nphyseq_lamps_AMF_2022)

lamps_2022_physeq_list <- list(
  lamps_16S_2022 = nphyseq_lamps_16S_2022,
  lamps_AMF_2022 = nphyseq_lamps_AMF_2022
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
