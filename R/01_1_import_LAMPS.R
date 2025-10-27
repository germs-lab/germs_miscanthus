###################################################################
# Import of R objects, .csv and .xlsx for Miscanthus overview
#
# Import sequence, phyloseqs and metadata files from:
# - LAMPS: 2018
#   - 16S (DNA & RNA)
#   - ITS (DNA & RNA)
# - LAMPS: 2022
#   - 16S (DNA & RNA)
#   - ITS (DNA & RNA)
#
# For LAMPS metadata: Clean and standardize metadata to match
# sequence file  names. See README in ~/data/input/LAMPS/README.md
#
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
dna_seq_names <- readr::read_lines(
  "data/input/LAMPS/DNA_amplicon_sequencing/raw_sequences/dna_seq_filenames.tsv"
) %>%
  stringr::str_trim() %>%
  purrr::discard(~ .x == "dna_seq_filenames.tsv")

# Clean sequence file names and correct inverted fields
dna_seq_names_clean <- dna_seq_names %>%
  str_remove("_CABBI\\.R[12]\\.fastq$") %>%
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
rna_seq_names <- readr::read_lines(
  "data/input/LAMPS/RNA_amplicon_sequencing/raw_sequences/rna_seq_filenames.tsv"
) %>%
  stringr::str_trim() %>%
  purrr::discard(~ .x == "rna_seq_filenames.tsv")


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
lamps_2018_16S_tax <- tax_table(ps_16S_LAMPS) %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sequence") %>%
  rename_with(str_to_lower, .cols = everything()) # Clean up needed after importing from .csv
rownames(lamps_2018_16S_tax) <- lamps_2018_16S_tax %>%
  pull(sequence)


base::setdiff(
  rownames(otu_table(ps_16S_LAMPS)),
  rownames(lamps_2018_16S_tax)
)

# Updated phyloseq
# Check class for replacement
class(tax_table(ps_16S_LAMPS))
class(tax_table(as.matrix(lamps_2018_16S_tax)))

tax_table(ps_16S_LAMPS) <- tax_table(as.matrix(lamps_2018_16S_tax))

lamps_2018_16S_physeq <- regen_physeq(
  ps_16S_LAMPS,
  complete_16S_metadata,
  rownames = "rownames"
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
    project = "LAMPS_2018",
    rownames = sample_id
  ) %>%
  select(!samples) %>%
  relocate(sequence_id, .before = sample_id) %>%
  relocate(plant, .before = year) %>%
  relocate(plot, .before = nitrogen_conc) %>%
  relocate(target_region, .after = nucleotide)


# Sequence file names
its_seq_names <- readr::read_lines(
  "data/input/LAMPS/ITS_sequencing/raw_sequences/its_seq_filenames.tsv"
) %>%
  stringr::str_trim() %>%
  purrr::discard(~ .x == "its_seq_filenames.tsv")


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

# Matching metadata to Millican created phyloseq
load("data/input/LAMPS/ITS_sequencing/manuscript_all_data.rda")
rm(df.meta)

# Mismatches
base::setdiff(
  {
    unique(its_seq_names_clean %>% str_remove("_[^_]+$"))
  },
  sample_names(ps.f)
)


# Remaking the phyloseq object like we want to
ps.f_otu_transposed <- phyloseq::phyloseq(
  tax_table(ps.f),
  t(as.matrix(otu_table(ps.f, taxa_are_rows = TRUE)))
)
# Clean up to taxonomy table
lamps_2018_ITS_tax <- tax_table(ps.f_otu_transposed) %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sequence") %>%
  rename_with(str_to_lower, .cols = everything())
rownames(lamps_2018_ITS_tax) <- lamps_2018_ITS_tax %>%
  pull(sequence)


tax_table(ps.f_otu_transposed) <- tax_table(as.matrix(lamps_2018_ITS_tax))

# Meta clean
its_meta_distinct <- its_metadata %>%
  distinct(sample_id, .keep_all = TRUE)


lamps_2018_ITS_physeq <- regen_physeq(
  ps.f_otu_transposed,
  its_meta_distinct,
  rownames = "rownames"
)

# Saving
lamps_2018_physeq_list <- list(
  lamps_2018_16S_physeq = lamps_2018_16S_physeq,
  lamps_2018_ITS_physeq = lamps_2018_ITS_physeq
)


save(
  lamps_2018_physeq_list,
  file = "data/output/processed/rdata/phyloseq/lamps_2018_physeq_list.rda"
)

# Save phyloseqs independently
purrr::iwalk(
  lamps_2018_physeq_list,
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
    sequence_id = sample_id,
    rownames = sample_id
  ) %>% # Sampling year
  relocate(sequence_id, .before = sample_id)

# Fix taxa orientation and names
otu_table(lamps_2022_16S) <- otu_table(t(otu_table(lamps_2022_16S)))
otu_table(lamps_2022_AMF) <- otu_table(t(otu_table(lamps_2022_AMF)))

# Update taxa names
list(lamps_2022_16S, lamps_2022_AMF) %>%
  purrr::map(
    ~ {
      taxa_names(.x) <- gsub("^ASV([0-9]+)$", "ASV_\\1", taxa_names(.x))
      .x
    }
  ) %>%
  set_names(c("lamps_2022_16S", "lamps_2022_AMF")) %>%
  list2env(envir = .GlobalEnv)


# Let's reconstruct the phyloseq object with simpler metadata
nphyseq_lamps_2022_16S <- regen_physeq(
  lamps_2022_16S,
  sample_metadata = lamps_metadata_2022,
  rownames = "rownames"
)

sample_names(nphyseq_lamps_2022_16S)


nphyseq_lamps_2022_AMF <- regen_physeq(
  lamps_2022_AMF,
  sample_metadata = lamps_metadata_2022,
  rownames = "rownames"
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
