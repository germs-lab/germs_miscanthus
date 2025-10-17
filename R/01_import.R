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
    save(
      list = .y,
      file = file.path(
        "data/output/processed/rdata/phyloseq/",
        paste0(.y, ".rda")
      )
    )
  }
)


#------------------------------------------------
# LAMPS: 2018
#------------------------------------------------

# Cleaning up meta data to match up sequence file names
# Reconstruction of file names. See README in ~/data/input/LAMPS/README.md

# DNA amplicons ---------------------
# Proper metadata worlflow
build_proper_metadata <- function(
  csv,
  xlsx,
  .distinct = TRUE,
  rna = FALSE
) {
  meta_csv <- read_csv(
    csv,
    col_types = cols(.default = col_character()),
    show_col_types = TRUE
  ) %>%
    janitor::clean_names() %>% # Everything should be character/categorical
    select(sample_id, sample_date, soil_type) %>%
    as_tibble()

  meta_xls <- read_xlsx(
    xlsx,
    sheet = "useful_metadata",
    .name_repair = "unique",
    col_types = "text"
  ) %>%
    janitor::clean_names() %>%
    rename(sample_id = id) %>%
    mutate(
      plant = case_when(
        plant == "Corn" ~ "C",
        plant == "Miscanthus" ~ "M",
        TRUE ~ plant
      ),
      plot = paste0("P", plot),
      nitrogen_conc = paste0("N", nitrogen_conc)
    ) %>%
    select(!samples)

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
      )
    ) %>%
    tidyr::unite(
      "sequence_id",
      c(sample_id, plant, year, plot, nitrogen_conc, replicate, sample_date),
      sep = "_",
      remove = FALSE,
      na.rm = TRUE
    ) %>%
    relocate(sequence_id, .before = sample_id)

  if (rna) {
    meta <- meta_xls %>%
      left_join(., meta_csv, by = "sample_id") %>%
      relocate(c(year:nitrogen_conc), .after = plant) %>%
      relocate(sample_date, .after = replicate) %>%
      mutate(
        # For re- samples, copy sample_date from the original sample
        sample_date = case_when(
          str_equal(sample_date, "20180429") ~ "20180430", # potential typo in seuqence file naming
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

      mutate(sequence_id = paste(sequence_id, "RNA", sep = "_"))
  }

  if (.distinct) {
    meta <- meta %>%
      distinct(sequence_id, .keep_all = TRUE)
  }

  # Basic checks
  stopifnot(!any(is.na(meta$sample_id)))
  if (anyDuplicated(meta$sample_id)) {
    warning("Duplicate sample_id detected.")
  }
  if (any(is.na(meta$plant))) {
    warning("Missing plant for some rows after join.")
  }

  meta
}
#---------------

csv_path <- "data/input/LAMPS/DNA_amplicon_sequencing/bacterial-metadata.csv"
xlsx_path <- "data/input/LAMPS/DNA_amplicon_sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx"

dna_proper_metadata <- build_proper_metadata(
  csv = csv_path,
  xlsx = xlsx_path,
  .distinct = FALSE
)

# Compare sequence_id with actual sequence file names.
# Actual files
sequence_file_names <- fs::dir_ls(
  "data/input/LAMPS/DNA_amplicon_sequencing/raw_sequences"
)

clean_sequence_file_names <- sequence_file_names %>%
  basename() %>%
  str_remove("_CABBI\\.R[12]\\.fastq$") %>%
  ifelse(
    # Correcting inverted fields in sequence file names.
    str_detect(., "_\\d{8}_BULK_DNA$"), # ends with date_BULK_DNA
    str_replace(., "^(.+)_(\\d{8})_BULK_DNA$", "\\1_BULK_DNA_\\2"),
    . # keep as-is if already correct
  )

# Check matches
clean_sequence_file_names %in%
  {
    dna_proper_metadata %>%
      pull(sequence_id)
  }

setdiff(clean_names, {
  dna_proper_metadata %>%
    pull(sequence_id)
})

# No mismatches. All sequence files have corresponding medata sequence_id.

# The proper metadata was generated to reconstruct the existing sequence file names in raw_sequences/

# Save
write_csv(
  dna_proper_metadata,
  "data/input/LAMPS/DNA_amplicon_sequencing/dna_proper_metadata.csv",
  na = ""
)


# RNA amplicons---------------------
csv_path <- "data/input/LAMPS/DNA_amplicon_sequencing/bacterial-metadata.csv"
xlsx_path <- "data/input/LAMPS/RNA_amplicon_sequencing/0.metadata_CABBI_LAMPS_RNA.xlsx"


rna_proper_metadata <- build_proper_metadata(
  csv = csv_path,
  xlsx = xlsx_path,
  .distinct = FALSE,
  rna = TRUE
)

duplicated_df <- rna_proper_metadata %>%
  slice(rep(1:n(), each = 2)) # Since the initial df is not in duplicate as are the sample _R1 and _R2, then we duplicate to make sure all matches.

rna_sequence_file_names <- fs::dir_ls(
  "data/input/LAMPS/RNA_amplicon_sequencing/raw_sequences"
)

clean_rna_sequence_file_names <- rna_sequence_file_names %>%
  basename() %>%
  str_remove("_CABBI\\_L001\\_R[12]\\_001\\.fastq\\.gz$") %>%
  ifelse(
    # Correcting inverted fields in sequence file names.
    str_detect(., "_(N\\d+)_(P\\d+)_"), # detects N0, N200, N400, etc. before P
    str_replace(., "^(.+)_(N\\d+)_(P\\d+)_(.+)$", "\\1_\\3_\\2_\\4"),
    . # keep as-is if already correct
  )

# Check matches
clean_rna_sequence_file_names %in%
  {
    duplicated_df %>%
      pull(sequence_id)
  }

# There was a lot of going back and forth here. Missing sequence sample replicates cause it to not follow A, B, C replicate order all the time.

setdiff(clean_rna_sequence_file_names, {
  duplicated_df %>%
    pull(sequence_id)
})

# No mismatches after proper field entering in "replicate" in Excel.

# Save
write_csv(
  rna_proper_metadata,
  "data/input/LAMPS/RNA_amplicon_sequencing/rna_proper_metadata.csv",
  na = ""
)
