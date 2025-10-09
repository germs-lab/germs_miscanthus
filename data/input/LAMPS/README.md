# Building proper metadata in R from two sources

This document explains how the final, “proper” metadata table was built entirely in R by combining two required inputs:

- Lab analysis from Micheal Millican
  - GERMS-DATAMAN/Millican_GERMS_Research/LAMPS_microbiome/analysis_data/bacterial-metada.csv
- Sequencing master metadata (Box)
  - GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA amplicon sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx
 or 
 DOE-CABBI/bolivar/germs_miscanthus/data/input/LAMPS/DNA_amplicon_sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx

Both files were necessary, although they had more or less the same data: the CSV provided the lab-side identifiers and treatments; the Excel file provided canonical sample annotations (e.g., plant type) needed to standardize names and downstream mappings.

This applies to DNA and RNA sequences.

## Why two files?
The sequence file names in `~Dna amplicon sequencing/raw_sequences` as the results of the sequencing run. They were probably built using `bacterial-metada.csv` where it had sampling date as one of the fields but the Excel file did not. I am using both as a source to reconstruct the sequence names and properly annotate the samples.
The two files contain overlapping but not identical information:
- bacterial-metada.csv: holds experiment- and lab-level fields (sample IDs, plot, N-rate, etc.).
- 0.metadata_CABBI_LAMPS_DNA.xlsx: holds sequencing metadata and canonical attributes (e.g., plant species) to normalize and disambiguate samples.

## R workflow (overview)

1) Read both files (CSV and Excel).
2) Trim/clean join keys (e.g., `sample_id`), align column names as needed.
3) Left-join Excel columns into the CSV by the shared key.
4) Normalize fields (e.g., map plant to tokens C/M; handle “Reextracted-” prefix).
5) Validate: check uniqueness, missing values, row counts.
6) Write the combined metadata for reuse.

## Minimal R example

```r
library(readr)
library(readxl)
library(dplyr)
library(stringr)

csv_path  <- "GERMS-DATAMAN/Millican_GERMS_Research/LAMPS_microbiome/analysis_data/bacterial-metada.csv"
xlsx_path <- "GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA amplicon sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx"

# Adjust the column names below to match your files
build_proper_metadata <- function(csv = csv_path, xlsx = xlsx_path) {
  meta_csv <- read_csv(csv, show_col_types = FALSE) %>%
    rename(sample_id = sample_id) %>%           # ensure consistent naming
    mutate(
      sample_id = str_trim(sample_id),
      plot      = str_trim(plot),
      n_rate    = str_trim(n_rate)
    )

  meta_xls <- read_xlsx(xlsx, .name_repair = "unique") %>%
    rename(sample_id = sample_id, plant = plant) %>%
    mutate(
      sample_id = str_trim(sample_id),
      plant     = str_trim(plant)
    ) %>%
    select(sample_id, plant)

  meta <- meta_csv %>%
    left_join(meta_xls, by = "sample_id") %>%
    mutate(
      # Normalize re-extraction marker
      sample_prefix = if_else(str_starts(sample_id, "Reextracted-"),
                              str_replace(sample_id, "^Reextracted-", "re-"),
                              sample_id),
      # Map plant to tokens used downstream
      plant_token = case_when(
        str_to_lower(plant) == "corn" ~ "C",
        str_to_lower(plant) == "miscanthus" ~ "M",
        TRUE ~ plant
      )
    )

  # Basic checks
  stopifnot(!any(is.na(meta$sample_id)))
  if (anyDuplicated(meta$sample_id)) warning("Duplicate sample_id detected.")
  if (any(is.na(meta$plant))) warning("Missing plant for some rows after join.")
  meta
}

proper_metadata <- build_proper_metadata()

# Optionally, save for reuse
# write_csv(proper_metadata, "analysis_data/proper_metadata.csv")
```

Tip: The Excel path contains spaces. In shell commands, quote it or escape spaces:
- "GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA amplicon sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx"
- GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA\ amplicon\ sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx

## Outputs

- A single, cleaned metadata table (`proper_metadata`) suitable for downstream R analyses and consistent filename construction.
- If file naming is needed, add a column that encodes your naming schema using the normalized fields (e.g., `sample_prefix`, `plant_token`, `plot`, `n_rate`).

## Reproducibility notes

- Keep the two source files versioned (CSV) or archived (Excel from Box).
- Record any column renames made during import so the join remains stable if schemas change.