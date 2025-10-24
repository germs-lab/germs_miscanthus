# Data Provenance and Preparation for "germs_miscanthus"


# LAMPS 2018
## Updating metadata for 16S (DNA and RNA) and ITS (DNA and RNA) samples

This document explains how the final, “updated” metadata table for LAMPS 2018 data was built by combining two required inputs:

Both files were necessary, although they had more or less the same data. Metadata and sequence files could very well have been matched using sample_id(id) since each is unique, but "replicate" information was important as well as "sample_date" which the RNA metadata was missing.

This applies to DNA (DNA amplicon sequencing) and RNA (RNA amplicon sequencing)sequences.


## Data and Metadata provenance

### 16S (DNA and RNA)

- **Metadata**:
  - Lab analysis from Micheal Millican: "GERMS-DATAMAN/Millican_GERMS_Research/LAMPS_microbiome/analysis_data/bacterial-metadata.csv"
  - Sequencing master metadata (Box): "GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA amplicon sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx"

- **Raw sequences**: 
  - "GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA\ amplicon\ sequencing/Raw\ sequences"
  - "GERMS-DATAMAN/DOE-CABBI/LAMPS/RNA\ amplicon\ sequencing/Raw\ sequences"
  - "GERMS-DATAMAN/Millican_GERMS_Research/LAMPS_microbiome/raw_data/"

**To avoid downloading samples, a list of all the samples was imported to "~/germs_miscanthus/data/input/LAMPS" into their respective subdirectories "DNA_*" and "RNA_*"from where it was used for sample corroboration and analyses.**

```{r}
its_seq_names <- fs::dir_ls(
  "/home/baponte/cybox/GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA\ amplicon\ sequencing/Raw sequences"
) %>% # rclone Box config, this may break
  basename() 

write_tsv(
  as.data.frame(its_seq_names),
file = "/home/baponte/cybox/GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA\ amplicon\ sequencing/dna_seq_filenames.tsv", col_names = FALSE
)
```

- **.Rdata (phyloseq objects)**: "GERMS-DATAMAN/DOE-CABBI/LAMPS/16S\ rRNA\ phyloseq\ data/ps_16S_LAMPS.rds"
**NOTE**: The metadata for this phyloseq object was updated to contain sequence_id, target_region and some other clean ups to increase readability for this project. **Now imported to "~germs_miscanthus/data/input/LAMPS/16S_rRNA_phyloseq_data"**

### ITS (DNA and RNA)

This is an interesting import. The files were previously worked on by M. Millican. Now we are going to use the phyloseq object he created. We have 152 sequence samples with 127 being unique. In his phyloseq we end up with 124: "2032" "2040" "5108" as removed. 
 
  - **Metadata**: "GERMS-DATAMAN/DOE-CABBI/LAMPS/ITS\ sequencing/0.metadata.xlsx"
  - **Raw sequences**: "GERMS-DATAMAN/DOE-CABBI/LAMPS/ITS\ sequencing/Raw\ sequences"

  **To avoid downloading samples, a list of all the samples was imported to "~/germs_miscanthus/data/input/LAMPS/ITS_sequencing/raw_sequences" from where it was used for sample corroboration and analyses.**

  - **.Rdata (phyloseq objects)**: "GERMS-DATAMAN/Millican_GERMS_Research/LAMPS_microbiome/manuscript_all_data.Rdata" ("ps.f")



## Why two files?
The sequence file names in "~Dna amplicon sequencing/raw_sequences" are the results of the sequencing run. They were probably built using "bacterial-metada.csv" where it had sampling date as one of the fields but the "0.metadata_CABBI_LAMPS_DNA.xlsx" file did not. I am using both as a source to reconstruct the sequence names and properly annotate the samples.
The two files contain overlapping but not identical information:
- "bacterial-metada.csv": holds experiment- and lab-level fields (sample IDs, plot, N-rate, etc.).
- "0.metadata_CABBI_LAMPS_DNA.xlsx": holds sequencing metadata and canonical attributes (e.g., plant species) to normalize and disambiguate samples.

## R workflow (overview)

1) Read both files (CSV and Excel).
2) Trim/clean join keys (e.g., "sample_id"), align column names as needed.
3) Left-join Excel columns into the CSV by the shared key.
4) Normalize fields (e.g., map plant to tokens C/M; handle “Reextracted-” prefix).
5) Correct inverted field orders in sequence file names
6) Validate metadata against actual sequence files: check uniqueness, missing values, row counts.
7) Write the combined metadata for reuse.


## Sample Name Standardization
During data processing, we encountered sample names with inconsistent field ordering that required standardization:

Issue 1: Nitrogen and Plot Order

Inconsistent: 2095_M_2016_N0_P39_B_20180430_RNA (N before P)
Standardized: 2095_M_2016_P39_N0_B_20180430_RNA (P before N)

Issue 2: BULK_DNA Position

Inconsistent: 3109_C_2018_P14_N0_20180514_BULK_DNA (BULK_DNA after date)
Standardized: 3109_C_2018_P14_N0_BULK_DNA_20180514 (BULK_DNA before date)

Note:
We also have "Bulk": 2163_C_2018_P14_N0_Bulk_20180429_CABBI.R1.fastq (Sentence case)

Rationale: These inconsistencies arose from different naming conventions used across sequencing batches or sample preparation workflows. Standardizing the field order ensures:

### Dropped strings
We eliminated the following trailing string after unique sample names.

DNA samples: "_CABBI.R[12]"
RNA samples: "*_CABBI_L001_R[12]_001"


The standardization process in "01_import.R" conditionally transforms only the incorrectly formatted names while preserving already correctly formatted ones.

## Outputs

- A single, cleaned metadata table ("*_updated_metadata") suitable for downstream R analyses and consistent filename construction.


## Reproducibility notes

- Keep the two source files versioned (CSV) or archived (Excel from Box).
- Record any column renames made during import so the join remains stable if schemas change.


# LAMPS 2023 (AKA LAMPS_EPS)

## Sequence and Metadata provenance

  - Metadata: "GERMS-DATAMAN/DOE-CABBI/Phillip/Projects/LAMPS_EPS/data/LAMPS_EPS_metadata.xlsx"
  - Raw sequences: Original sequences in 
      - 16S: "GERMS-DATAMAN/Argonne_Sequencing_Results/Mar2023_16S"
      - AMF: "GERMS-DATAMAN/Argonne_Sequencing_Results/Feb2025_AMF"
  - .RData (phyloseq objects): "GERMS-DATAMAN/DOE-CABBI/Phillip/LAMPS_EPS/data/"

  **Now in "~germs_miscanthus/data/input/LAMPS/2023/data". This is now categorized by sampling year.**

