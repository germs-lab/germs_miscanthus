# Building proper metadata in R from two sources

This document explains how the final, “proper” metadata table was built entirely in R by combining two required inputs:

Both files were necessary, although they had more or less the same data. Metadata and sequence files could very well have been matched using sample_id(id) since each is unique, but "replicate" information was important as well as "sample_date" which the RNA metadata was missing.

This applies to DNA (DNA amplicon sequencing) and RNA (RNA amplicon sequencing)sequences.


## Metadata provenance
- Lab analysis from Micheal Millican
  - `GERMS-DATAMAN/Millican_GERMS_Research/LAMPS_microbiome/analysis_data/bacterial-metada.csv`
- Sequencing master metadata (Box)
  - `GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA amplicon sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx`
 or 
 `DOE-CABBI/bolivar/germs_miscanthus/data/input/LAMPS/DNA_amplicon_sequencing/0.metadata_CABBI_LAMPS_DNA.xlsx`



**Files are imported to `~/germs_miscanthus/data/input/LAMPS` from where they'll be loaded for further analyses.**

## Why two files?
The sequence file names in `~Dna amplicon sequencing/raw_sequences` are the results of the sequencing run. They were probably built using `bacterial-metada.csv` where it had sampling date as one of the fields but the `0.metadata_CABBI_LAMPS_DNA.xlsx` file did not. I am using both as a source to reconstruct the sequence names and properly annotate the samples.
The two files contain overlapping but not identical information:
- `bacterial-metada.csv`: holds experiment- and lab-level fields (sample IDs, plot, N-rate, etc.).
- `0.metadata_CABBI_LAMPS_DNA.xlsx`: holds sequencing metadata and canonical attributes (e.g., plant species) to normalize and disambiguate samples.

## R workflow (overview)

1) Read both files (CSV and Excel).
2) Trim/clean join keys (e.g., `sample_id`), align column names as needed.
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

DNA samples: `_CABBI.R[12]`
RNA samples: `*_CABBI_L001_R[12]_001`


The standardization process in `01_import.R` conditionally transforms only the incorrectly formatted names while preserving already correctly formatted ones.

## Outputs

- A single, cleaned metadata table (`*_proper_metadata`) suitable for downstream R analyses and consistent filename construction.


## Reproducibility notes

- Keep the two source files versioned (CSV) or archived (Excel from Box).
- Record any column renames made during import so the join remains stable if schemas change.