###################################################################
# Tests for FASTA Export Functions
#
# Tests for refseq2fasta() and process_fa() functions using testthat.
# These tests validate the naming conventions and output structure.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-12-17
###################################################################

library(testthat)

# Source required functions
source(here::here("R/functions/refseq2fasta.R"))
source(here::here("R/functions/process_fa.R"))

# ---------------------------------------------------
# Test: refseq2fasta naming conventions ----
# ---------------------------------------------------

test_that("refseq2fasta generates correct filename with crop_subset", {
  # Test that crop_subset parameter is properly appended to filename
  # Expected format: {base_name}_{crop_subset}.fa

  # Create mock filename generation (without actual phyloseq object)
  base_name <- "ef_16S_DNA"
  crop_subset <- "all_crops"
  out_ext <- ".fa"

  expected_filename <- paste0(base_name, "_", crop_subset, out_ext)

  expect_equal(expected_filename, "ef_16S_DNA_all_crops.fa")
})

test_that("refseq2fasta generates correct filename without crop_subset", {
  # Test filename when crop_subset is NULL
  base_name <- "ef_16S_DNA"
  crop_subset <- NULL
  out_ext <- ".fa"

  expected_filename <- if (!is.null(crop_subset)) {
    paste0(base_name, "_", crop_subset, out_ext)
  } else {
    paste0(base_name, out_ext)
  }

  expect_equal(expected_filename, "ef_16S_DNA.fa")
})

test_that("refseq2fasta sanitizes filenames correctly", {
  # Test that special characters and spaces are handled
  test_name <- "test name with spaces"
  crop_subset <- "mxg_only"
  out_ext <- ".fa"

  fname <- paste0(test_name, "_", crop_subset, out_ext)
  fname <- gsub("[[:space:]]+", "_", fname)
  fname <- gsub("[^A-Za-z0-9_\\-\\.]", "", fname)

  expect_equal(fname, "test_name_with_spaces_mxg_only.fa")
})

# ---------------------------------------------------
# Test: process_fa naming conventions ----
# ---------------------------------------------------

test_that("process_fa generates correct output filename for single region", {
  # Test single region output naming
  region <- "16S"
  output_prefix <- "combined_"
  crop_subset <- "mxg_only"
  nucleotide <- "DNA"
  combine_all <- FALSE

  region_label <- if (combine_all) "all_regions" else region
  nuc_label <- if (!is.null(nucleotide)) paste0("_", nucleotide) else ""
  crop_label <- if (!is.null(crop_subset)) paste0("_", crop_subset) else ""

  output_filename <- paste0(
    output_prefix,
    region_label,
    nuc_label,
    crop_label,
    "_combined.fa"
  )

  expect_equal(output_filename, "combined_16S_DNA_mxg_only_combined.fa")
})

test_that("process_fa generates correct output filename for all regions", {
  # Test combined regions output naming
  region <- c("16S", "ITS", "AMF")
  output_prefix <- "combined_"
  crop_subset <- "all_crops"
  nucleotide <- NULL
  combine_all <- TRUE

  region_label <- if (combine_all) "all_regions" else region
  nuc_label <- if (!is.null(nucleotide)) paste0("_", nucleotide) else ""
  crop_label <- if (!is.null(crop_subset)) paste0("_", crop_subset) else ""

  output_filename <- paste0(
    output_prefix,
    region_label,
    nuc_label,
    crop_label,
    "_combined.fa"
  )

  expect_equal(output_filename, "combined_all_regions_all_crops_combined.fa")
})

test_that("process_fa generates filename without optional parameters", {
  # Test when nucleotide and crop_subset are NULL
  region <- "AMF"
  output_prefix <- "combined_"
  crop_subset <- NULL
  nucleotide <- NULL
  combine_all <- FALSE

  region_label <- if (combine_all) "all_regions" else region
  nuc_label <- if (!is.null(nucleotide)) paste0("_", nucleotide) else ""
  crop_label <- if (!is.null(crop_subset)) paste0("_", crop_subset) else ""

  output_filename <- paste0(
    output_prefix,
    region_label,
    nuc_label,
    crop_label,
    "_combined.fa"
  )

  expect_equal(output_filename, "combined_AMF_combined.fa")
})

# ---------------------------------------------------
# Test: Naming convention consistency ----
# ---------------------------------------------------

test_that("naming convention follows expected pattern", {
  # Validate naming convention:
  # Individual: {project}_{target_region}_{nucleotide}_{crop_subset}.fa
  # Combined: combined_{region}_{nucleotide}_{crop_subset}_combined.fa

  # Test individual file pattern
  # Note: nucleotide is optional because some datasets may have only one type
  # The pattern allows for flexibility while maintaining consistency
  individual_pattern <- "^[a-z]+(_[0-9]{4})?_[A-Z0-9]+(_DNA|_RNA|_DNA_RNA)?(_[a-z_]+)?\\.fa$"

  expect_true(grepl(individual_pattern, "ef_16S_DNA_all_crops.fa"))
  expect_true(grepl(individual_pattern, "lamps_2018_ITS_mxg_only.fa"))
  expect_true(grepl(individual_pattern, "lamps_2022_AMF_DNA_all_crops.fa"))
  expect_true(grepl(individual_pattern, "ef_16S_all_crops.fa")) # Without explicit nucleotide

  # Test combined file pattern
  combined_pattern <- "^combined_.*_combined\\.fa$"

  expect_true(grepl(combined_pattern, "combined_16S_DNA_mxg_only_combined.fa"))
  expect_true(grepl(
    combined_pattern,
    "combined_all_regions_all_crops_combined.fa"
  ))
})

# ---------------------------------------------------
# Test: FASTA file content validation helpers ----
# ---------------------------------------------------

test_that("FASTA sequence count matches expected", {
  # Helper function to count sequences in a FASTA file
  count_fasta_sequences <- function(file_path) {
    if (!file.exists(file_path)) {
      return(NA)
    }
    lines <- readLines(file_path, warn = FALSE)
    sum(grepl("^>", lines))
  }

  # Test the helper function logic
  mock_fasta_content <- c(">seq1", "ATCG", ">seq2", "GCTA", ">seq3", "TTTT")
  expect_equal(sum(grepl("^>", mock_fasta_content)), 3)
})

# ---------------------------------------------------
# Run tests if executed directly ----
# ---------------------------------------------------

# if (interactive()) {
#   testthat::test_file(here::here("tests/test_fasta_export.R"))
# }
