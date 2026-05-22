###################################################################
# From Phyloseq to FASTA
#
# SECTION 1: Script for building main phyloseq objest used throughout the project. Use of lists to organize phyloseq objects and projects.
# Output:
# - main_physeq_list = list of phyloseq object with all the information (e.g. crops, nucleotides, target regions)
# - main_mxg_physeq_list = list of phyloseq objects containing only MXG crop with 16S_DNA, ITS and AMF target regions and nucleotides
# - main_16S_physeq_list = list of phyloseq object subsetted to 16S target region and DNA nucleotide for all crops across projects.

# SECTION 2: Extracts and concatenates sequences and exports them to FASTA. In the process, the renamed and concatenated sequences are reassigned to the corresponding phyloseq object.

# SECTION 3: Exporting and concatenating all otu tables for 16S DNA and fungal (ITS and AMF) phyloseqs to create a new sequence (OTU) table for reassigning taxonomy downstream.
# - 16S: Resulting union of otu tables yields 56395 AVSs, the same as the proccesd FASTA files from Section 2

# Author: Bolívar Aponte Rolón
# Date: 2025-10-23
###################################################################

source("R/utils/00_setup.R")

load("data/output/rdata/phyloseq/lamps_2018_physeq_list.rda")
load("data/output/rdata/phyloseq/lamps_2022_physeq_list.rda")
load("data/output/rdata/phyloseq/ef_physeq_list.rda")


# SECTION 1: Subset phyloseq ----

## Main list for downstream analyses ----
# All projects/datasets with all target regions, nucleotides and crops.
main_physeq_list <- list(
  ef_physeq_list = ef_physeq_list,
  lamps_2018_physeq_list = lamps_2018_physeq_list,
  lamps_2022_physeq_list = lamps_2022_physeq_list
)

save(
  main_physeq_list,
  file = "data/output/rdata/main_physeq_list_02.rda"
)

## Only 16S data for all projects ----
## The idea here is to use this list of object for developing the pipeline. We can then explore fungal datasets.

# Subset to only 16S DNA data
## lamps_2018_16S is the only one that needs it.
lamps_2018_16S_DNA <- lamps_2018_physeq_list$lamps_2018_16S_physeq %>%
  subset_samples(., nucleotide == "DNA") %>%
  filter_taxa(., function(x) sum(x > 0) > 0, TRUE)


main_16S_physeq_list <- list(
  ef_16S_DNA = ef_physeq_list$ef_16S_physeq,
  lamps_2018_16S_DNA = lamps_2018_16S_DNA,
  lamps_2022_16S_DNA = lamps_2022_physeq_list$lamps_2022_16S_physeq
)

## MXG physeq lists ----
# To subset the actual phyloseq object by MXG
## With 16S, only DNA nucleotide

crop_patterns <- c("MXG", "M", "Miscanthus")

main_mxg_physeq_list <- purrr::imap(
  main_physeq_list,
  function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      # Subset to crops first
      ps_subset <- subset_samples(physeq_obj, crop %in% crop_patterns) %>%
        filter_taxa(function(x) sum(x) > 0, TRUE)

      # Subset to only 16S DNA nucleotide (check if nucleotide column exists)
      if ("nucleotide" %in% colnames(sample_data(ps_subset))) {
        ps_subset <- ps_subset %>%
          subset_samples(nucleotide == "DNA") %>%
          filter_taxa(function(x) sum(x > 0) > 0, TRUE)
      }

      return(ps_subset)
    }) %>%
      # Keep only non-NULL objects and rename
      purrr::discard(is.null) %>%
      set_names(names(.) %>% str_replace("_physeq$", "_DNA"))
  }
) %>%
  set_names(names(.) %>% str_remove("_physeq_list$"))


# Need this one to keep DNA and RNA nucleotides for FASTA
mxg_lamps_2018 <- purrr::map(
  lamps_2018_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "M")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
) |>
  set_names(names(lamps_2018_physeq_list) |> str_remove("_physeq$"))


# save( # Replaced in 05_rebuild_and_transform.R
#   main_16S_mxg_physeq_list,
#   file = "data/output/rdata/main_16S_mxg_physeq_list_02.rda"
# )

# SECTION 2: Phyloseq to FASTA ----
#
# FASTA Naming Convention:
# {project}_{target_region}_{nucleotide}_{crop_subset}.fa
#
# Components:
# - project: ef, lamps_2018, lamps_2022
# - target_region: 16S, ITS, AMF
# - nucleotide: DNA, RNA (or DNA_RNA for combined)
# - crop_subset: all_crops, mxg_only
#
# Examples:
# - ef_16S_DNA_all_crops.fa        (Energy Farm, 16S, DNA, all crops)
# - lamps_2018_ITS_DNA_mxg_only.fa (LAMPS 2018, ITS, DNA, Miscanthus only)
# - lamps_2022_AMF_DNA_all_crops.fa (LAMPS 2022, AMF, DNA, all crops)
#
# Combined file naming:
# combined_{crop_subset}_{region}_{nucleotide}_combined.fa
# - combined_mxg_only_16S_DNA_combined.fa (MXG 16S DNA sequences combined)
# - combined_all_crops_all_regions_DNA_combined.fa (All crops, all regions)

fasta_output_dir <- "data/output/sequences"

## 2.1 FASTA Export: All Crops ----
# Export sequences from all projects with all crops included

# Energy Farm - all target regions (16S_DNA, AMF)
refseq2fasta(
  ef_physeq_list,
  crop_subset = "all_crops",
  out_dir = fasta_output_dir
)

# LAMPS 2018 - all target regions (16S has DNA+RNA, ITS has DNA)
# Note: lamps_2018_16S includes both DNA and RNA nucleotides
refseq2fasta(
  lamps_2018_physeq_list,
  crop_subset = "all_crops",
  out_dir = fasta_output_dir
)

refseq2fasta(
  lamps_2018_16S_DNA,
  crop_subset = "all_crops",
  out_dir = fasta_output_dir
)

# LAMPS 2022 - all target regions (16S_DNA, AMF)
refseq2fasta(
  lamps_2022_physeq_list,
  crop_subset = "all_crops",
  out_dir = fasta_output_dir
)


## 2.2 FASTA Export: Miscanthus (MXG) Only ----
# Export sequences subsetted to Miscanthus crop only
# This list only contains DNA nucleotide for all target regions.

# Energy Farm MXG subset
refseq2fasta(
  main_mxg_physeq_list$ef,
  crop_subset = "mxg_only",
  out_dir = fasta_output_dir
)

# LAMPS 2018 MXG subset
refseq2fasta(
  main_mxg_physeq_list$lamps_2018,
  crop_subset = "mxg_only",
  out_dir = fasta_output_dir
)


## LAMPS 2018 MXG subset with 16S DNA+RNA
refseq2fasta(
  mxg_lamps_2018$lamps_2018_16S,
  crop_subset = "DNA_RNA_mxg_only",
  out_dir = fasta_output_dir
)

# LAMPS 2022 MXG subset
refseq2fasta(
  main_mxg_physeq_list$lamps_2022,
  crop_subset = "mxg_only",
  out_dir = fasta_output_dir
)

## 2.3 FASTA Export: 16S DNA Only (All Crops) ----
# Export 16S DNA sequences from the curated 16S-only list

refseq2fasta(
  main_16S_physeq_list,
  crop_subset = "all_crops",
  out_dir = fasta_output_dir
)
# | Project | File Name | Equivalence |
# |---------|-----------|-------------|
# | LAMPS 2022 | lamps_2022_16S_DNA_all_crops.fa | Same as lamps_2022_16S_all_crops.fa |
# | Energy Farm | ef_16S_DNA_all_crops.fa | Same as ef_16S_all_crops.fa |
# |
# | Note: 16S_DNA is the only target region in these projects |

# SECTION 2.4: FASTA Concatenation ----

# Combine FASTA files across projects for downstream analysis
# process_fa() reads existing FASTA files, deduplicates, and exports combined files

target_regions <- c("16S", "ITS", "AMF")

## 2.4.1 Combined FASTA: MXG All Regions ----
# Combine all MXG sequences across projects and all target regions

results_mxg_all <- process_fa(
  region = target_regions,
  path = fasta_output_dir,
  exclude_str = "_RNA",
  crop_subset = "mxg_only",
  output_prefix = "combined_",
  combine_all = TRUE,
  rename_headers = FALSE
)

## 2.4.2 Combined FASTA: MXG by Region ----
# Combine MXG sequences across projects, one file per target region

results_mxg_by_region <- purrr::map(
  target_regions,
  function(region) {
    process_fa(
      region = region,
      nucleotide = "DNA",
      exclude_str = "_RNA",
      path = fasta_output_dir,
      crop_subset = "mxg_only",
      output_prefix = "combined_",
      combine_all = FALSE,
      rename_headers = FALSE
    )
  }
) %>%
  purrr::set_names(target_regions)

# Summary of MXG by region
purrr::map_dfr(
  results_mxg_by_region,
  ~ data.frame(
    region = .x$region,
    n_files = .x$n_files,
    n_unique_sequences = .x$n_sequences,
    output_file = basename(.x$output_path)
  ),
  .id = "target"
)

## 2.4.3 Combined FASTA: 16S DNA All Crops ----
# Combine 16S DNA sequences across all projects and crops

results_16S_DNA_all_crops <- process_fa(
  region = "16S",
  path = fasta_output_dir,
  vec_files = c(
    "ef_16S_all_crops.fa",
    "lamps_2018_16S_DNA_all_crops.fa",
    "lamps_2022_16S_all_crops.fa"
  ),
  exclude_str = "_RNA",
  crop_subset = "all_crops",
  output_prefix = "combined_",
  extra_id = "DNA",
  rename_headers = FALSE
)


# SECTION 3: OTU_TABLE Exports ----

# Export all otu tables for 16S DNA phyloseqs

ef_16S_otu_table <- export_otu_table(main_16S_physeq_list$ef_16S_DNA)

lamps_2018_16S_otu_table <- export_otu_table(
  main_16S_physeq_list$lamps_2018_16S_DNA
)

lamps_2022_16S_otu_table <- export_otu_table(
  main_16S_physeq_list$lamps_2022_16S_DNA
)

# Export OTU tables for Fungal datasets only MXG crop

ef_mxg_AMF <- export_otu_table(main_mxg_physeq_list$ef$ef_AMF)

lamps_2018_mxg_ITS <- export_otu_table(
  main_mxg_physeq_list$lamps_2018$lamps_2018_ITS
)
lamps_2022_mxg_AMF <- export_otu_table(
  main_mxg_physeq_list$lamps_2022$lamps_2022_AMF
)

## Do we have shared ASV sequences? ----
# union, intersect and setdiff discard any duplicated values in the arguments

A <- colnames(ef_16S_otu_table)
B <- colnames(lamps_2018_16S_otu_table)
C <- colnames(lamps_2022_16S_otu_table)

# basic intersections/unions
A_and_B <- intersect(A, B)
A_and_C <- intersect(A, C)
B_and_C <- intersect(B, C)
all_three <- Reduce(intersect, list(A, B, C)) # A ∩ B ∩ C
union_all <- Reduce(union, list(A, B, C)) # A ∪ B ∪ C

# # "only" sets (in one set and not any other)
# A_only <- setdiff(A, union(B, C))
# B_only <- setdiff(B, union(A, C))
# C_only <- setdiff(C, union(A, B))

# # pairwise-only (in exactly two sets)
# A_and_B_only <- setdiff(A_and_B, C) # in A and B but NOT in C
# A_and_C_only <- setdiff(A_and_C, B)
# B_and_C_only <- setdiff(B_and_C, A)

# # assemble into a named list
# seq_combinations <- list(
#   A_only = A_only,
#   B_only = B_only,
#   C_only = C_only,
#   A_and_B_only = A_and_B_only,
#   A_and_C_only = A_and_C_only,
#   B_and_C_only = B_and_C_only,
#   all_three = all_three,
#   union_all = union_all,
#   A_and_B = A_and_B,
#   A_and_C = A_and_C,
#   B_and_C = B_and_C
# )

## Use the intersection to by = NULL. From the documentation: If NULL, the default, *_join() will perform a natural join, using all variables in common across x and y.

# The combined DT should be the union of all three OTU tables (21833+40189+5132) - minus their intersection (10759)

# Inline: join_otu_tables — full-joins three OTU tables, fills NAs with 0, returns integer matrix
join_otu_tables <- function(table_a, table_b, table_c) {
  combined_seqtab_DT <- suppressMessages({
    full_join(full_join(table_a, table_b, by = NULL), table_c, by = NULL)
  }) %>%
    data.table::as.data.table(.)

  for (j in seq_len(ncol(combined_seqtab_DT))) {
    data.table::set(
      combined_seqtab_DT,
      which(is.na(combined_seqtab_DT[[j]])),
      j,
      0
    )
  }

  combined_seqtab_DT[,
    (names(combined_seqtab_DT)) := lapply(.SD, function(x) {
      if (is.numeric(x)) as.integer(x) else x
    })
  ]

  combined_seqtab_DT %>%
    column_to_rownames(., var = "sample_id") %>%
    as.matrix()
}

## New 16S DNA joined otu table (aka seqtab)
new_16S_DNA_seqtab <- join_otu_tables(
  ef_16S_otu_table,
  lamps_2018_16S_otu_table,
  lamps_2022_16S_otu_table
)

save(new_16S_DNA_seqtab, file = "data/output/rdata/new_16S_DNS_seqtab_02.rda")


## New fungal joined otu tables

new_fungal_seqtab <- join_otu_tables(
  ef_mxg_AMF,
  lamps_2018_mxg_ITS,
  lamps_2022_mxg_AMF
)

save(new_fungal_seqtab, file = "data/output/rdata/new_fungal_seqtab_02.rda")

# NOTE ----
# This new seqtab matrix is meant to be used to reassign taxonomy. Another approach is to take the taxonomy tables and reassign sequence ID to a concatenated table of all projects then subset by project.
