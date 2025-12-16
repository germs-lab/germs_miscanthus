###################################################################
# From Phyloseq to FASTA
#
# SECTION 1: Script for building main phyloseq objest used throughout the project. Use of lists to organize phyloseq objects and projects. Output:
# - main_physeq_list = list of phyloseq object with all the information (e.g. crops, nucleotides, target regions)
# - main_mxg_physeq_list = list of phyloseq objects containing only MXG crop with all target regions and nucleotides
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

# ---------------------------------------------------
# SECTION 1: Subset phyloseq ----
# ---------------------------------------------------

## Main list for downstream analyses ----
# All projects/datasets with all target regions
main_physeq_list <- list(
  ef_physeq_list = ef_physeq_list,
  lamps_2018_physeq_list = lamps_2018_physeq_list,
  lamps_2022_physeq_list = lamps_2022_physeq_list
)

save(
  main_physeq_list,
  file = "data/output/rdata/main_physeq_list.rda"
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
save(
  main_16S_physeq_list,
  file = "data/output/rdata/main_16S_physeq_list.rda"
)

## Only MXG ----
# To subset the actual phyloseq object by MXG
mxg_ef <- purrr::imap(
  ef_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "MXG")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
) |>
  set_names(names(ef_physeq_list) |> str_remove("_physeq$"))

mxg_lamps_2018 <- purrr::map(
  lamps_2018_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "M")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
) |>
  set_names(names(lamps_2018_physeq_list) |> str_remove("_physeq$"))

mxg_lamps_2022 <- purrr::map(
  lamps_2022_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "Miscanthus")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
) |>
  set_names(names(lamps_2022_physeq_list) |> str_remove("_physeq$"))


## MXG physeq lists ----
## 16S is only DNA nucleotide
main_mxg_physeq_list <- list(
  ef = mxg_ef,
  lamps_2018 = list(
    lamps_2018_16S = lamps_2018_16S_DNA, # Making sure we only have 16S_DNA
    lamps_2018_ITS = mxg_lamps_2018$lamps_2018_ITS
  ),
  lamps_2022 = mxg_lamps_2022
)

save(
  main_mxg_physeq_list,
  file = "data/output/rdata/main_mxg_physeq_list.rda"
)


# ---------------------------------------------------
# SECTION 2: Phyloseq to FASTA ----
# ---------------------------------------------------

## LAMPS ----
# Full with all target regions
refseq2fasta(
  lamps_2018_physeq_list,
  out_dir = "data/output/sequences"
)

refseq2fasta(
  lamps_2022_physeq_list,
  out_dir = "data/output/sequences"
)

## MXG with all target regions ----
refseq2fasta(
  mxg_lamps_2018,
  extra_id = "_mxg",
  out_dir = "data/output/sequences"
)

refseq2fasta(
  mxg_lamps_2022,
  extra_id = "_mxg",
  out_dir = "data/output/sequences"
)

## Energy Farm Collab ----
# Fullwith all target regions
refseq2fasta(
  ef_physeq_list,
  out_dir = "data/output/sequences"
)

# MXG with all target regions
refseq2fasta(
  mxg_ef,
  extra_id = "_mxg",
  out_dir = "data/output/sequences"
)

# Only 16S DNA data ----
refseq2fasta(
  main_16S_physeq_list,
  out_dir = "data/output/sequences"
)
# NOTE ----
# In the case of ef_16S.fa = ef_16S_DNA.fa and lamps_2022_16S = lamps_16S_DNA.fa, the resulting files with the "_DNA" suffix should be the same as when output with "all target regions" since they only contain 1 nucleotide.

# ---------------------------------------------------
# FASTA Exports ----
# ---------------------------------------------------

# Concatenate and export per target region

target_regions <- c("16S", "ITS", "AMF")


# Export and summary
# process_fa() takes preexisting FASTA files when given file names, combines and reassigns sequence tags.

# MXG, regions and nucleotides together, exported as individual files
results <- purrr::map(
  target_regions,
  function(region) {
    process_fa(
      region,
      path = "data/output/sequences/",
      combined_suffix = "_combined_asv_renamed.fa",
      prefix = "mxg_",
      target_suffix = "_mxg.fa",
      new_headers = FALSE
    )
  }
) %>%
  set_names(target_regions)
purrr::map_dfr(results, ~ data.frame(.x), .id = "region")


# MXG, target regions and nucleotides together exported as combined file
result <- list(
  all_regions = process_fa(
    region = target_regions,
    path = "data/output/sequences/",
    combined_suffix = "_combined_asv_renamed.fa",
    prefix = "mxg_",
    target_suffix = "_mxg.fa",
    new_headers = FALSE,
    .all = "all"
  )
)
purrr::map_dfr(result, ~ data.frame(.x), .id = "region")


# 16S region and DNA nucleotide - all crops
results_16S_DNA <- purrr::map(
  "16S",
  function(region) {
    process_fa(
      region,
      path = "data/output/sequences/",
      combined_suffix = "_DNA_combined_asv_renamed.fa",
      prefix = "all_",
      target_suffix = "_DNA.fa",
      .all = "16S"
    )
  }
) %>%
  set_names("16S")
purrr::map_dfr(results_16S_DNA, ~ data.frame(.x), .id = "region")

# ---------------------------------------------------
# SECTION 3: OTU_TABLE Exports
# ---------------------------------------------------
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

# TODO
# MAke in to function and create new seqtab for fungal data
combined_16S_seqtab_DT <- full_join(
  {
    full_join(
      ef_16S_otu_table,
      lamps_2018_16S_otu_table,
      by = NULL
    )
  },
  lamps_2022_16S_otu_table,
  by = NULL
) %>%
  data.table::as.data.table(.) # We need some speed with this big data


# Replace resulting NAs
for (j in seq_len(ncol(combined_16S_seqtab_DT))) {
  data.table::set(
    combined_16S_seqtab_DT,
    which(is.na(combined_16S_seqtab_DT[[j]])),
    j,
    0
  )
}

# Ensuring types
combined_16S_seqtab_DT[,
  (names(combined_16S_seqtab_DT)) := lapply(.SD, function(x) {
    if (is.numeric(x)) as.integer(x) else x
  })
]

new_16S_DNA_seqtab <- combined_16S_seqtab_DT %>%
  column_to_rownames(., var = "sample_id") %>%
  as.matrix()

save(new_16S_DNA_seqtab, file = "data/output/rdata/new_16S_DNS_seqtab.rda")

# NOTE ----
# This new seqtab matrix is meant to be used to reassign taxonomy. Another approach is to take the taxonomy tables and reassign sequence ID to a concatenated table of all projects then subset by project.
