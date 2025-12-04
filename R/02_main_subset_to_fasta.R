###################################################################
# From Phyloseq to FASTA
#
# Script for building main phyloseq objest used throughout the project. Use of lists to organize phyloseq objects and projects.

# It also extracts and concatenates sequences and exports them to FASTA. In the process, the renamed and concatenated sequences are reassigned to the phyloseq object.

# Author: Bolívar Aponte Rolón
# Date: 2025-10-23
###################################################################

source("R/utils/00_setup.R")

load("data/output/rdata/phyloseq/lamps_2018_physeq_list.rda")
load("data/output/rdata/phyloseq/lamps_2022_physeq_list.rda")
load("data/output/rdata/phyloseq/ef_physeq_list.rda")

# ---------------------------------------------------
# Subset phyloseq by plant type
# ---------------------------------------------------

# To subset the actual phyloseq object by MXG
mxg_lamps_2018 <- purrr::map(
  lamps_2018_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "M")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
)

mxg_lamps_2022 <- purrr::map(
  lamps_2022_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "Miscanthus")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
)

mxg_ef <- purrr::map(
  ef_physeq_list,
  ~ {
    ps_subset <- subset_samples(.x, crop == "MXG")
    filter_taxa(ps_subset, function(x) sum(x > 0) > 0, TRUE)
  }
)

# Main list for downstream analyses
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
# ---------------------------------------------------
# Only MXG data from projects with all target regions
# ---------------------------------------------------

main_mxg_physeq_list <- list(
  mxg_ef = mxg_ef,
  mxg_lamps_2018 = mxg_lamps_2018,
  mxg_lamps_2022 = mxg_lamps_2022
)

save(
  main_mxg_physeq_list,
  file = "data/output/rdata/main_mxg_physeq_list.rda"
)

# ---------------------------------------------------
# Only 16S data for all projects
# ---------------------------------------------------

# Subset to only 16S DNA data
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

# Phyloseq to FASTA
# create output dir

# LAMPS --------------------------
# Full with all target regions
refseq2fasta(
  lamps_2018_physeq_list,
  out_dir = "data/output/sequences"
)

refseq2fasta(
  lamps_2022_physeq_list,
  out_dir = "data/output/sequences"
)

# MXG with all target regions
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

# Energy Farm Collab --------------------------
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


# Only 16S DNA data --------------------------
refseq2fasta(
  main_16S_physeq_list,
  out_dir = "data/output/sequences"
)
# In the case of ef_16S.fa = ef_16S_DNA.fa and lamps_2022_16S = lamps_16S_DNA.fa, the resulting files with the "_DNA" suffix should be the same as when output with "all target regions" since they only contain 1 nucleotide.

# ---------------------------------------------------
# Concatenate and export per target region
# ---------------------------------------------------

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

# -----------------------------------------------------------------
# Concatenating all otu tables for 16S DNA phyloseqs
# -----------------------------------------------------------------
ef_otu_table <- main_16S_physeq_list$ef_16S_DNA %>%
  otu_table() %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sample_id")

lamps_2018_otu_table <- main_16S_physeq_list$lamps_2018_16S_DNA %>%
  otu_table() %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sample_id")

lamps_2022_otu_table <- main_16S_physeq_list$lamps_2022_16S_DNA %>%
  otu_table() %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sample_id")


combined_seqtab <- full_join(
  {
    full_join(
      ef_otu_table,
      lamps_2018_otu_table,
      join_by("sample_id")
    )
  },
  lamps_2022_otu_table,
  join_by("sample_id")
)


# Replace resulting NAs
combined_seqtab_DT <- data.table::as.data.table(combined_seqtab)
for (j in seq_len(ncol(combined_seqtab_DT))) {
  data.table::set(
    combined_seqtab_DT,
    which(is.na(combined_seqtab_DT[[j]])),
    j,
    0
  )
}

# Ensuring types
combined_seqtab_DT[,
  (names(combined_seqtab_DT)) := lapply(.SD, function(x) {
    if (is.numeric(x)) as.integer(x) else x
  })
]

new_16S_DNA_seqtab <- combined_seqtab_DT %>%
  column_to_rownames(., var = "sample_id") %>%
  as.matrix()

save(new_16S_DNA_seqtab, file = "data/output/rdata/new_16S_DNS_seqtab.rda")
# TODO
# This new seqtab matrix is meant to be used to reassign taxonomy. Another approach is to take the taxonomy tables and reassing sequence ID to a concatenated table of all projects then subset by project.

!duplicated(colnames(new_16S_DNA_seqtab))
