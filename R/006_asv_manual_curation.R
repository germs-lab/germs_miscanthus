#####################################################################
# Manual Curation of ASV blasted against NCBI nt database
#
# Curates BLASTn tabular results (outfmt 6) for amplicon ASVs into
# conservative taxonomic assignments.
#
# Filters out uninformative hits, enforces coverage/identity thresholds,
# resolves multiple hits via a tie window, and assigns the lowest
# reliable rank (species/genus/higher).
#
# Last Common Ancestor-like decision
# Purpose: When multiple top hits are essentially tied,
# assign the most specific taxonomic rank that all those hits consistently support.
#
# When it runs: After prefiltering, deduping, and defining the “tie set”
# of equally good hits (within a small bitscore/identity/coverage window).
#
# Decision logic:
#
# Species: If all tied hits agree on the same, fully named species
# (not “Genus sp.”) and the chosen hit meets stricter species thresholds
# (≈≥99% identity and ≥95% coverage), assign species.
#
# Genus: If species disagree but all genera agree, assign genus (species left blank).
#
# Higher: If tied hits span multiple genera (or evidence is insufficient), assign
# a higher/ambiguous rank.
#
# What’s reported: The script still picks a single “chosen” hit
# (for accession/title details), but assignment_rank and taxon_name come
# from the LCA consensus across the tie set. It also emits n_hits_tie and a
# decision_reason indicating which branch was taken.

# Scope/limits: This is a name-based consensus (from ssciname), not a full taxonomy-tree LCA using staxids; it doesn’t resolve synonyms or lineage conflicts beyond genus.

# Author: Bolívar Aponte Rolón
# Date: 2025-10-16
#####################################################################

# Load required libraries
source("R/utils/000_setup.R")

# Import results from BLAST
blast_results <- read_tsv(
  "data/output/processed/sequences/sabr_2023_corn_signif_asv_blast.tsv"
)

# Pre-checks
# Let's find duplicates, meaning they had multiple hits and we need to resolve them.
one_hits <- blast_results |>
  group_by(qseqid) |>
  filter(n() <= 1) |>
  ungroup()

multi_hits <- blast_results |>
  group_by(qseqid) |>
  filter(n() > 1) |>
  ungroup()

#--------------------------------------------------------
# Curation approach
#--------------------------------------------------------

# Parameters --------------------
MIN_QCOV <- as.numeric(0.90)
MIN_PIDENT_GENUS <- as.numeric(97)
MIN_PIDENT_SPECIES <- as.numeric(99)
MIN_QCOV_SPECIES <- as.numeric(0.95)
TIE_BITSCORE <- as.numeric(2)
TIE_PIDENT <- as.numeric(0.2) # percent points
TIE_QCOV <- as.numeric(0.05) # fraction (5%)

# Title/record exclusion patterns (case-insensitive)
EXCLUDE_REGEX <- "uncultured|Uncultured|environmental|metagenom|clone|unidentified|
  mixed culture|synthetic construct|vector|plasmid|chloroplast|mitochondr|
  plastid|soil|Soil|marine|Marine|bacterium|Bacterium"

# Expected BLAST outfmt 6 fields used here
EXPECTED_FIELDS <- c(
  "qseqid",
  "sseqid",
  "pident",
  "length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore",
  "qlen",
  "staxids",
  "ssciname",
  "stitle"
)

# Field validation
if (!all(EXPECTED_FIELDS %in% names(blast_results))) {
  if (ncol(blast_results) == length(EXPECTED_FIELDS)) {
    names(blast_results) <- EXPECTED_FIELDS
  } else {
    stop(
      "Input does not have expected headers and column count doesn't match expected fields."
    )
  }
}

# Pipeline --------------------

# Normalize types, compute coverage, extract accession, flag species/genus
df1 <- blast_results %>%
  mutate(
    pident = as.numeric(pident),
    bitscore = as.numeric(bitscore),
    qlen = as.numeric(qlen),
    length = as.numeric(length),
    evalue = as.numeric(evalue),
    qcov = if_else(qlen > 0, length / qlen, NA_real_),

    # Extract accession from sseqid if pipe-delimited; else keep sseqid
    sacc = case_when(
      str_detect(sseqid, "\\|") ~
        {
          parts <- str_split(sseqid, "\\|", n = 5, simplify = TRUE)
          if (ncol(parts) >= 4) parts[, 4] else sseqid
        },
      TRUE ~ sseqid
    ),

    # Flags for prioritization and assignment
    species_named = str_detect(
      stitle,
      "^[A-Za-z][A-Za-z-]*\\s+[a-z][A-Za-z0-9._-]*$"
    ) &
      !str_detect(stitle, "\\bsp\\.?\\b"),
    genus = str_replace_na(str_match(stitle, "^([A-Za-z][A-Za-z-]*)\\b")[,
      2
    ]),
    species = str_replace_na(str_match(
      stitle,
      "^[A-Za-z][A-Za-z-]*\\s+([a-z][A-Za-z0-9._-]*)\\b"
    )[, 2]),
    has_type = str_detect(stitle, "(?i)strain"),
    is_refseq_rrna = str_detect(sacc, "^NR_\\d+"),
    has_binomial = species_named
  ) %>%
  # Exclude uninformative or contaminant titles
  filter(!str_detect(stitle, EXCLUDE_REGEX)) %>%
  # Basic sanity
  filter(!is.na(pident), !is.na(qcov), !is.na(bitscore)) %>%
  # Amplicon QC thresholds (genus-level floor)
  filter(qcov >= MIN_QCOV, pident >= MIN_PIDENT_GENUS) %>%
  # Collapse duplicates (per qseqid, sacc, staxids) keeping best bitscore
  group_by(qseqid, sacc, staxids) %>%
  arrange(desc(bitscore), desc(pident), desc(qcov), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# Determine top metrics per query (for tie window relative to the top hit per our sort)
df2 <- df1 %>%
  group_by(qseqid) %>%
  arrange(desc(bitscore), desc(pident), desc(qcov), .by_group = TRUE) %>%
  mutate(
    top_bitscore = dplyr::first(bitscore),
    top_pident = dplyr::first(pident),
    top_qcov = dplyr::first(qcov)
  ) %>%
  ungroup()

# Tie set: within Δbitscore ≤ TIE_BITSCORE OR within both pident/qcov windows
tie_set <- df2 %>%
  filter(
    (bitscore >= (top_bitscore - TIE_BITSCORE)) |
      (abs(pident - top_pident) <= TIE_PIDENT &
        abs(qcov - top_qcov) <= TIE_QCOV)
  )

# Summaries over tie set (for LCA-style decision)
tie_summ <- tie_set %>%
  mutate(genus_in_tie = genus) %>%
  group_by(qseqid) %>%
  summarise(
    n_hits_tie = n(),
    n_species_named_distinct = n_distinct(
      species[species_named],
      na.rm = TRUE
    ),
    n_genus_distinct = n_distinct(genus_in_tie, na.rm = TRUE),
    .groups = "drop"
  )

# Chose single best within tie set (priority ordering)
chosen_asv <- tie_set %>%
  group_by(qseqid) %>%
  arrange(
    desc(has_type),
    desc(is_refseq_rrna),
    desc(has_binomial),
    desc(qcov),
    desc(pident),
    desc(bitscore),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    qseqid,
    sacc,
    sseqid,
    staxids,
    ssciname,
    stitle,
    pident,
    qcov,
    bitscore,
    species_named,
    genus,
    species
  ) %>%
  rename(
    chosen_sacc = sacc,
    chosen_sseqid = sseqid,
    chosen_staxids = staxids,
    chosen_ssciname = ssciname,
    chosen_stitle = stitle
  )

# Counts of considered hits per query (after prefiltering/collapse)
n_considered <- df1 %>%
  dplyr::count(qseqid, name = "n_hits_considered")

# Join summaries
final_asv <- chosen_asv %>%
  dplyr::left_join(tie_summ, by = "qseqid") %>%
  dplyr::left_join(n_considered, by = "qseqid") %>%
  mutate(
    species_ok = species_named &
      pident >= MIN_PIDENT_SPECIES &
      qcov >= MIN_QCOV_SPECIES,
    assignment_rank = case_when(
      !is.na(n_species_named_distinct) &
        n_species_named_distinct == 1 &
        species_ok ~
        "species",
      !is.na(n_genus_distinct) & n_genus_distinct == 1 & !is.na(genus) ~
        "genus",
      TRUE ~ "higher"
    ),
    taxon_name = case_when(
      assignment_rank == "species" ~ species,
      assignment_rank == "genus" ~ genus,
      TRUE ~ NA_character_
    ),
    decision_reason = case_when(
      assignment_rank == "species" ~ "clear_best_or_consensus_species",
      assignment_rank == "genus" ~ "tie_LCA_genus",
      TRUE ~ "tie_LCA_higher_or_conflict"
    )
  ) %>%
  select(
    qseqid,
    chosen_sacc,
    chosen_sseqid,
    chosen_staxids,
    chosen_ssciname,
    chosen_stitle,
    pident,
    qcov,
    bitscore,
    n_hits_considered,
    n_hits_tie,
    assignment_rank,
    taxon_name,
    decision_reason
  ) %>%
  arrange(qseqid)

#

# Write output --------------------

readr::write_tsv(
  final_asv,
  file = "data/output/processed/sequences/sabr_2023_corn_curated_asv_blast.tsv"
)
