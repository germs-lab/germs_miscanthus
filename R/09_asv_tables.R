##########################################################################
# ASV Taxonomic Tables: Shared ASVs Between Projects
#
# Creates formatted gt() tables with taxonomic information for shared ASVs
# identified at different occupancy thresholds across projects.
#
# Author: Bolívar Aponte Rolón
# Date: 2026-01-26
##########################################################################

source("R/utils/00_setup.R")
library(gt)
library(gtExtras)

# Load shared ASV data from 07_abund_occupancy.R
load("data/output/rdata/shared_asvs_by_threshold_07.rda")

# Get full taxonomy from phyloseq objects
get_asv_taxonomy <- function(asv_ids, physeq_obj) {
  tax_df <- tax_table(physeq_obj) %>%
    as.data.frame() %>%
    rownames_to_column("asv") %>%
    filter(asv %in% asv_ids) %>%
    distinct()

  return(tax_df)
}

# SECTION 1: Taxonomic tables for all shared ASVs in Miscanthus in all 3 projects ----

cat("\n", rep("=", 40), "\n")
cat("Creating Taxonomic Tables for Shared ASVs\n")
cat(rep("=", 40), "\n\n")

# Use the first phyloseq object for taxonomy reference
# All phyloseq objects should have consistent taxonomy since they were
# reassigned together in 05_rebuild_and_transform.R
ref_physeq <- main_16S_physeq_list[[1]]

# Verify taxonomy table exists
if (is.null(tax_table(ref_physeq, errorIfNULL = FALSE))) {
  stop("Reference phyloseq object does not contain a taxonomy table")
}

# Create tables for ASVs shared across all projects at each threshold
mxg_shared_asv_tables <- purrr::imap(
  mxg_shared_asvs_by_threshold,
  function(sharing_data, threshold_name) {
    # Extract ASVs shared by all projects
    all_projects_asvs <- sharing_data$all_projects

    if (length(all_projects_asvs) == 0) {
      cat("No ASVs shared across all projects at", threshold_name, "\n")
      return(NULL)
    }

    # Get taxonomy for these ASVs
    tax_df <- get_asv_taxonomy(all_projects_asvs, ref_physeq)

    # Count how many projects each ASV is in
    asv_counts <- sharing_data$asv_project_counts[all_projects_asvs]

    # Combine with counts
    tax_summary <- tax_df %>%
      mutate(n_projects = asv_counts[asv]) %>%
      select(asv, kingdom, phylum, class, order, family, genus, n_projects) %>%
      arrange(phylum, class, order, family, genus)

    # Create gt table
    gt_table <- tax_summary %>%
      gt() %>%
      tab_header(
        title = "Shared Miscanthus ASVs Across All Projects",
        subtitle = paste0("ASVs at ", threshold_name, " occupancy threshold")
      ) %>%
      cols_label(
        asv = "ASV ID",
        kingdom = "Kingdom",
        phylum = "Phylum",
        class = "Class",
        order = "Order",
        family = "Family",
        genus = "Genus",
        n_projects = "# Projects"
      ) %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      )

    # Only add alternating row colors if there are at least 2 rows
    if (nrow(tax_summary) >= 2) {
      gt_table <- gt_table %>%
        tab_style(
          style = cell_fill(color = "#E8F4F8"),
          locations = cells_body(rows = seq(2, nrow(tax_summary), by = 2))
        )
    }

    gt_table <- gt_table %>%
      cols_width(
        asv ~ px(150),
        kingdom ~ px(100),
        phylum ~ px(150),
        class ~ px(150),
        order ~ px(150),
        family ~ px(150),
        genus ~ px(150),
        n_projects ~ px(80)
      ) %>%
      tab_options(
        table.font.size = 12,
        heading.title.font.size = 16,
        heading.subtitle.font.size = 14,
        column_labels.font.weight = "bold"
      )

    # Save table as HTML
    gtsave(
      gt_table,
      filename = paste0(
        "data/output/figures/mxg_shared_asvs_",
        threshold_name,
        "_taxonomy_table.html"
      )
    )

    cat(
      "Created table for",
      threshold_name,
      ":",
      nrow(tax_summary),
      "shared ASVs\n"
    )

    return(list(
      data = tax_summary,
      gt_table = gt_table
    ))
  }
)

# SECTION 2: Taxonomic tables for Miscanthus ASVs shared by 2+ projects ----

shared_2plus_tables <- purrr::imap(
  mxg_shared_asvs_by_threshold,
  function(sharing_data, threshold_name) {
    # Extract ASVs shared by 2 or more projects
    two_plus_asvs <- sharing_data$two_or_more

    if (length(two_plus_asvs) == 0) {
      cat("No ASVs shared by 2+ projects at", threshold_name, "\n")
      return(NULL)
    }

    # Get taxonomy for these ASVs
    tax_df <- get_asv_taxonomy(two_plus_asvs, ref_physeq)

    # Count how many projects each ASV is in
    asv_counts <- sharing_data$asv_project_counts[two_plus_asvs]

    # Combine with counts
    tax_summary <- tax_df %>%
      mutate(n_projects = asv_counts[asv]) %>%
      select(asv, kingdom, phylum, class, order, family, genus, n_projects) %>%
      arrange(desc(n_projects), phylum, class, order, family, genus)

    # Create gt table
    gt_table <- tax_summary %>%
      gt() %>%
      tab_header(
        title = "Miscanthus ASVs Shared by Two or More Projects",
        subtitle = paste0("ASVs at ", threshold_name, " occupancy threshold")
      ) %>%
      cols_label(
        asv = "ASV ID",
        kingdom = "Kingdom",
        phylum = "Phylum",
        class = "Class",
        order = "Order",
        family = "Family",
        genus = "Genus",
        n_projects = "# Projects"
      ) %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
      ) %>%
      tab_style(
        style = cell_fill(color = "#E8F4F8"),
        locations = cells_body(rows = seq(2, nrow(tax_summary), 2))
      ) %>%
      data_color(
        columns = n_projects,
        fn = scales::col_numeric(
          palette = c("#FFF7BC", "#FEC44F", "#D95F0E"),
          domain = c(2, length(sharing_data$core_by_project))
        )
      ) %>%
      cols_width(
        asv ~ px(150),
        kingdom ~ px(100),
        phylum ~ px(150),
        class ~ px(150),
        order ~ px(150),
        family ~ px(150),
        genus ~ px(150),
        n_projects ~ px(80)
      ) %>%
      tab_options(
        table.font.size = 12,
        heading.title.font.size = 16,
        heading.subtitle.font.size = 14,
        column_labels.font.weight = "bold"
      )

    # Save table as HTML
    gtsave(
      gt_table,
      filename = paste0(
        "data/output/figures/mxg_shared_asvs_2plus_",
        threshold_name,
        "_taxonomy_table.html"
      )
    )

    cat(
      "Created 2+ table for",
      threshold_name,
      ":",
      nrow(tax_summary),
      "ASVs\n"
    )

    return(list(
      data = tax_summary,
      gt_table = gt_table
    ))
  }
)

# SECTION 3: Summary statistics table ----

# Create a summary table showing counts at each threshold
summary_data <- purrr::imap_dfr(
  mxg_shared_asvs_by_threshold,
  function(sharing_data, threshold_name) {
    tibble(
      threshold = gsub("threshold_", "", threshold_name),
      all_projects = length(sharing_data$all_projects),
      two_or_more = length(sharing_data$two_or_more),
      total_unique = length(sharing_data$asv_project_counts)
    )
  }
)

summary_gt <- summary_data %>%
  gt() %>%
  tab_header(
    title = "Summary of Shared Miscanthus ASVs Across Thresholds",
    subtitle = "Number of ASVs identified at different occupancy thresholds"
  ) %>%
  cols_label(
    threshold = "Occupancy Threshold",
    all_projects = "All Projects",
    two_or_more = "2+ Projects",
    total_unique = "Total Unique ASVs"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_fill(color = "#E8F4F8"),
    locations = cells_body(rows = seq(2, nrow(summary_data), by = 2))
  ) %>%
  cols_width(
    threshold ~ px(150),
    all_projects ~ px(120),
    two_or_more ~ px(120),
    total_unique ~ px(150)
  ) %>%
  tab_options(
    table.font.size = 14,
    heading.title.font.size = 18,
    heading.subtitle.font.size = 14,
    column_labels.font.weight = "bold"
  )

# Save summary table
gtsave(
  summary_gt,
  filename = "data/output/figures/mxg_shared_asvs_summary_table.html"
)


cat("\n", rep("=", 40), "\n")
cat("All taxonomic tables saved to data/output/figures/\n")
cat(rep("=", 40), "\n")

# Display the first table for viewing
# if (!is.null(mxg_shared_asv_tables[[1]])) {
#   cat("\nExample table (first threshold):\n")
#   print(mxg_shared_asv_tables[[1]]$gt_table)
# }
