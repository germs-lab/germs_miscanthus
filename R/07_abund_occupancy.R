###########################################################################
# Abundance Occupancy Analysis on Miscanthus 16S
#
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-24
##########################################################################

source("R/utils/00_setup.R")

# Getting our cores
main_16S_core_list <- purrr::map(
  main_16S_physeq_list,
  function(physeq_obj) {
    # purrr::imap(project_list, function(physeq_obj, physeq_name) {
    get_prevalent_rare(
      physeq_obj,
      thresholds = c(90, 80, 70, 60, 30),
      detection = 0 / 100,
      include.lowest = FALSE
    )
  }
)
# }
#)

# Abundance-Occupancy Plots - All crops ----
abundance_occ_plots <- purrr::imap(
  main_16S_physeq_list,
  function(physeq_obj, physeq_name) {
    # Get summaries
    summarized_df_list <- summarize_abundance_occupancy(physeq_obj)

    # Define thresholds to test
    thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)

    # Plots at different thresholds - overall
    overall_plots <- purrr::map(thresholds, function(threshold) {
      plot_abundance_occupancy(
        summarized_df_list$overall,
        threshold = threshold,
        title = paste0(physeq_name, " - Overall (threshold:", threshold, ")")
      )
    }) %>%
      set_names(paste0("threshold_", thresholds))

    # Plots at different thresholds - by crop
    by_crop_plots <- purrr::imap(
      summarized_df_list$by_crop,
      function(crop_df, crop_name) {
        purrr::map(thresholds, function(threshold) {
          plot_abundance_occupancy(
            crop_df,
            threshold = threshold,
            title = paste0(
              physeq_name,
              " - ",
              crop_name,
              " (threshold:",
              threshold,
              ")"
            )
          )
        }) %>%
          set_names(paste0("threshold_", thresholds))
      }
    )

    return(list(
      summary_df = summarized_df_list,
      overall_plots = overall_plots,
      by_crop_plots = by_crop_plots
    ))
  }
)

abundance_occ_plots$ef_16S$overall_plots$threshold_0.6
abundance_occ_plots$ef_16S$by_crop_plots$MXG$threshold_0.6
abundance_occ_plots$ef_16S$by_crop_plots$SB$threshold_0.6
abundance_occ_plots$ef_16S$by_crop_plots$ZM$threshold_0.6

abundance_occ_plots$lamps_2018_16S$overall_plots$threshold_0.6
abundance_occ_plots$lamps_2018_16S$by_crop_plots$M$threshold_0.6
abundance_occ_plots$lamps_2018_16S$by_crop_plots$C$threshold_0.6

abundance_occ_plots$lamps_2022_16S$overall_plots$threshold_0.6
abundance_occ_plots$lamps_2022_16S$by_crop_plots$Miscanthus$threshold_0.6
abundance_occ_plots$lamps_2022_16S$by_crop_plots$Maize$threshold_0.6
abundance_occ_plots$lamps_2022_16S$by_crop_plots$Grass$threshold_0.6


# Identifying core ASVs at 60% threshold
ef_16S_core_asvs_60 <- rownames(
  summarize_abundance_occupancy(main_16S_physeq_list$ef_16S_DNA) %>%
    .$overall %>%
    filter(occupancy_prop >= 0.6)
)
lamps_2018_16S_core_asvs_60 <- rownames(
  summarize_abundance_occupancy(main_16S_physeq_list$lamps_2018_16S_DNA) %>%
    .$overall %>%
    filter(occupancy_prop >= 0.6)
)

(c(ef_16S_core_asvs_60, lamps_2018_16S_core_asvs_60))
diffs <- base::setdiff(
  ef_16S_core_asvs_60,
  lamps_2018_16S_core_asvs_60
)

# SECTION: Shared ASVs Between Projects at Different Thresholds ----

# Define thresholds to test
# These represent the proportion of samples where an ASV must be present
# to be considered "core" (e.g., 0.6 = present in 60% of samples)
# Modify this vector to analyze different thresholds
thresholds <- c(0.0, 0.6, 0.7, 0.8, 0.9)

# Get project names from the phyloseq list (focusing on DNA projects)
project_names <- grep("_DNA$", names(main_16S_physeq_list), value = TRUE)

# Extract core ASVs for each project at each threshold
core_asvs_by_threshold <- purrr::map(thresholds, function(threshold) {
  project_asvs <- purrr::map(project_names, function(proj_name) {
    rownames(
      summarize_abundance_occupancy(main_16S_physeq_list[[proj_name]]) %>%
        .$overall %>%
        filter(occupancy_prop >= threshold)
    )
  }) %>%
    set_names(gsub("_DNA$", "", project_names))

  return(project_asvs)
}) %>%
  set_names(paste0("threshold_", thresholds))

# Calculate shared ASVs between all projects at each threshold
shared_asvs_by_threshold <- purrr::imap(
  core_asvs_by_threshold,
  function(asv_list, threshold_name) {
    # Get all unique ASVs across all projects
    all_asvs <- unique(unlist(asv_list))

    # Count how many projects each ASV appears in
    asv_counts <- purrr::map_int(all_asvs, function(asv) {
      sum(purrr::map_lgl(asv_list, ~ asv %in% .x))
    })
    names(asv_counts) <- all_asvs

    # Create a structured list with different sharing levels
    sharing_levels <- list(
      all_projects = names(asv_counts)[asv_counts == length(asv_list)],
      two_or_more = names(asv_counts)[asv_counts >= 2],
      core_by_project = asv_list,
      asv_project_counts = asv_counts
    )

    return(sharing_levels)
  }
)

# Save the shared ASV data
save(
  core_asvs_by_threshold,
  shared_asvs_by_threshold,
  file = "data/output/rdata/shared_asvs_by_threshold_07.rda"
)

# Print summary statistics
cat("\n", rep("=", 40), "\n")
cat("Shared ASVs Between Projects - Summary\n")
cat(rep("=", 40), "\n\n")

purrr::iwalk(shared_asvs_by_threshold, function(sharing_data, threshold_name) {
  cat("Threshold:", threshold_name, "\n")
  cat(
    "  ASVs above threshold shared across ALL projects:",
    length(sharing_data$all_projects),
    "\n"
  )
  cat(
    "  ASVs above threshold shared by 2+ projects:",
    length(sharing_data$two_or_more),
    "\n"
  )
  cat("  ASVs above threshold by project:\n")
  purrr::iwalk(sharing_data$by_project, function(asvs, proj) {
    cat("    ", proj, ":", length(asvs), "ASVs\n")
  })
  cat("\n")
})
