###########################################################################
# Abundance Occupancy Analysis on Miscanthus 16S
#
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-24
##########################################################################

source("R/utils/00_setup.R")


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

# Getting our "core" ASVs at different thresholds
main_16S_core_list <- purrr::map(
  main_16S_physeq_list,
  function(physeq_obj) {
    get_prevalent_rare(
      physeq_obj,
      thresholds = c(0, 60, 70, 80, 90),
      detection = 0 / 100,
      include.lowest = FALSE
    )
  }
)

# Modify this vector to analyze different thresholds
thresholds <- c(0.0, 0.6, 0.7, 0.8, 0.9)

# Get project names from the phyloseq list (focusing on DNA projects)
project_names <- grep("_DNA$", names(main_16S_core_list), value = TRUE)

# PROJECT-LEVEL ANALYSIS ----
# Extract core ASVs for each project at each threshold (overall across all crops)
core_asvs_by_threshold <- purrr::map(thresholds, function(threshold) {
  project_asvs <- purrr::map(project_names, function(proj_name) {
    taxa_names(main_16S_core_list[[proj_name]][[paste0(
      "prevalent_",
      as.character(threshold * 100)
    )]])
  }) %>%
    set_names(gsub("_DNA$", "", project_names))

  return(project_asvs)
}) %>%
  set_names(paste0("threshold_", thresholds))


# Calculate shared ASVs between all projects at each threshold (PROJECT-LEVEL)
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

# CROP-LEVEL ANALYSIS ----
# Extract core ASVs for each crop within each project at each threshold
core_asvs_by_crop_by_threshold <- purrr::map(thresholds, function(threshold) {
  # For each project, get ASVs by crop
  all_crop_asvs <- purrr::map(project_names, function(proj_name) {
    # Get abundance/occupancy data by crop
    crop_data <- summarize_abundance_occupancy(main_16S_physeq_list[[
      proj_name
    ]])$by_crop

    # For each crop, get ASVs above threshold
    crop_asvs <- purrr::imap(crop_data, function(crop_df, crop_name) {
      rownames(crop_df %>% filter(occupancy_prop >= threshold))
    })

    # Create proper naming (e.g., ef_MXG, lamps_2018_M)
    proj_prefix <- gsub("_DNA$", "", proj_name)
    if (grepl("^lamps", proj_prefix, ignore.case = TRUE)) {
      # For LAMPS projects, keep the year part
      names(crop_asvs) <- paste0(proj_prefix, "_", names(crop_asvs))
    } else {
      # For other projects (e.g., ef)
      names(crop_asvs) <- paste0(proj_prefix, "_", names(crop_asvs))
    }

    return(crop_asvs)
  })

  # Flatten the list to get all crops across all projects
  all_crop_asvs <- unlist(all_crop_asvs, recursive = FALSE)

  return(all_crop_asvs)
}) %>%
  set_names(paste0("threshold_", thresholds))


# Calculate shared ASVs between all crops at each threshold (CROP-LEVEL)
shared_asvs_by_crop_by_threshold <- purrr::imap(
  core_asvs_by_crop_by_threshold,
  function(asv_list, threshold_name) {
    # Get all unique ASVs across all crops
    all_asvs <- unique(unlist(asv_list))

    # Count how many crops each ASV appears in
    asv_counts <- purrr::map_int(all_asvs, function(asv) {
      sum(purrr::map_lgl(asv_list, ~ asv %in% .x))
    })
    names(asv_counts) <- all_asvs

    # Create a structured list with different sharing levels
    sharing_levels <- list(
      all_crops = names(asv_counts)[asv_counts == length(asv_list)],
      two_or_more = names(asv_counts)[asv_counts >= 2],
      core_by_crop = asv_list,
      asv_crop_counts = asv_counts
    )

    return(sharing_levels)
  }
)

# Save the shared ASV data
save(
  core_asvs_by_threshold,
  shared_asvs_by_threshold,
  core_asvs_by_crop_by_threshold,
  shared_asvs_by_crop_by_threshold,
  file = "data/output/rdata/shared_asvs_by_threshold_07.rda"
)

# Print summary statistics - PROJECT LEVEL
cat("\n", rep("=", 40), "\n")
cat("Shared ASVs Between Projects - PROJECT LEVEL Summary\n")
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
  purrr::iwalk(sharing_data$core_by_project, function(asvs, proj) {
    cat("    ", proj, ":", length(asvs), "ASVs\n")
  })
  cat("\n")
})

# Print summary statistics - CROP LEVEL
cat("\n", rep("=", 60), "\n")
cat("Shared ASVs Between Crops - CROP LEVEL Summary\n")
cat(rep("=", 60), "\n\n")

purrr::iwalk(
  shared_asvs_by_crop_by_threshold,
  function(sharing_data, threshold_name) {
    cat("Threshold:", threshold_name, "\n")
    cat(
      "  ASVs above threshold shared across ALL crops:",
      length(sharing_data$all_crops),
      "\n"
    )
    cat(
      "  ASVs above threshold shared by 2+ crops:",
      length(sharing_data$two_or_more),
      "\n"
    )
    cat("  ASVs above threshold by crop (showing first 10):\n")
    crop_list <- sharing_data$core_by_crop
    if (length(crop_list) > 10) {
      crop_list <- crop_list[1:10]
      show_more <- TRUE
    } else {
      show_more <- FALSE
    }
    purrr::iwalk(crop_list, function(asvs, crop) {
      cat("    ", crop, ":", length(asvs), "ASVs\n")
    })
    if (show_more) {
      cat("    ... and", length(sharing_data$core_by_crop) - 10, "more crops\n")
    }
    cat("\n")
  }
)
