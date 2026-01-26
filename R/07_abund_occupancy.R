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

# TODO
# Polish in to a cohesive workflow for all 16S phyloseqs
# 2025-11-24

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


test <- summarize_abundance_occupancy(
  main_16S_physeq_list$ef_16S
)
test$overall

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
