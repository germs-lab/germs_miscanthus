###########################################################################
# Abundance Occupancy Analysis on Miscanthus
#
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-24
# Modified: 2026-01-27
##########################################################################

source("R/utils/00_setup.R")


# Abundance-Occupancy Plots ----

## Quick helper function
multiple_abundance_occupancy_plots <- function(
  physeq_obj,
  physeq_name,
  thresholds = c(0.5, 0.6, 0.7, 0.8, 0.9)
) {
  # Get summaries
  summarized_df_list <- summarize_abundance_occupancy(physeq_obj)

  # Define thresholds to test
  thresholds <- thresholds

  # Plots at different thresholds - overall
  overall_plots <- purrr::map(thresholds, function(threshold) {
    plot_abundance_occupancy(
      summarized_df_list$overall,
      threshold = threshold,
      title = paste0(
        physeq_name,
        " - Overall (threshold:",
        threshold,
        ")"
      )
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

##  16S All crops ----
abundance_occ_plots_16S <- purrr::imap(
  main_16S_physeq_list,
  function(physeq_obj, physeq_name) {
    multiple_abundance_occupancy_plots(physeq_obj, physeq_name)
  }
)

abundance_occ_plots_16S$ef_16S$overall_plots$threshold_0.6

abundance_occ_plots$lamps_2018_16S$overall_plots$threshold_0.6

abundance_occ_plots$lamps_2022_16S$overall_plots$threshold_0.6

## All target regions - All crops ----
abundance_occ_plots_all <- purrr::imap(
  main_physeq_list,
  function(project_list, project_name) {
    purrr::imap(
      project_list,
      function(physeq_obj, physeq_name) {
        multiple_abundance_occupancy_plots(
          physeq_obj,
          physeq_name
        )
      }
    )
  }
)


## Multi-panel Abundance-Occupancy Plots by Target Region ----

# Define target regions and crop mappings
target_regions <- c("16S", "ITS", "AMF")
crop_map <- list(
  ef = "MXG",
  lamps_2018 = "M",
  lamps_2022 = "Miscanthus"
)
# Function to safely extract a plot
safe_extract_plot <- function(
  project_list,
  project_name,
  region,
  crop_name,
  threshold = "threshold_0.6"
) {
  project_list[[paste0(
    project_name,
    "_",
    region,
    "_physeq"
  )]]$by_crop_plots[[crop_name]][[threshold]]
}


# Generate plots for each region
region_plots <- purrr::imap(target_regions, function(region, index) {
  # Extract plots for this region from all three projects
  ef_plot <- safe_extract_plot(
    abundance_occ_plots_all$ef_physeq_list,
    "ef",
    region,
    crop_map$ef
  )

  lamps_2018_plot <- safe_extract_plot(
    abundance_occ_plots_all$lamps_2018_physeq_list,
    "lamps_2018",
    region,
    crop_map$lamps_2018
  )

  lamps_2022_plot <- safe_extract_plot(
    abundance_occ_plots_all$lamps_2022_physeq_list,
    "lamps_2022",
    region,
    crop_map$lamps_2022
  )

  # Collect available plots
  available_plots <- list(
    EF = ef_plot,
    `LAMPS 2018` = lamps_2018_plot,
    `LAMPS 2022` = lamps_2022_plot
  )

  # Remove NULL plots
  available_plots <- purrr::discard(available_plots, is.null)

  # Check if we have at least one plot
  if (length(available_plots) == 0) {
    message("Skipping ", region, " - no plots available across any project")
    return(NULL)
  }

  # Add clean titles to available plots
  plots_with_titles <- purrr::imap(
    available_plots,
    function(plot, project_name) {
      plot +
        ggtitle(project_name) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    }
  )

  # Combine into multi-panel plot
  combined <- patchwork::wrap_plots(
    plots_with_titles,
    ncol = length(plots_with_titles)
  ) +
    patchwork::plot_annotation(
      title = paste0(
        "Abundance-Occupancy Curves: Miscanthus ",
        region,
        " (Threshold 0.6) - ",
        length(plots_with_titles),
        " project(s)"
      ),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
    )

  message(
    "Created plot for ",
    region,
    " with ",
    length(plots_with_titles),
    " project(s)"
  )
  return(combined)
}) |>
  purrr::discard(is.null) %>%
  set_names(paste0(target_regions[!sapply(., is.null)], "_region_plot"))


# Display individual region plots
purrr::walk(region_plots, print)

# Save all plots
purrr::iwalk(region_plots, function(plot, region_name) {
  ggsave(
    filename = paste0(
      "data/output/figures/abundance_occupancy_",
      region_name,
      ".png"
    ),
    plot = plot,
    width = 300,
    height = 250,
    dpi = 300,
    units = "mm"
  )
})
