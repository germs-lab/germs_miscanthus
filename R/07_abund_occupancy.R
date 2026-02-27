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

# Display plots
purrr::walk(abundance_occ_plots_all, function(project_plots) {
  purrr::walk(project_plots, function(physeq_plots) {
    print(physeq_plots$overall_plots$threshold_0.6)
    # Filter for Miscanthus crops only
    miscanthus_crops <- c("MXG", "M", "Miscanthus")
    miscanthus_plots <- physeq_plots$by_crop_plots[
      names(physeq_plots$by_crop_plots) %in% miscanthus_crops
    ]

    purrr::walk(miscanthus_plots, function(crop_plots) {
      print(crop_plots$threshold_0.6)
    })
  })
})

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
        "- ",
        length(plots_with_titles),
        " project(s)"
      ),
      subtitle = "Occupancy threshold: 0.6",
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

# SECTION 2: Abundance-Occupancy by time and location ----

# EF we have "timing" which is a proxy for sampling_date as the main varibale for time.

# LAMPS 2018 we have "year" which has 2016-2018, but this seems to be time of planting because "sampling_date" has the same date `20180425`.

# LAMPS 2022 has "site" which has `Central`, `Southeast` and `Northwest` and only one year of sampling (2022) so we can't do a time series but we can do a location comparison.

# Subset helper function
subset_physeq_by_variable <- function(physeq_obj, variable, levels) {
  keep_df <- sample_data(physeq_obj)[[variable]] == levels

  ps_subset <- prune_samples(keep_df, physeq_obj) |>
    filter_taxa(function(x) sum(x) > 0, TRUE)

  summary_df <- summarize_abundance_occupancy(ps_subset)
  ugly_names <- deparse(substitute(physeq_obj))

  if (grepl("^lamps", ugly_names, ignore.case = TRUE)) {
    # LAMPS:  extract first two parts (lamps_2018 or lamps_2022)
    pretty_name <- gsub("^([^_]+_[^_]+).*", "\\1", ugly_names)
  } else {
    # (ef, mxg, etc. )
    pretty_name <- gsub("^([^_]+).*", "\\1", ugly_names)
  }

  plot_abundance_occupancy(
    summary_df$overall,
    threshold = 0.6,
    title = paste(
      stringr::str_to_upper(pretty_name),
      "Miscanthus",
      "-",
      levels,
      sep = " "
    )
  )
}


# EF: Subset by "timing" (temporal proxy)
ef_mxg_by_timing <- main_16S_mxg_physeq_list$ef_16S |>
  subset_samples(crop == "MXG")

timing_levels <- unique(sample_data(ef_mxg_by_timing)$timing)

ef_timing_plots <- purrr::map(timing_levels, function(time_point) {
  subset_physeq_by_variable(ef_mxg_by_timing, "timing", time_point)
}) |>
  set_names(paste0("timing_", timing_levels))

# LAMPS 2018: Subset by "year" (but this is likely time of planting, not sampling)
lamps_2018_mxg_by_year <- main_16S_mxg_physeq_list$lamps_2018_16S |>
  subset_samples(crop == "M")
year_levels <- unique(sample_data(lamps_2018_mxg_by_year)$year)

lamps_2018_year_plots <- purrr::map(year_levels, function(year) {
  subset_physeq_by_variable(lamps_2018_mxg_by_year, "year", year)
}) |>
  set_names(paste0("year_", year_levels))

# LAMPS 2022: Subset by "site"
lamps_2022_mxg_by_site <- main_16S_mxg_physeq_list$lamps_2022_16S |>
  subset_samples(crop == "Miscanthus")
site_levels <- unique(sample_data(lamps_2022_mxg_by_site)$site)

lamps_2022_site_plots <- purrr::map(site_levels, function(site) {
  subset_physeq_by_variable(lamps_2022_mxg_by_site, "site", site)
}) |>
  set_names(paste0("site_", site_levels))


# Arrange plots into multi-panels
timing_panel <- patchwork::wrap_plots(
  ef_timing_plots,
  ncol = length(ef_timing_plots)
) +
  patchwork::plot_annotation(
    title = "Abundance-Occupancy Curves: Miscanthus 16S - EF by Timing",
    subtitle = "Occupancy threshold: 0.6",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  )
year_panel <- patchwork::wrap_plots(
  lamps_2018_year_plots,
  ncol = length(lamps_2018_year_plots)
) +
  patchwork::plot_annotation(
    title = "Abundance-Occupancy Curves: Miscanthus 16S - LAMPS 2018 by Year",
    subtitle = "Occupancy threshold: 0.6",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  )
site_panel <- patchwork::wrap_plots(
  lamps_2022_site_plots,
  ncol = length(lamps_2022_site_plots)
) +
  patchwork::plot_annotation(
    title = "Abundance-Occupancy Curves: Miscanthus 16S - LAMPS 2022 by Site",
    subtitle = "Occupancy threshold: 0.6",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  )


# Display panels
panel_list <- list(
  EF = timing_panel,
  LAMPS_2018 = year_panel,
  LAMPS_2022 = site_panel
)

purrr::walk(panel_list, print)

# Save panels
purrr::iwalk(
  panel_list,
  function(panel, project_name) {
    panel_type <- case_when(
      identical(panel, timing_panel) ~ "timing_panel",
      identical(panel, year_panel) ~ "year_panel",
      identical(panel, site_panel) ~ "site_panel"
    )

    ggsave(
      filename = paste0(
        "data/output/figures/abundance_occupancy_",
        project_name,
        "_",
        panel_type,
        ".png"
      ),
      width = 300,
      height = 250,
      dpi = 300,
      units = "mm"
    )
  }
)
