##########################################################################
# UpSet Plots: Exploring ASV, Phylum, and Class Intersections
#
# Explores shared microbial taxa across Miscanthus experiments using UpSet
# plots to visualize intersections at ASV, Phylum, and Class levels.
#
# Author: Bolívar Aponte Rolón
# Date: 2025-12-11
##########################################################################

source("R/utils/00_setup.R")
library(ComplexUpset)

# SECTION 1: Prepare data for all taxonomic levels ----

asv_data <- create_asv_upset_data(main_16S_physeq_list)
phylum_presence <- create_tax_presence_list(main_16S_physeq_list, "phylum")
class_presence <- create_tax_presence_list(main_16S_physeq_list, "class")

# SECTION 2: Summary tables ----

asv_summary <- asv_data$presence_list |>
  purrr::flatten() |>
  map_int(length) |>
  enframe(name = "crop_project", value = "n_asvs") |>
  arrange(desc(n_asvs))

phylum_summary <- phylum_presence |>
  purrr::flatten() |>
  map_int(length) |>
  enframe(name = "crop_project", value = "n_phyla") |>
  arrange(desc(n_phyla))

class_summary <- class_presence |>
  purrr::flatten() |>
  map_int(length) |>
  enframe(name = "crop_project", value = "n_classes") |>
  arrange(desc(n_classes))

combined_summary <- asv_summary |>
  left_join(phylum_summary, by = "crop_project") |>
  left_join(class_summary, by = "crop_project")

cat("\n", rep("=", 50), "\n")
cat("Taxonomic Richness by Crop-Project\n")
cat(rep("=", 50), "\n")
print(combined_summary)

# SECTION 3: UpSet DataFrames ----

asv_upset_df <- create_upset_df(
  asv_data$presence_list,
  asv_data$attributes,
  tax_level = "asv"
)

phylum_upset_df <- create_tax_upset_df(phylum_presence) |>
  rename(phylum = taxon)

class_upset_df <- create_tax_upset_df(class_presence) |>
  rename(class = taxon)

# SECTION 4: UpSet plots - All crops ----

cat("\n", rep("=", 50), "\n")
cat("ASV Intersections Across All Crops\n")
cat(rep("=", 50), "\n")

set_cols <- names(purrr::flatten(asv_data$presence_list))
upset(
  asv_upset_df,
  intersect = set_cols,
  n_intersections = 20,
  base_annotations = list(
    "Intersection size" = intersection_size(
      bar_number_threshold = 1,
      fill = "#1f77b4",
      color = "#1f77b4"
    )
  ),
  queries = list(
    upset_query(
      intersect = set_cols,
      color = 'red',
      fill = 'red',
      only_components = c('intersections_matrix', 'Intersection size')
    )
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(fill = "#2ca02c")
  ),
  sort_sets = "descending",
  sort_intersections = "descending",
  themes =
) +
  labs(title = "ASV Intersections Across All Crops") +
  theme(plot.title = element_text(hjust = 0.5))

cat("\n", rep("=", 50), "\n")
cat("Phylum Intersections Across All Crops\n")
cat(rep("=", 50), "\n")

phylum_set_cols <- names(purrr::flatten(phylum_presence))
upset(
  phylum_upset_df,
  intersect = phylum_set_cols,
  n_intersections = 25,
  base_annotations = list(
    "Intersection size" = intersection_size(
      bar_number_threshold = 1,
      fill = "#ff7f0e"
    )
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(fill = "#d62728")
  ),
  sort_sets = "descending",
  sort_intersections = "descending"
) +
  labs(title = "Phylum Intersections Across All Crops") +
  theme(plot.title = element_text(hjust = 0.5))

cat("\n", rep("=", 50), "\n")
cat("Class Intersections Across All Crops\n")
cat(rep("=", 50), "\n")

class_set_cols <- names(purrr::flatten(class_presence))
upset(
  class_upset_df,
  intersect = class_set_cols,
  n_intersections = 25,
  base_annotations = list(
    "Intersection size" = intersection_size(
      bar_number_threshold = 1,
      fill = "#9467bd"
    )
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(fill = "#8c564b")
  ),
  sort_sets = "descending",
  sort_intersections = "descending"
) +
  labs(title = "Class Intersections Across All Crops") +
  theme(plot.title = element_text(hjust = 0.5))

# SECTION 5: UpSet plots - Miscanthus only ----

mxg_asv_list <- get_mxg_sets(asv_data$presence_list)
mxg_phylum_list <- get_mxg_sets(phylum_presence)
mxg_class_list <- get_mxg_sets(class_presence)

label_crop_project <- function(x) {
  case_when(
    grepl("^ef_", x) ~ "Energy Farm",
    grepl("^lamps_2018_", x) ~ "LAMPS 2018",
    grepl("^lamps_2022_", x) ~ "LAMPS 2022",
    TRUE ~ x
  )
}


cat("\n", rep("=", 50), "\n")
cat("ASV Intersections - Miscanthus Only\n")
cat(rep("=", 50), "\n")


if (length(mxg_asv_list) >= 2) {
  mxg_asv_df <- create_binary_df_from_flat(mxg_asv_list) |>
    rename(asv = item)
  mxg_phylum_df <- create_binary_df_from_flat(mxg_phylum_list)

  # Get phylum for each ASV from the taxonomy table
  asv_phylum_map <- asv_data$attributes |>
    select(asv, phylum) |>
    distinct() |>
    filter(!is.na(phylum))

  top_20_phyla <- asv_phylum_map |>
    count(phylum, name = "n_asvs") |>
    slice_max(order_by = n_asvs, n = 20, with_ties = FALSE)

  # Join phylum to ASV binary df
  mxg_asv_df_annotated <- mxg_asv_df |>
    left_join(
      asv_data$attributes |> select(asv, phylum) |> distinct(),
      by = "asv"
    ) |>
    filter(phylum %in% top_20_phyla$phylum)

  # Miscanthus UpSet plot with phylum annotation
  mxg_upset <- upset(
    mxg_asv_df_annotated,
    intersect = names(mxg_asv_list),
    base_annotations = list(
      "Intersection size" = intersection_size(
        bar_number_threshold = 1,
        #fill = "#1f77b4",
        text = list(
          vjust = -0.5
        )
      )
    ),
    annotations = list(
      "Phylum" = (ggplot(mapping = aes(fill = phylum)) +
        geom_bar(stat = "count", position = "fill") +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_viridis_d(option = "turbo", na.value = "grey50") +
        labs(y = "Relative Abundance", fill = "Top 20 Phyla") +
        theme(
          legend.position = "right",
          legend.title = element_text(face = "bold")
        ))
    ),
    queries = list(
      upset_query(
        intersect = c('ef_MXG', 'lamps_2018_M', 'lamps_2022_Miscanthus'),
        color = '#47c204ff',
        fill = '#47c204ff',
        only_components = c('intersections_matrix', 'Intersection size')
      )
    ),
    width_ratio = 0.1,
    min_size = 10,
    wrap = TRUE,
    labeller = label_crop_project,
    set_sizes = FALSE,
    guides = 'over'
  ) +
    labs(title = "ASVs in Top 20 Phyla in Miscanthus", ) +
    theme(
      plot.title = element_text(hjust = 0, size = 14, face = "bold"),
      panel.grid = element_blank(),
      panel.grid.major.x = element_blank()
    )

  mxg_upset
}

## Saving UpSet plot ----
ggsave(
  filename = "data/output/figures/mxg_upset_abundance.svg",
  plot = mxg_upset,
  width = 250,
  height = 225,
  units = "mm",
  dpi = 600
)

# SECTION 6: MXG Intersection summary tables ----

get_intersection_details <- function(upset_df, set_cols, id_col = "asv") {
  upset_df |>
    mutate(
      n_sets = rowSums(across(all_of(set_cols))),
      sets_present = apply(
        across(all_of(set_cols)),
        1,
        function(x) paste(set_cols[x == 1], collapse = ", ")
      )
    ) |>
    select(all_of(id_col), n_sets, sets_present) |>
    arrange(desc(n_sets))
}

# Get MXG ASV intersection details
mxg_intersection_details <- get_intersection_details(
  mxg_asv_df_annotated,
  names(mxg_asv_df_annotated)[-c(1, 5)],
  "asv"
)

shared_mxg_asv_summary <- mxg_intersection_details |>
  group_by(n_sets) |>
  summarize(count = n(), .groups = "drop") |>
  mutate(
    category = case_when(
      n_sets == 1 ~ "Unique to one crop",
      n_sets == 2 ~ "Shared by 2 crops",
      n_sets == 3 ~ "Shared by 3 crops",
      TRUE ~ paste0("Shared by ", n_sets, " crops")
    )
  )

# Details of shared ASVs (shared by all 3 crops)
mxg_shared_asvs <- mxg_intersection_details |>
  filter(n_sets == max(n_sets)) |>
  left_join(
    asv_data$attributes |> select(asv, phylum, class, genus) |> distinct(),
    by = "asv"
  )

count_shared_asvs <- mxg_shared_asvs |>
  count(phylum, sort = TRUE) |>
  mutate(pct = n / sum(n) * 100)

# Compare unique vs shared phylum composition
compare_shared_unique <- mxg_asv_df_annotated |>
  mutate(
    n_projects = ef_MXG + lamps_2018_M + lamps_2022_Miscanthus,
    sharing = case_when(
      n_projects == 3 ~ "Core (all 3)",
      n_projects == 2 ~ "Shared (2)",
      TRUE ~ "Unique (1)"
    )
  ) |>
  count(sharing, phylum) |>
  pivot_wider(names_from = sharing, values_from = n, values_fill = 0)

cat("\n", rep("=", 50), "\n")
cat("ASV Sharing Summary\n")
cat(rep("=", 50), "\n")
print(shared_mxg_asv_summary)

cat("\n", rep("=", 40), "\n")
cat("Shared ASVs in Miscanthus\n")
cat(rep("=", 40), "\n")
print(count_shared_asvs)

cat("\n", rep("=", 40), "\n")
cat("Comparison of Shared vs Unique ASVs by Phylum:\n")
print(compare_shared_unique)

cat("\n", rep("=", 40), "\n")
cat("Details of Shared ASVs:\n")
print(mxg_shared_asvs)

# SECTION 7: Shared ASVs by Threshold - UpSet Plots ----

# Load shared ASV data from 07_abund_occupancy.R
load("data/output/rdata/shared_asvs_by_threshold_07.rda")

# Ensure asv_data is available (should be loaded from SECTION 1)
if (!exists("asv_data")) {
  asv_data <- create_asv_upset_data(main_16S_physeq_list)
}

## VALIDATION: Check that threshold_0.0 matches create_asv_upset_data() ----
cat("\n", rep("=", 60), "\n")
cat("VALIDATION: Comparing threshold_0.0 with create_asv_upset_data()\n")
cat(rep("=", 60), "\n\n")

# Get the crop-level ASVs from create_asv_upset_data (flattened)
upset_data_crop_asvs <- purrr::flatten(asv_data$presence_list)

# Get the crop-level ASVs from threshold_0.0
threshold_0_crop_asvs <- shared_asvs_by_crop_by_threshold$threshold_0.0$core_by_crop

# Compare counts
cat("Crop-level comparison:\n")
all_crop_names <- union(
  names(upset_data_crop_asvs),
  names(threshold_0_crop_asvs)
)

validation_results <- purrr::map_dfr(all_crop_names, function(crop_name) {
  upset_count <- length(upset_data_crop_asvs[[crop_name]])
  threshold_count <- length(threshold_0_crop_asvs[[crop_name]])
  match <- upset_count == threshold_count

  tibble(
    crop = crop_name,
    upset_data_count = upset_count,
    threshold_0_count = threshold_count,
    match = match
  )
})

print(validation_results)

# Overall summary
all_match <- all(validation_results$match, na.rm = TRUE)
if (all_match) {
  cat("\n✓ VALIDATION PASSED: All crop-level counts match!\n")
} else {
  cat("\n✗ VALIDATION FAILED: Some crop-level counts do not match.\n")
  cat("Mismatches:\n")
  print(validation_results %>% filter(!match))
}

# Project-level validation
cat("\n\nProject-level comparison:\n")
cat("Checking if project-level threshold_0.0 equals union of crop ASVs:\n\n")

project_validation <- purrr::map_dfr(
  names(shared_asvs_by_threshold$threshold_0.0$core_by_project),
  function(proj_name) {
    # Get project-level ASVs at threshold 0.0
    project_asvs <- shared_asvs_by_threshold$threshold_0.0$core_by_project[[
      proj_name
    ]]

    # Get corresponding crop-level ASVs and take their union
    # Match crops to this project (e.g., ef_MXG, ef_SB for ef_16S)
    crop_pattern <- paste0("^", proj_name, "_")
    matching_crops <- grep(
      crop_pattern,
      names(threshold_0_crop_asvs),
      value = TRUE
    )

    if (length(matching_crops) > 0) {
      crop_union <- unique(unlist(threshold_0_crop_asvs[matching_crops]))

      # Compare
      project_count <- length(project_asvs)
      union_count <- length(crop_union)
      match <- project_count == union_count

      # Also check if the actual ASVs are the same
      same_asvs <- setequal(project_asvs, crop_union)

      tibble(
        project = proj_name,
        project_level_count = project_count,
        crop_union_count = union_count,
        counts_match = match,
        asvs_identical = same_asvs
      )
    } else {
      tibble(
        project = proj_name,
        project_level_count = length(project_asvs),
        crop_union_count = NA,
        counts_match = NA,
        asvs_identical = NA
      )
    }
  }
)

print(project_validation)

if (all(project_validation$asvs_identical, na.rm = TRUE)) {
  cat(
    "\n✓ VALIDATION PASSED: Project-level ASVs equal union of crop-level ASVs!\n"
  )
} else {
  cat(
    "\n✗ VALIDATION FAILED: Project-level ASVs do not equal union of crop-level ASVs.\n"
  )
}

cat("\n", rep("=", 60), "\n\n")

cat("\n", rep("=", 40), "\n")
cat("UpSet Plots for Shared ASVs at Different Thresholds\n")
cat(rep("=", 40), "\n\n")

# Create UpSet plots for each threshold
threshold_upset_plots <- purrr::imap(
  shared_asvs_by_threshold,
  function(sharing_data, threshold_name) {
    # Get the project-level ASV lists
    asv_list <- sharing_data$core_by_project

    # Create binary dataframe for upset plot
    upset_df <- create_binary_df_from_flat(asv_list) |>
      rename(asv = item)

    # Add taxonomic information
    asv_tax_info <- asv_data$attributes |>
      select(asv, phylum, class, genus) |>
      distinct()

    upset_df_annotated <- upset_df |>
      left_join(asv_tax_info, by = "asv")

    # Create UpSet plot
    set_cols <- names(asv_list)
    n_projects <- length(set_cols)

    plot_obj <- upset(
      upset_df_annotated,
      intersect = set_cols,
      n_intersections = 20,
      base_annotations = list(
        "Intersection size" = intersection_size(
          bar_number_threshold = 1,
          fill = "#1f77b4",
          text = list(vjust = -0.5)
        )
      ),
      queries = list(
        upset_query(
          intersect = set_cols,
          color = 'darkgreen',
          fill = 'darkgreen',
          only_components = c('intersections_matrix', 'Intersection size')
        )
      ),
      set_sizes = upset_set_size(
        geom = geom_bar(fill = "#2ca02c")
      ),
      sort_sets = "descending",
      sort_intersections = "descending",
      width_ratio = 0.15
    ) +
      labs(
        title = paste0("Shared ASVs Between Projects (", threshold_name, ")")
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid = element_blank()
      )

    # Save plot
    ggsave(
      filename = paste0(
        "data/output/figures/shared_asvs_",
        threshold_name,
        "_upset.svg"
      ),
      plot = plot_obj,
      width = 250,
      height = 200,
      units = "mm",
      dpi = 600
    )

    return(plot_obj)
  }
)

# Display the plots
purrr::walk(threshold_upset_plots, print)

cat("\n", rep("=", 60), "\n")
cat("UpSet plots saved to data/output/figures/\n")
cat(rep("=", 60), "\n")
