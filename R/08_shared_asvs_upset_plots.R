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

cat("\n", rep("=", 40), "\n")
cat("Taxonomic Richness by Crop-Project\n")
cat(rep("=", 40), "\n")
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

cat("\n", rep("=", 40), "\n")
cat("ASV Intersections Across All Crops\n")
cat(rep("=", 40), "\n")

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

# cat("\n", rep("=", 40), "\n")
# cat("Phylum Intersections Across All Crops\n")
# cat(rep("=", 40), "\n")

# phylum_set_cols <- names(purrr::flatten(phylum_presence))
# upset(
#   phylum_upset_df,
#   intersect = phylum_set_cols,
#   n_intersections = 25,
#   base_annotations = list(
#     "Intersection size" = intersection_size(
#       bar_number_threshold = 1,
#       fill = "#ff7f0e"
#     )
#   ),
#   set_sizes = upset_set_size(
#     geom = geom_bar(fill = "#d62728")
#   ),
#   sort_sets = "descending",
#   sort_intersections = "descending"
# ) +
#   labs(title = "Phylum Intersections Across All Crops") +
#   theme(plot.title = element_text(hjust = 0.5))

# cat("\n", rep("=", 40), "\n")
# cat("Class Intersections Across All Crops\n")
# cat(rep("=", 40), "\n")

# class_set_cols <- names(purrr::flatten(class_presence))
# upset(
#   class_upset_df,
#   intersect = class_set_cols,
#   n_intersections = 25,
#   base_annotations = list(
#     "Intersection size" = intersection_size(
#       bar_number_threshold = 1,
#       fill = "#9467bd"
#     )
#   ),
#   set_sizes = upset_set_size(
#     geom = geom_bar(fill = "#8c564b")
#   ),
#   sort_sets = "descending",
#   sort_intersections = "descending"
# ) +
#   labs(title = "Class Intersections Across All Crops") +
#   theme(plot.title = element_text(hjust = 0.5))

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


cat("\n", rep("=", 40), "\n")
cat("ASV Intersections - Miscanthus Only\n")
cat(rep("=", 40), "\n")

# MXG upset
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
    labs(
      title = "ASVs in Top 20 Phyla in Miscanthus",
      subtitle = "16S data only"
    ) +
    theme(
      plot.title = element_text(hjust = 0, size = 14, face = "bold"),
      panel.grid = element_blank(),
      panel.grid.major.x = element_blank()
    )

  mxg_upset
}

## Saving UpSet plot ----
if (!file.exists("data/output/figures/mxg_upset_abundance.svg")) {
  ggsave(
    filename = "data/output/figures/mxg_upset_abundance.svg",
    plot = mxg_upset,
    width = 250,
    height = 225,
    units = "mm",
    dpi = 600
  )
}

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

cat("\n", rep("=", 40), "\n")
cat("ASV Sharing Summary\n")
cat(rep("=", 40), "\n")
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

# SECTION 7: Shared ASVs in Miscanthus Between Projects at Different Thresholds ----

# Define thresholds to test
# These represent the proportion of samples where an ASV must be present
# to be considered "core" (e.g., 0.6 = present in 60% of samples)

# Getting our "core" ASVs in MXG at different thresholds
mxg_prevalence_physeq_16S <- purrr::map(
  main_16S_mxg_physeq_list,
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
project_names <- grep("_DNA$", names(mxg_prevalence_physeq_16S), value = TRUE)

## PROJECT-MXG-THRESHOLD-LEVEL ANALYSIS ----
# Extract core Miscanthus ASVs for each project at each threshold

mxg_prevalent_asvs_by_project_by_threshold <- purrr::map(
  thresholds,
  function(threshold) {
    project_asvs <- purrr::map(project_names, function(proj_name) {
      physeq_object <- mxg_prevalence_physeq_16S

      taxa_names(physeq_object[[proj_name]][[paste0(
        "prevalent_",
        as.character(threshold * 100)
      )]])
    }) %>%
      set_names(gsub("_DNA$", "", project_names))

    return(project_asvs)
  }
) %>%
  set_names(paste0("threshold_", thresholds))


# Calculate shared ASVs between all projects at each threshold (PROJECT-LEVEL)
mxg_shared_asvs_by_project_by_threshold <- purrr::imap(
  mxg_prevalent_asvs_by_project_by_threshold,
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

## PROJECT-CROP-THRESHOLD-LEVEL ANALYSIS ----
# Needs to be calculated with "main_16S_physeq_list" to include all crops

# Extract core ASVs for each crop within each project at each threshold
all_crops_prevalent_asvs_by_threshold <- purrr::map(
  thresholds,
  function(threshold) {
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
  }
) %>%
  set_names(paste0("threshold_", thresholds))


## All crops core ASVs at different thresholds
# Calculate shared ASVs between all crops at each threshold (CROP-LEVEL)
all_crops_shared_asvs_by_threshold <- purrr::imap(
  all_crops_prevalent_asvs_by_threshold,
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
if (!file.exists("data/output/rdata/shared_asvs_by_threshold_08.rda")) {
  save(
    mxg_prevalent_asvs_by_project_by_threshold,
    mxg_shared_asvs_by_project_by_threshold,
    all_crops_prevalent_asvs_by_threshold,
    all_crops_shared_asvs_by_threshold,
    file = "data/output/rdata/shared_asvs_by_threshold_08.rda"
  )
}

# SECTION 8: Summary Statistics of Shared ASVs at Different Thresholds ----
# Print summary statistics - PROJECT-MXG-THRESHOLD LEVEL
cat("\n", rep("=", 40), "\n")
cat("Shared ASVs Between Projects - PROJECT LEVEL Summary\n")
cat(rep("=", 40), "\n\n")

purrr::iwalk(
  mxg_shared_asvs_by_project_by_threshold,
  function(sharing_data, threshold_name) {
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
  }
)

# Print summary statistics - PROJECT-CROP LEVEL
cat("\n", rep("=", 40), "\n")
cat("Shared ASVs Between Crops - CROP LEVEL Summary\n")
cat(rep("=", 40), "\n\n")

purrr::iwalk(
  all_crops_shared_asvs_by_threshold,
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
