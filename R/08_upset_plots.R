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

# SECTION 1: Helper functions for UpSet data preparation ----

create_binary_df_from_flat <- function(flat_list) {
  all_items <- unique(unlist(flat_list))
  set_names <- names(flat_list)
  binary_df <- map_dfc(set_names, function(set_name) {
    tibble(!!set_name := as.integer(all_items %in% flat_list[[set_name]]))
  })
  binary_df$item <- all_items
  binary_df |> relocate(item)
}

list_to_binary_df <- function(presence_list) {
  flattened <- purrr::flatten(presence_list)
  create_binary_df_from_flat(flattened)
}

create_upset_df <- function(presence_list, attributes_df, tax_level = "asv") {
  binary_df <- list_to_binary_df(presence_list)

  tax_info <- attributes_df |>
    select(asv, phylum, class, genus) |>
    distinct()

  if (tax_level == "asv") {
    upset_df <- binary_df |>
      rename(asv = item) |>
      left_join(tax_info, by = "asv")
  } else {
    upset_df <- binary_df |>
      rename(!!tax_level := item)
  }

  upset_df
}

create_tax_upset_df <- function(presence_list) {
  binary_df <- list_to_binary_df(presence_list) |>
    rename(taxon = item)
  binary_df
}

create_tax_presence_list <- function(physeq_list, tax_level = "phylum") {
  imap(physeq_list, function(psq, project_name) {
    crops <- sample_data(psq)$crop |> unique()

    crop_taxa <- map(crops, function(crop_name) {
      psq_subset <- prune_samples(sample_data(psq)$crop == crop_name, psq)
      tax_df <- tax_table(psq_subset) |>
        as.data.frame() |>
        rownames_to_column("asv")

      present_asvs <- taxa_names(psq_subset)[taxa_sums(psq_subset) > 0]
      tax_df |>
        filter(asv %in% present_asvs) |>
        pull(!!sym(tax_level)) |>
        unique() |>
        na.omit() |>
        as.character()
    })

    if (grepl("^lamps", project_name, ignore.case = TRUE)) {
      prefix <- gsub("^([^_]+_[^_]+).*", "\\1", project_name)
    } else {
      prefix <- gsub("^([^_]+).*", "\\1", project_name)
    }
    set_names(crop_taxa, paste0(prefix, "_", crops))
  })
}

get_mxg_sets <- function(presence_list) {
  flattened <- purrr::flatten(presence_list)
  mxg_patterns <- c("MXG", "Miscanthus", "_M$")
  mxg_idx <- grep(paste(mxg_patterns, collapse = "|"), names(flattened))
  flattened[mxg_idx]
}

# SECTION 2: Prepare data for all taxonomic levels ----

asv_data <- create_asv_upset_data(main_16S_physeq_list)
phylum_presence <- create_tax_presence_list(main_16S_physeq_list, "phylum")
class_presence <- create_tax_presence_list(main_16S_physeq_list, "class")

# SECTION 3: Summary tables ----

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

# SECTION 4: UpSet DataFrames ----

asv_upset_df <- create_upset_df(
  asv_data$presence_list,
  asv_data$attributes,
  tax_level = "asv"
)

phylum_upset_df <- create_tax_upset_df(phylum_presence) |>
  rename(phylum = taxon)

class_upset_df <- create_tax_upset_df(class_presence) |>
  rename(class = taxon)

# SECTION 5: UpSet plots - All crops ----

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

# SECTION 6: UpSet plots - Miscanthus only ----

mxg_asv_list <- get_mxg_sets(asv_data$presence_list)
mxg_phylum_list <- get_mxg_sets(phylum_presence)
mxg_class_list <- get_mxg_sets(class_presence)

cat("\n", rep("=", 50), "\n")
cat("ASV Intersections - Miscanthus Only\n")
cat(rep("=", 50), "\n")

if (length(mxg_asv_list) >= 2) {
  mxg_asv_df <- create_binary_df_from_flat(mxg_asv_list)
  upset(
    mxg_asv_df,
    intersect = names(mxg_asv_list),
    base_annotations = list(
      "Intersection size" = intersection_size(
        bar_number_threshold = 1,
        fill = "#1f77b4"
      )
    ),
    set_sizes = upset_set_size(geom = geom_bar(fill = "#2ca02c")),
    sort_sets = "descending",
    sort_intersections = "descending"
  ) +
    labs(title = "ASV Intersections - Miscanthus Only") +
    theme(plot.title = element_text(hjust = 0.5))
}

cat("\n", rep("=", 50), "\n")
cat("Phylum Intersections - Miscanthus Only\n")
cat(rep("=", 50), "\n")

if (length(mxg_phylum_list) >= 2) {
  mxg_phylum_df <- create_binary_df_from_flat(mxg_phylum_list)
  upset(
    mxg_phylum_df,
    intersect = names(mxg_phylum_list),
    base_annotations = list(
      "Intersection size" = intersection_size(
        bar_number_threshold = 1,
        fill = "#ff7f0e"
      )
    ),
    set_sizes = upset_set_size(geom = geom_bar(fill = "#d62728")),
    sort_sets = "descending",
    sort_intersections = "descending"
  ) +
    labs(title = "Phylum Intersections - Miscanthus Only") +
    theme(plot.title = element_text(hjust = 0.5))
}

cat("\n", rep("=", 50), "\n")
cat("Class Intersections - Miscanthus Only\n")
cat(rep("=", 50), "\n")

if (length(mxg_class_list) >= 2) {
  mxg_class_df <- create_binary_df_from_flat(mxg_class_list)
  upset(
    mxg_class_df,
    intersect = names(mxg_class_list),
    base_annotations = list(
      "Intersection size" = intersection_size(
        bar_number_threshold = 1,
        fill = "#9467bd"
      )
    ),
    set_sizes = upset_set_size(geom = geom_bar(fill = "#8c564b")),
    sort_sets = "descending",
    sort_intersections = "descending"
  ) +
    labs(title = "Class Intersections - Miscanthus Only") +
    theme(plot.title = element_text(hjust = 0.5))
}

# SECTION 7: Intersection summary tables ----

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

asv_intersection_details <- get_intersection_details(
  asv_upset_df,
  set_cols,
  "asv"
)

shared_asv_summary <- asv_intersection_details |>
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

cat("\n", rep("=", 50), "\n")
cat("ASV Sharing Summary\n")
cat(rep("=", 50), "\n")
print(shared_asv_summary)

core_asvs <- asv_intersection_details |>
  filter(n_sets == max(n_sets)) |>
  left_join(
    asv_data$attributes |> select(asv, phylum, class, genus) |> distinct(),
    by = "asv"
  )

cat("\n", rep("=", 50), "\n")
cat("Core ASVs (Present in Most Crops)\n")
cat(rep("=", 50), "\n")
print(core_asvs)
