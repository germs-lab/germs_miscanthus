##########################################################################

# Author: Bolívar Aponte Rolón
# Date: 2025-12-11
##########################################################################

source("R/utils/00_setup.R")


# Explore beta_diversity_results structure
str(main_16S_physeq_list, max.level = 2)
library(UpSetR)
movies <- read.csv(
  system.file("extdata", "movies.csv", package = "UpSetR"),
  header = T,
  sep = ";"
)

listInput <- list(
  EF = rownames(otu_table(main_16S_physeq_list$ef_16S_DNA)),
  LAMPS_2018 = rownames(otu_table(main_16S_physeq_list$lamps_2018_16S_DNA))
)


upset(fromList(listInput), order.by = "freq")

# TODO
# Focus on main_16S_physeq_list then subset to only the MXG within the 16S_DNA

# Check sample data columns for each project
main_16S_physeq_list |>
  map(~ sample_data(.x) |> colnames()) |>
  enframe(name = "project", value = "columns")


# TODO
# Fix this pipeline
library(UpSetR)
# Pipeline 2: Enhanced ASV presence with biological attributes

# Create the list
asv_data <- create_asv_upset_data(main_16S_physeq_list)

# Check lenght of list
length(asv_data$presence_list$ef_16S_DNA$ef_MXG)

# Quick summary
map_int(asv_data$presence_list, length) |>
  enframe(name = "crop_project", value = "n_asvs")


# Basic UpSet plot: all crop-project combinations
# Create a matrix suitable for boxplot.summary
phylum_matrix <- asv_data$attributes |>
  select(asv, project, crop, phylum) |>
  distinct() |>
  mutate(set_name = paste0(substr(project, 1, 4), "_", crop)) |>
  pivot_wider(
    names_from = set_name,
    values_from = phylum,
    values_fill = NA
  ) |>
  #column_to_rownames("asv") |>
  as.matrix()
# Create a phylum summary for each intersection
intersection_phylum <- asv_data$attributes |>
  group_by(project, crop) |>
  summarize(
    total_asvs = n_distinct(asv),
    dominant_phylum = names(sort(table(phylum), decreasing = TRUE)[1]),
    n_phyla = n_distinct(phylum),
    .groups = "drop"
  ) |>
  mutate(set_name = paste0(substr(project, 1, 4), "_", crop))

intersection_phylum

upset(
  fromList(flatten(asv_data$presence_list)),
  order.by = "freq",
  nsets = 8,
  nintersects = 20,
  boxplot.summary = phylum_matrix
)


upset(
  fromList(flatten(asv_data$presence_list)),
  order.by = "freq",
  nsets = 8,
  nintersects = 20,
  boxplot.summary = top(asv_data$attributes, "phylum", 5)
)


# SCRAPS ----

# View attributes

# Summary by phylum across projects
asv_data$attributes |>
  group_by(phylum, project, crop) |>
  summarize(n_asv = n_distinct(asv), .groups = "drop") |>
  arrange(phylum, project)


# Create attribute summary: for each intersection, show phylum composition
# Since UpSetR attribute.plots is complex, let's create a detailed summary table instead

# Get all intersections and their phylum composition
intersection_summary <- asv_data$attributes |>
  group_by(project, crop) |>
  summarize(
    total_asvs = n_distinct(asv),
    n_phyla = n_distinct(phylum),
    dominant_phylum = names(sort(table(phylum), decreasing = TRUE)[1]),
    top_2_phyla = paste(
      names(sort(table(phylum), decreasing = TRUE)[1:2]),
      collapse = ", "
    ),
    .groups = "drop"
  ) |>
  mutate(combo = paste0(substr(project, 1, 4), "_", crop)) |>
  arrange(desc(total_asvs))

intersection_summary

# Count shared ASVs between crops
shared_summary <- asv_data$attributes |>
  group_by(asv) |>
  summarize(
    n_crops = n_distinct(crop),
    crops = paste(sort(unique(crop)), collapse = ", "),
    phylum = dplyr::first(phylum),
    .groups = "drop"
  ) |>
  group_by(n_crops) |>
  summarize(
    n_asvs = n(),
    .groups = "drop"
  ) |>
  mutate(
    interpretation = case_when(
      n_crops == 1 ~ "Crop-specific ASVs",
      n_crops == 2 ~ "Shared between 2 crops",
      n_crops == 3 ~ "Shared between 3 crops",
      TRUE ~ "Shared across multiple crops"
    )
  )

shared_summary
