##########################################################################

# Author: Bolívar Aponte Rolón
# Date: 2025-12-11
##########################################################################

source("R/utils/00_setup.R")

library(UpSetR)

# TODO
# Focus on main_16S_physeq_list then subset to only the MXG within the 16S_DNA

# TODO
# Fix this pipeline

# Create ASV list ----
asv_data <- create_asv_upset_data(main_16S_physeq_list)

# Check lenght of list
length(asv_data$presence_list$ef_16S_DNA$ef_MXG)
length(asv_data$presence_list$lamps_2018_16S_DNA$lamps_2018_M)
# They should all have different number of ASVs

# Table summaries ----

## Quick summary ----
asv_data$presence_list %>%
  flatten() %>%
  map_int(., length) %>%
  enframe(name = "crop_project", value = "n_asvs") %>%
  arrange(desc(n_asvs))

## Phylum summary for each intersection -----
intersection_phylum <- asv_data$attributes |>
  group_by(project, crop) |>
  summarize(
    total_asvs = n_distinct(asv),
    dominant_phylum = names(sort(table(phylum), decreasing = TRUE)[1]),
    n_phyla = n_distinct(phylum),
    n_class = n_distinct(class),

    n_genera = n_distinct(genus),
    .groups = "drop"
  ) |>
  mutate(set_name = paste0(crop))

intersection_phylum

# Plotting matrix with attributes of interest ----
# Create the upset matrix
upset_matrix <- fromList(flatten(asv_data$presence_list))

# Recreate the flattened list to get ASV order matching the matrix rows
flattened_asvs <- flatten(asv_data$presence_list)
all_asvs <- unique(unlist(flattened_asvs))

# Create dataframe linking matrix rows to ASVs
asv_row_map <- tibble(
  matrix_row = seq_along(all_asvs),
  asv = all_asvs
)

# Join with taxonomy info and summarize
upset_df <- asv_row_map |>
  left_join(
    asv_data$attributes |>
      select(asv, phylum, genus) |>
      distinct(),
    by = "asv"
  ) |>
  group_by(matrix_row) |>
  summarize(
    asv = dplyr::first(asv),
    phylum = dplyr::first(phylum),
    genus = dplyr::first(genus),
    n_phyla = n_distinct(phylum),
    n_genera = if_else(all(is.na(genus)), 0L, n_distinct(genus)),
    .groups = "drop",
  ) |>
  bind_cols(as_tibble(upset_matrix, rownames = NA)) |>
  as.data.frame()

upset_df

# Basic UpSet plot: all crop-project combinations

movies <- read.csv(
  system.file("extdata", "movies.csv", package = "UpSetR"),
  header = T,
  sep = ";"
)
str(movies)
str(upset_df)

upset(
  fromList(flatten(asv_data$presence_list)),
  order.by = "freq",
  nsets = 8,
  nintersects = 20,
)

upset(
  upset_df[, 7:14],
  order.by = "freq",
  nsets = 8,
  nintersects = 20,
  boxplot.summary = "n_phyla"
)

ggplot(intersection_phylum, aes(x = crop, y = n_genera)) +
  geom_boxplot() +
  labs(title = "Distribution of Genera per ASV", y = "Number of Genera")
facet_wrap(~crop)


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
  mutate(combo = paste0(crop)) |>
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
