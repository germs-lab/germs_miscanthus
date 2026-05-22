###########################################################################
# Alpha and Beta Diversity Analysis for Nested Phyloseq Objects
#
# This script performs alpha and beta diversity analyses on nested phyloseq
# lists with the structure: main_mxg_physeq_list$project_list$sequencing_type_physeq
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-02
##########################################################################

source("R/utils/00_setup.R")

# Pre-process
main_mxg_physeq_list <- subset_to_miscanthus(main_physeq_list) # Contains  16S, ITS and AMF

# SECTION 1: Alpha Diversity Analysis ----

# Inline: calculate_alpha_diversity — alpha diversity for a flat named list of phyloseq objects
alpha_diversity_results_16S <- purrr::imap(
  main_16S_physeq_list,
  function(physeq_obj, physeq_name) {
    alpha_div_full <- alpha_workflow(physeq_obj = physeq_obj)
    alpha_div_full %>%
      mutate(
        project = gsub("([^_]+).*", "\\1", physeq_name),
        physeq_name = physeq_name,
        sample_date = if ("sample_date" %in% names(.)) {
          as.character(sample_date)
        } else {
          NA_character_
        }
      ) %>%
      relocate(c(observed:inv_simpson), .before = project)
  }
)

# Inline: calculate_alpha_diversity_nested — alpha diversity for a nested list (project > physeq)
main_alpha_diversity_results <- purrr::imap(
  main_physeq_list,
  function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      alpha_div_full <- alpha_workflow(physeq_obj = physeq_obj)
      alpha_div_full %>%
        mutate(
          project = project_name,
          region = gsub(".*_([^_]+)_physeq$", "\\1", physeq_name),
          physeq_name = physeq_name,
          sample_date = if ("sample_date" %in% names(.)) {
            as.character(sample_date)
          } else {
            NA_character_
          }
        ) %>%
        relocate(c(observed:inv_simpson), .before = project)
    })
  }
)

# Flatten lists for easier plotting

alpha_diversity_df_16S <- alpha_diversity_results_16S %>%
  purrr::map_dfr(~.x) %>%
  select(-c(sample_id.y)) %>%
  relocate(sequence_id, .before = sample_id) %>%
  relocate(c(observed:physeq_name), .after = site)


main_alpha_diversity_df_nested <- main_alpha_diversity_results %>%
  purrr::list_flatten(name_spec = "{outer}_{inner}") %>%
  purrr::map_dfr(~.x) %>%
  select(-c(sample_id.y, target_region)) %>%
  relocate(sequence_id, .before = sample_id) %>%
  relocate(c(observed:physeq_name), .after = site)

# View summary
summary(alpha_diversity_df_16S[, c(
  "observed",
  "shannon",
  "simpson",
  "inv_simpson"
)])


# SECTION 2: Alpha Diversity Plots ----

## Define comparison mappings ----
comparison_map <- list(
  crop = list(
    lamps_2018 = list(c("C", "M")),
    lamps_2022 = list(
      c("Miscanthus", "Grass"),
      c("Miscanthus", "Maize"),
      c("Maize", "Grass")
    ),
    ef = list(c("MXG", "SB"), c("MXG", "ZM"), c("ZM", "SB"))
  )
)

## Helpers ----
get_comparisons <- function(data, physeq_name, comp_map) {
  # Extract project/dataset identifier from physeq_name

  if (grepl("^lamps", physeq_name, ignore.case = TRUE)) {
    # LAMPS:  extract first two parts (lamps_2018 or lamps_2022)
    dataset_key <- gsub("^([^_]+_[^_]+).*", "\\1", physeq_name)
  } else {
    # (ef, mxg, etc. )
    dataset_key <- gsub("^([^_]+).*", "\\1", physeq_name)
  }

  return(comp_map$crop[[dataset_key]])
}

alpha_plot_workflow <- function(alpha_data, physeq_name, comp_map) {
  # Get the appropriate comparisons
  comparisons <- get_comparisons(
    alpha_data,
    physeq_name,
    comp_map = comparison_map
  )

  # Determine grouping variable
  group_var <- names(alpha_data)[
    names(alpha_data) %in% names(comparison_map)
  ][1]

  comparison_setting <- {
    if (!is.null(comparisons)) {
      ggpubr::stat_compare_means(
        comparisons = comparisons,
        method = "t.test",
        label = "p.signif"
        # symnum.args = list(
        #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
        #   symbols = c("****", "***", "**", "*", "ns")
        # )
      )
    }
  }

  theme_settings <- list(
    theme_minimal(),
    theme(
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  )

  # Observed richness plot
  p_observed <- create_alpha_plot(
    alpha_data,
    group_var = "crop",
    "observed",
    title = "Observed Richness",
    subtitle = physeq_name,
    y_title_add = "index"
  ) +
    comparison_setting +
    theme_settings

  # Shannon diversity plot
  p_shannon <- create_alpha_plot(
    alpha_data,
    group_var = "crop",
    "shannon",
    title = "Shannon Diversity",
    subtitle = physeq_name,
    y_title_add = "index"
  ) +
    comparison_setting +
    theme_settings

  # Simpson diversity plot
  p_simpson <- create_alpha_plot(
    alpha_data,
    group_var = "crop",
    "simpson",
    title = "Simpson Diversity",
    subtitle = physeq_name,
    y_title_add = "index"
  ) +
    comparison_setting +
    theme_settings

  all_div_plots <- ggpubr::ggarrange(
    p_observed,
    p_shannon,
    p_simpson,
    labels = c("A", "B", "C"),
    nrow = 1,
    ncol = 3
  )

  # Return list of plots
  list(
    all_div_plots = all_div_plots
  )
}


## Alpha diversity plots for 16S ----
alpha_diversity_plots_16S <-
  purrr::imap(alpha_diversity_results_16S, function(alpha_data, physeq_name) {
    alpha_plot_workflow(
      alpha_data = alpha_data,
      physeq_name = physeq_name,
      comp_map = comparison_map
    )
  })


## Alpha plots for main_physeq ----
main_alpha_diversity_plots <- purrr::imap(
  main_alpha_diversity_results,
  function(project_list, project_name) {
    purrr::imap(project_list, function(alpha_data, physeq_name) {
      alpha_plot_workflow(
        alpha_data = alpha_data,
        physeq_name = physeq_name,
        comp_map = comparison_map
      )
    })
  }
)

cat("\n", rep("=", 40), "\n")
cat("Alpha Diversity by Project and Crop\n")
cat(rep("=", 40), "\n")
print(alpha_diversity_plots_16S$ef_16S_DNA$all_div_plots)
print(main_alpha_diversity_plots$ef_physeq_list$ef_16S_physeq$all_div_plots) # These shouldn't be too different

print(alpha_diversity_plots_16S)
print(main_alpha_diversity_plots)


# SECTION 3: Beta Diversity Analysis ----

# # Calculate beta diversity (Bray-Curtis dissimilarity)
beta_diversity_results_16S <- calculate_beta_diversity(
  main_16S_physeq_list,
  method = "PCoA",
  distance = "bray"
)

mxg_beta_diversity_results <- calculate_beta_diversity(
  main_mxg_physeq_list,
  method = "PCoA",
  distance = "bray",
  nested = TRUE
)


combined_df <- purrr::imap_dfr(
  phyloseq_list,
  function(project_list, project_name) {
    purrr::imap_dfr(
      project_list,
      function(physeq_obj, physeq_name) {
        physeq_filtered <- physeq_obj %>%
          prune_samples(sample_sums(.) > 0, .) %>%
          prune_taxa(taxa_sums(.) > 0, .) %>%
          psmelt() %>%
          select(OTU, Sample, Abundance) |>
          mutate(project = project_name, nucleotide = )
      }
    )
  }
)
test <- phyloseq::psmelt(main_mxg_physeq_list$ef$ef_16S_DNA) |>
  select(OTU, Sample, Abundance) |>
  mutate(project = "ef", nucleotide = "16S")
# group_by(project, crop) %>%
# summarise(n = n(), .groups = "drop")

main_beta_diversity_results <- calculate_beta_diversity(
  main_physeq_list,
  method = "PCoA",
  distance = "bray",
  nested = TRUE
)

# SECTION 5: Alpha Diversity Summary Statistics

# Summary table of alpha diversity by project and crop
alpha_summary <- alpha_diversity_df_16S %>%
  group_by(project, crop) %>%
  summarise(
    n = n(),
    mean_observed = mean(observed, na.rm = TRUE),
    sd_observed = sd(observed, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    mean_simpson = mean(simpson, na.rm = TRUE),
    sd_simpson = sd(simpson, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n", rep("=", 40), "\n")
cat("Alpha Diversity Summary by Project and Crop\n")
cat(rep("=", 40), "\n")
print(alpha_summary)

# SECTION 5: Beta Diversity Plots (PCoA)
## Beta diversity plots for each phyloseq object
beta_diversity_plots_16S <-
  purrr::imap(beta_diversity_results_16S, function(beta_data, physeq_name) {
    beta_plot_workflow(beta_data = beta_data, physeq_name = physeq_name)
  })

mxg_beta_diversity_plots <- purrr::imap(
  mxg_beta_diversity_results,
  function(project_list, project_name) {
    purrr::imap(project_list, function(beta_data, physeq_name) {
      beta_plot_workflow(
        beta_data = beta_data,
        physeq_name = physeq_name,
        color = "crop"
      )
    })
  }
)

combined_df <- purrr::imap_dfr(
  main_mxg_physeq_list,
  function(project_list, project_name) {
    purrr::imap_dfr(
      project_list,
      function(physeq_obj, physeq_name) {
        physeq_filtered <- physeq_obj %>%
          prune_samples(sample_sums(.) > 0, .) %>%
          prune_taxa(taxa_sums(.) > 0, .) %>%
          psmelt() %>%
          mutate(project = project_name, physeq_name = physeq_name)
      }
    )
  }
)

main_beta_diversity_plots <-
  purrr::imap(
    main_beta_diversity_results,
    function(project_list, project_name) {
      purrr::imap(project_list, function(beta_data, physeq_name) {
        beta_plot_workflow(beta_data = beta_data, physeq_name = physeq_name)
      })
    }
  )

# Beta diversity plots are stored in the nested list

cat("\n", rep("=", 40), "\n")
cat("Beta Diversity Project and Crop\n")
cat(rep("=", 40), "\n")
print(beta_diversity_plots_16S)
print(mxg_beta_diversity_plots$ef)

# Checkpoint
alpha_beta_plots <- list(
  alpha_diversity_plots = alpha_diversity_plots_16S,
  beta_diversity_plots = beta_diversity_plots_16S
)
# save(
#   alpha_beta_plots,
#   file = "data/output/rdata/alpha_beta_plots.rda"
# )
