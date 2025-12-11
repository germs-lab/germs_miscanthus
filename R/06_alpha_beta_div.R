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

#--------------------------------------------------------
# SECTION 1: Alpha Diversity Analysis ----
#--------------------------------------------------------

# Calculate alpha diversity
alpha_diversity_results_16S <- calculate_alpha_diversity(main_16S_physeq_list)

main_alpha_diversity_results <- calculate_alpha_diversity_nested(
  main_physeq_list
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

#--------------------------------------------------------
# SECTION 2: Alpha Diversity Plots ----
#--------------------------------------------------------

# Define comparison mappings
comparison_map <- list(
  crop = list(
    lamps_2018_physeq = list(treatment = list(c("C", "M"))),
    lamps_2022_physeq = list(
      plant_type = list(
        c("Miscanthus", "Grass"),
        c("Miscanthus", "Maize"),
        c("Maize", "Grass")
      )
    ),
    ef_physeq = list(
      genotype = list(c("MXG", "SB"), c("MXG", "ZM"), c("ZM", "SB"))
    )
  )
)
# Helper
get_comparisons <- function(data, physeq_name, comparison_map) {
  # Extract project/dataset identifier from physeq_name
  dataset_key <- gsub(".*_([^_]+)_physeq$", "\\1", physeq_name)

  # Get comparisons for this dataset if they exist
  if (!is.null(comparison_map$crop[[dataset_key]])) {
    # Find which variable exists in the data
    available_vars <- names(data)[
      names(data) %in% names(comparison_map$crop[[dataset_key]])
    ]

    if (length(available_vars) > 0) {
      return(comparison_map$crop[[dataset_key]][[available_vars[1]]])
    }
  }

  return(NULL)
}

# Alpha diversity plots for each phyloseq object
alpha_diversity_plots_16S <-
  # purrr::imap(
  # alpha_diversity_results,
  # function(project_list, project_name) {
  purrr::imap(alpha_diversity_results_16S, function(alpha_data, physeq_name) {
    # Get appropriate comparisons for this dataset
    comparisons <- get_comparisons(alpha_data, physeq_name, comparison_map)
    # Determine grouping variable
    group_var <- names(alpha_data)[
      names(alpha_data) %in% names(comparison_map)
    ][1]

    comparison_setting <- {
      if (is.null(comparisons)) {
        ggpubr::stat_compare_means(
          comparisons = comparisons,
          method = "wilcox.test",
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
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 6),
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
  })
#   }
# )

main_alpha_diversity_plots <- purrr::imap(
  main_alpha_diversity_results,
  function(project_list, project_name) {
    purrr::imap(project_list, function(alpha_data, physeq_name) {
      # Get appropriate comparisons for this dataset
      comparisons <- get_comparisons(alpha_data, physeq_name, comparison_map)
      # Determine grouping variable
      group_var <- names(alpha_data)[
        names(alpha_data) %in% names(comparison_map)
      ][1]

      comparison_setting <- {
        if (is.null(comparisons)) {
          ggpubr::stat_compare_means(
            comparisons = comparisons,
            method = "wilcox.test",
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
          plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 6),
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
    })
  }
)

cat("\n", rep("=", 40), "\n")
cat("Alpha Diversity by Project and Crop\n")
cat(rep("=", 40), "\n")
print(alpha_diversity_plots_16S$ef_16S_DNA$all_div_plots)
print(main_alpha_diversity_plots$ef_physeq_list$ef_16S_physeq$all_div_plots) #These should be too different

print(alpha_diversity_plots_16S)
print(main_alpha_diversity_plots)
# TODO
# Fix comparisons in alpha diversity - priority: normal

#--------------------------------------------------------
# SECTION 3: Beta Diversity Analysis ----
#--------------------------------------------------------

# # Calculate beta diversity (Bray-Curtis dissimilarity)
beta_diversity_results_16S <- calculate_beta_diversity(
  main_16S_physeq_list,
  method = "PCoA",
  distance = "bray"
)

main_beta_diversity_results <- calculate_alpha_diversity_nested(
  main_physeq_list,
  method = "PCoA",
  distance = "bray"
)


#--------------------------------------------------------
# SECTION 4: Beta Diversity Plots (PCoA)
#--------------------------------------------------------

# Create beta diversity plots for each phyloseq object
beta_diversity_plots <-
  # purrr::imap(
  # beta_diversity_results,
  # function(project_list, project_name) {
  purrr::imap(beta_diversity_results, function(beta_data, physeq_name) {
    plot_title <- gsub("_physeq$", " ", physeq_name) %>%
      str_to_upper(.)
    # PCoA plot colored by crop
    p_pcoa <- plot_ordination(
      beta_data$physeq,
      beta_data$ordination,
      type = "samples",
      color = "crop"
    ) +
      geom_point(size = 3, alpha = 0.7) +
      labs(
        title = paste("PCoA (Bray-Curtis) -", plot_title),
        color = "Crop"
      ) +
      theme_bw() +
      theme(legend.position = "right")

    return(p_pcoa)
  })
# }
# )

# Beta diversity plots are stored in the nested list

cat("\n", rep("=", 40), "\n")
cat("Beta Diversity Project and Crop\n")
cat(rep("=", 40), "\n")
print(beta_diversity_plots)


# Checkpoint
alpha_beta_plots <- list(
  alpha_diversity_plots = alpha_diversity_plots,
  beta_diversity_plots = beta_diversity_plots
)
# save(
#   alpha_beta_plots,
#   file = "data/output/rdata/alpha_beta_plots.rda"
# )

#--------------------------------------------------------
# SECTION 5: PERMANOVA Analysis
#--------------------------------------------------------

# # Perform PERMANOVA for each phyloseq object
# set.seed(123)
# permanova_results <- purrr::imap(
#   main_mxg_physeq_list,
#   function(project_list, project_name) {
#     purrr::imap(project_list, function(physeq_obj, physeq_name) {
#       # Get distance matrix
#       dist_matrix <- phyloseq::distance(physeq_obj, method = "bray")

#       # Get sample data
#       sampledf <- data.frame(sample_data(physeq_obj))

#       # Run PERMANOVA (adonis2)
#       perm_result <- vegan::adonis2(
#         dist_matrix ~ crop,
#         data = sampledf,
#         permutations = 999
#       )

#       # Return results with metadata
#       list(
#         physeq_name = physeq_name,
#         project = project_name,
#         permanova = perm_result,
#         n_samples = nsamples(physeq_obj)
#       )
#     })
#   }
# )

# # Display PERMANOVA results
# cat("\n", rep("=", 60), "\n")
# cat("PERMANOVA Results Summary\n")
# cat(rep("=", 60), "\n\n")

# purrr::iwalk(permanova_results, function(project_list, project_name) {
#   cat("\nPROJECT:", toupper(project_name), "\n")
#   cat(rep("-", 40), "\n")

#   purrr::iwalk(project_list, function(result, physeq_name) {
#     cat("\n", physeq_name, "(n =", result$n_samples, "samples)\n")
#     print(result$permanova)
#     cat("\n")
#   })
# })

#--------------------------------------------------------
# SECTION 6: Summary Statistics
#--------------------------------------------------------

# Summary table of alpha diversity by project and crop
alpha_summary <- alpha_diversity_df %>%
  group_by(project, region, crop) %>%
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
