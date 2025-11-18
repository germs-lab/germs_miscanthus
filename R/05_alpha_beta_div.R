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
# SECTION 1: Alpha Diversity Analysis
#--------------------------------------------------------

# Calculate alpha diversity
alpha_diversity_results <- calculate_alpha_diversity_nested(
  main_physeq_list
)

# Flatten the nested list for easier plotting
alpha_diversity_df <- alpha_diversity_results %>%
  purrr::list_flatten(name_spec = "{outer}_{inner}") %>%
  purrr::map_dfr(~.x) %>%
  select(-c(sample_id.y, target_region)) %>%
  relocate(sequence_id, .before = sample_id) %>%
  relocate(c(observed:physeq_name), .after = site)

# View summary
summary(alpha_diversity_df[, c(
  "observed",
  "shannon",
  "simpson",
  "inv_simpson"
)])

#--------------------------------------------------------
# SECTION 2: Alpha Diversity Plots
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

# Create alpha diversity plots for each phyloseq object
alpha_diversity_plots <- purrr::imap(
  alpha_diversity_results,
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
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
      )
      # Shannon diversity plot
      p_shannon <- ggplot(
        alpha_data,
        aes(x = .data[[group_var]], y = shannon, fill = .data[[group_var]])
      ) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(
          title = "Shannon Diversity",
          subtitle = physeq_name,
          x = group_var,
          y = "Shannon Index"
        ) +
        comparison_setting +
        theme_settings

      # Observed richness plot
      p_observed <- ggplot(
        alpha_data,
        aes(x = .data[[group_var]], y = observed, fill = .data[[group_var]])
      ) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(
          title = "Observed Richness",
          subtitle = physeq_name,
          x = group_var,
          y = "Observed ASVs"
        ) +
        comparison_setting +
        theme_settings

      # Simpson diversity plot
      p_simpson <- ggplot(
        alpha_data,
        aes(x = .data[[group_var]], y = simpson, fill = .data[[group_var]])
      ) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(
          title = "Simpson Diversity",
          subtitle = physeq_name,
          x = group_var,
          y = "Simpson Index"
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
# Alpha diversity plots
# my_comparisons <- list(c("MXG", "SB"), c("MXG", "ZM"), c("ZM", "SB"))

# # Determine grouping variable
# group_var <- names(alpha_diversity_results$ef_physeq_list$ef_16S_physeq)[
#   names(alpha_diversity_results$ef_physeq_list$ef_16S_physeq) %in%
#     names(comparison_map)
# ][1]

# comparison_setting <- {
#   if (is.null(comparisons)) {
#     ggpubr::stat_compare_means(
#       comparisons = comparisons,
#       method = "anova",
#       label = "p.signif"
#       # symnum.args = list(
#       #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
#       #   symbols = c("****", "***", "**", "*", "ns")
#       # )
#     )
#   }
# }
# p_shannon <- ggplot(
#   alpha_diversity_results$ef_physeq_list$ef_16S_physeq,
#   aes(x = .data[[group_var]], y = shannon, fill = .data[[group_var]])
# ) +
#   geom_boxplot(alpha = 0.7) +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   labs(
#     title = "Shannon Diversity",
#     #subtitle = physeq_name,
#     x = group_var,
#     y = "Shannon Index"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "none"
#   ) +
#   ggpubr::stat_compare_means(
#     comparisons = my_comparisons,
#     method = "anova",
#     label = "p.signif"
#     # symnum.args = list(
#     #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
#     #   symbols = c("****", "***", "**", "*", "ns")
#     # )
#   )

#   p_shannon

cat("\n", rep("=", 40), "\n")
cat("Alpha Diversity by Project and Crop\n")
cat(rep("=", 40), "\n")
print(alpha_diversity_plots)
# TODO
# Fix comparisons in alpha diversity

#--------------------------------------------------------
# SECTION 3: Beta Diversity Analysis
#--------------------------------------------------------

# Calculate beta diversity (Bray-Curtis dissimilarity)
beta_diversity_results <- calculate_beta_diversity_nested(
  main_physeq_list,
  method = "PCoA",
  distance = "bray"
)

#--------------------------------------------------------
# SECTION 4: Beta Diversity Plots (PCoA)
#--------------------------------------------------------

# Create beta diversity plots for each phyloseq object
beta_diversity_plots <- purrr::imap(
  beta_diversity_results,
  function(project_list, project_name) {
    purrr::imap(project_list, function(beta_data, physeq_name) {
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
  }
)

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
