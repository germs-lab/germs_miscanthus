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

# Calculate alpha diversity metrics for nested phyloseq objects
calculate_alpha_diversity_nested <- function(nested_list) {
  purrr::imap(nested_list, function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      # Calculate alpha diversity indices
      alpha_div <- estimate_richness(physeq_obj, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
      
      # Add sample data
      sample_df <- data.frame(sample_data(physeq_obj))
      alpha_div$sample_id <- rownames(alpha_div)
      
      # Merge with sample metadata
      alpha_div_full <- merge(alpha_div, sample_df, by.x = "sample_id", by.y = "row.names")
      
      # Add project and region information
      alpha_div_full$project <- project_name
      alpha_div_full$region <- gsub(".*_([^_]+)_physeq$", "\\1", physeq_name)
      alpha_div_full$physeq_name <- physeq_name
      
      return(alpha_div_full)
    })
  })
}

# Calculate alpha diversity
alpha_diversity_results <- calculate_alpha_diversity_nested(main_mxg_physeq_list)

# Flatten the nested list for easier plotting
alpha_diversity_df <- alpha_diversity_results %>%
  purrr::list_flatten(name_spec = "{outer}_{inner}") %>%
  purrr::map_dfr(~.x)

# View summary
summary(alpha_diversity_df[, c("Observed", "Shannon", "Simpson", "InvSimpson")])

#--------------------------------------------------------
# SECTION 2: Alpha Diversity Plots
#--------------------------------------------------------

# Create alpha diversity plots for each phyloseq object
alpha_diversity_plots <- purrr::imap(alpha_diversity_results, function(project_list, project_name) {
  purrr::imap(project_list, function(alpha_data, physeq_name) {
    
    # Shannon diversity plot
    p_shannon <- ggplot(alpha_data, aes(x = crop, y = Shannon, fill = crop)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = paste("Shannon Diversity -", physeq_name),
        x = "Crop",
        y = "Shannon Index"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # Observed richness plot
    p_observed <- ggplot(alpha_data, aes(x = crop, y = Observed, fill = crop)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = paste("Observed Richness -", physeq_name),
        x = "Crop",
        y = "Number of ASVs"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # Simpson diversity plot
    p_simpson <- ggplot(alpha_data, aes(x = crop, y = Simpson, fill = crop)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = paste("Simpson Diversity -", physeq_name),
        x = "Crop",
        y = "Simpson Index"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    # Return list of plots
    list(
      shannon = p_shannon,
      observed = p_observed,
      simpson = p_simpson
    )
  })
})

# Alpha diversity plots are stored in the nested list
# Access individual plots with: alpha_diversity_plots$project$physeq_name$plot_type
# Example: alpha_diversity_plots$mxg_ef$ef_16S_physeq$shannon

#--------------------------------------------------------
# SECTION 3: Beta Diversity Analysis
#--------------------------------------------------------

# Calculate beta diversity (ordination) for nested phyloseq objects
calculate_beta_diversity_nested <- function(nested_list, method = "bray") {
  purrr::imap(nested_list, function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      # Calculate ordination
      ord <- ordinate(physeq_obj, method = "PCoA", distance = method)
      
      # Return both ordination and phyloseq object for plotting
      list(
        ordination = ord,
        physeq = physeq_obj,
        physeq_name = physeq_name
      )
    })
  })
}

# Calculate beta diversity (Bray-Curtis dissimilarity)
beta_diversity_results <- calculate_beta_diversity_nested(main_mxg_physeq_list, method = "bray")

#--------------------------------------------------------
# SECTION 4: Beta Diversity Plots (PCoA)
#--------------------------------------------------------

# Create beta diversity plots for each phyloseq object
beta_diversity_plots <- purrr::imap(beta_diversity_results, function(project_list, project_name) {
  purrr::imap(project_list, function(beta_data, physeq_name) {
    
    # PCoA plot colored by crop
    p_pcoa <- plot_ordination(
      beta_data$physeq,
      beta_data$ordination,
      type = "samples",
      color = "crop"
    ) +
      geom_point(size = 3, alpha = 0.7) +
      labs(
        title = paste("PCoA (Bray-Curtis) -", beta_data$physeq_name),
        color = "Crop"
      ) +
      theme_bw() +
      theme(legend.position = "right")
    
    return(p_pcoa)
  })
})

# Beta diversity plots are stored in the nested list
# Access individual plots with: beta_diversity_plots$project$physeq_name
# Example: beta_diversity_plots$mxg_ef$ef_16S_physeq

#--------------------------------------------------------
# SECTION 5: PERMANOVA Analysis
#--------------------------------------------------------

# Perform PERMANOVA for each phyloseq object
set.seed(123)  # Set seed once for reproducibility
permanova_results <- purrr::imap(main_mxg_physeq_list, function(project_list, project_name) {
  purrr::imap(project_list, function(physeq_obj, physeq_name) {
    
    # Get distance matrix
    dist_matrix <- phyloseq::distance(physeq_obj, method = "bray")
    
    # Get sample data
    sampledf <- data.frame(sample_data(physeq_obj))
    
    # Run PERMANOVA (adonis2)
    perm_result <- vegan::adonis2(dist_matrix ~ crop, data = sampledf, permutations = 999)
    
    # Return results with metadata
    list(
      physeq_name = physeq_name,
      project = project_name,
      permanova = perm_result,
      n_samples = nsamples(physeq_obj)
    )
  })
})

# Display PERMANOVA results
cat("\n", rep("=", 60), "\n")
cat("PERMANOVA Results Summary\n")
cat(rep("=", 60), "\n\n")

purrr::iwalk(permanova_results, function(project_list, project_name) {
  cat("\nPROJECT:", toupper(project_name), "\n")
  cat(rep("-", 40), "\n")
  
  purrr::iwalk(project_list, function(result, physeq_name) {
    cat("\n", physeq_name, "(n =", result$n_samples, "samples)\n")
    print(result$permanova)
    cat("\n")
  })
})

#--------------------------------------------------------
# SECTION 6: Summary Statistics
#--------------------------------------------------------

# Summary table of alpha diversity by project and crop
alpha_summary <- alpha_diversity_df %>%
  group_by(project, region, crop) %>%
  summarise(
    n = n(),
    mean_observed = mean(Observed, na.rm = TRUE),
    sd_observed = sd(Observed, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    mean_simpson = mean(Simpson, na.rm = TRUE),
    sd_simpson = sd(Simpson, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n", rep("=", 60), "\n")
cat("Alpha Diversity Summary by Project and Crop\n")
cat(rep("=", 60), "\n")
print(alpha_summary)
