# Examples and Tests for p_iNEXT_improved.R
#
# This script demonstrates usage of the improved parallel iNEXT implementation
# and compares it with the original version.

library(iNEXT)
library(ggplot2)

# Source the improved version
source("R/functions/p_iNEXT_improved.R")

# ============================================================================
# Example 1: Using Spider Dataset (Built-in iNEXT Data)
# ============================================================================

# Load example data from iNEXT package
data(spider)

# Spider data is a list with multiple sites
# Convert to matrix format (samples as rows, species as columns)
spider_matrix <- do.call(rbind, spider)
rownames(spider_matrix) <- names(spider)

cat("\n=== Example 1: Spider Dataset ===\n")
cat(
  "Matrix dimensions:",
  nrow(spider_matrix),
  "samples x",
  ncol(spider_matrix),
  "species\n"
)

# Run improved p_iNEXT with combined output (compatible with ggiNEXT)
result_combined <- p_iNEXT(
  x = spider_matrix,
  q = c(0, 1, 2),
  nCores = 2,
  combine = TRUE, # Returns single iNEXT object
  verbose = TRUE
)

# Verify output structure
cat("\nOutput class:", class(result_combined), "\n")
cat("Output components:", names(result_combined), "\n")

# Plot with standard ggiNEXT
p1 <- ggiNEXT(
  result_combined,
  type = 1,
  facet.var = "Order.q",
  color.var = "Assemblage"
)
print(p1 + ggtitle("Spider Dataset - Rarefaction Curves (Improved p_iNEXT)"))

# ============================================================================
# Example 2: Compare with Original iNEXT
# ============================================================================

cat("\n=== Example 2: Comparison with Original iNEXT ===\n")

# Original iNEXT (for comparison)
result_original <- iNEXT(spider, q = c(0, 1, 2), datatype = "abundance")

# Improved p_iNEXT
result_improved <- p_iNEXT(
  x = spider_matrix,
  q = c(0, 1, 2),
  nCores = 2,
  combine = TRUE
)

# Both should produce similar plots
p_orig <- ggiNEXT(result_original, type = 1, facet.var = "Order.q")
p_impr <- ggiNEXT(result_improved, type = 1, facet.var = "Order.q")

cat(
  "Original iNEXT output structure matches improved version:",
  all(names(result_original) == names(result_improved)),
  "\n"
)

# ============================================================================
# Example 3: Using List Output (Individual iNEXT Objects)
# ============================================================================

cat("\n=== Example 3: List Output Format ===\n")

# Get individual iNEXT objects for each sample
result_list <- p_iNEXT(
  x = spider_matrix,
  q = c(0),
  nCores = 2,
  combine = FALSE, # Returns list of iNEXT objects
  verbose = TRUE
)

cat("\nReturned list with", length(result_list), "elements\n")
cat("Sample names:", names(result_list), "\n")

# Access individual sample
sample1 <- result_list[[1]]
cat("Sample 1 class:", class(sample1), "\n")

# Plot individual sample
p_sample1 <- ggiNEXT(sample1, type = 1)
print(p_sample1 + ggtitle(paste("Individual Sample:", names(result_list)[1])))

# ============================================================================
# Example 4: Custom Plotting with Flattened Data
# ============================================================================

cat("\n=== Example 4: Custom Plotting ===\n")

# Flatten results to data frame for custom plotting
result <- p_iNEXT(spider_matrix, q = c(0, 1, 2), nCores = 2)
df <- flatten_iNEXT_results(result)

cat("Flattened data frame dimensions:", nrow(df), "x", ncol(df), "\n")
cat("Column names:", paste(head(colnames(df), 10), collapse = ", "), "...\n")

# Custom ggplot
p_custom <- ggplot(
  df,
  aes(
    x = size_based.m,
    y = size_based.qD,
    color = sample,
    group = sample
  )
) +
  geom_line(size = 1.2) +
  facet_wrap(
    ~size_based.Order.q,
    labeller = labeller(
      size_based.Order.q = c("0" = "q = 0", "1" = "q = 1", "2" = "q = 2")
    )
  ) +
  labs(
    x = "Number of Individuals",
    y = "Species Diversity",
    title = "Custom Plot - Rarefaction Curves",
    color = "Site"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_custom)

# Or use the provided custom plotting function
p_custom2 <- gg_inext_custom(
  df,
  type = 1,
  facet.var = "Order.q",
  color.var = "Assemblage",
  se = TRUE
)
print(p_custom2 + ggtitle("Using gg_inext_custom()"))

# ============================================================================
# Example 5: Error Handling
# ============================================================================

cat("\n=== Example 5: Error Handling ===\n")

# Create matrix with an empty sample (all zeros)
spider_with_empty <- rbind(
  spider_matrix,
  EmptySample = rep(0, ncol(spider_matrix))
)

# The improved version handles this gracefully
result_with_handling <- p_iNEXT(
  x = spider_with_empty,
  q = c(0),
  nCores = 2,
  verbose = TRUE
)

cat("\nSuccessfully handled empty samples\n")

# ============================================================================
# Example 6: Different Parallelization Strategies
# ============================================================================

cat("\n=== Example 6: Parallelization Strategies ===\n")

# Sequential (no parallelization)
system.time({
  result_seq <- p_iNEXT(
    spider_matrix,
    q = c(0, 1, 2),
    nCores = 1,
    plan_strategy = "sequential",
    verbose = FALSE
  )
})

# Multisession (parallel with separate R sessions)
system.time({
  result_par <- p_iNEXT(
    spider_matrix,
    q = c(0, 1, 2),
    nCores = 2,
    plan_strategy = "multisession",
    verbose = FALSE
  )
})

cat(
  "\nBoth strategies produce same results:",
  identical(result_seq$DataInfo, result_par$DataInfo),
  "\n"
)

# ============================================================================
# Example 7: Combining Results Later
# ============================================================================

cat("\n=== Example 7: Combining Results Later ===\n")

# Run with combine=FALSE
results_separate <- p_iNEXT(
  spider_matrix,
  q = c(0, 1, 2),
  nCores = 2,
  combine = FALSE,
  verbose = FALSE
)

# Combine later using helper function
results_combined_later <- combine_iNEXT_list(results_separate)

cat("Combined later - class:", class(results_combined_later), "\n")
cat("Can use with ggiNEXT:", "iNEXT" %in% class(results_combined_later), "\n")

# ============================================================================
# Example 8: Working with Phyloseq Objects (Conceptual)
# ============================================================================

cat("\n=== Example 8: Phyloseq Workflow (Conceptual) ===\n")

# This example shows the typical workflow with phyloseq objects
# Uncomment and modify for actual use

# library(phyloseq)
#
# # Load your phyloseq object
# physeq <- readRDS("path/to/phyloseq.rds")
#
# # Extract OTU table (transpose so samples are rows)
# otu_mat <- as.matrix(t(otu_table(physeq)))
#
# # Calculate endpoint
# max_lib_size <- max(rowSums(otu_mat))
# endpoint <- max_lib_size * 2
#
# # Run parallel iNEXT
# result <- p_iNEXT(
#   x = otu_mat,
#   q = c(0, 1, 2),
#   endpoint = endpoint,
#   nboot = 100,
#   nCores = 4,
#   combine = TRUE,
#   verbose = TRUE
# )
#
# # Plot
# ggiNEXT(result, type = 1, facet.var = "Order.q") +
#   theme_bw() +
#   labs(title = "Rarefaction Curves for My Dataset")

cat("\nSee above commented code for phyloseq workflow\n")

# ============================================================================
# Performance Comparison (Optional)
# ============================================================================

cat("\n=== Performance Test (Optional) ===\n")
cat("To run performance test, uncomment the code below\n")

# # Create larger synthetic dataset for performance testing
# set.seed(123)
# n_samples <- 20
# n_species <- 100
#
# # Generate synthetic abundance data
# synthetic_data <- matrix(
#   rpois(n_samples * n_species, lambda = 5),
#   nrow = n_samples,
#   ncol = n_species
# )
# rownames(synthetic_data) <- paste0("Sample", 1:n_samples)
# colnames(synthetic_data) <- paste0("Species", 1:n_species)
#
# cat("\nPerformance test with", n_samples, "samples and", n_species, "species\n")
#
# # Sequential
# cat("Sequential processing...\n")
# time_seq <- system.time({
#   result_seq <- p_iNEXT(
#     synthetic_data,
#     q = c(0, 1, 2),
#     nCores = 1,
#     nboot = 50,
#     verbose = FALSE
#   )
# })
#
# # Parallel (2 cores)
# cat("Parallel processing (2 cores)...\n")
# time_par2 <- system.time({
#   result_par2 <- p_iNEXT(
#     synthetic_data,
#     q = c(0, 1, 2),
#     nCores = 2,
#     nboot = 50,
#     verbose = FALSE
#   )
# })
#
# # Parallel (4 cores)
# cat("Parallel processing (4 cores)...\n")
# time_par4 <- system.time({
#   result_par4 <- p_iNEXT(
#     synthetic_data,
#     q = c(0, 1, 2),
#     nCores = 4,
#     nboot = 50,
#     verbose = FALSE
#   )
# })
#
# cat("\nTiming Results:\n")
# cat("Sequential:", time_seq["elapsed"], "seconds\n")
# cat("Parallel (2 cores):", time_par2["elapsed"], "seconds - Speedup:",
#     round(time_seq["elapsed"] / time_par2["elapsed"], 2), "x\n")
# cat("Parallel (4 cores):", time_par4["elapsed"], "seconds - Speedup:",
#     round(time_seq["elapsed"] / time_par4["elapsed"], 2), "x\n")

cat("\n=== All Examples Completed Successfully ===\n")
