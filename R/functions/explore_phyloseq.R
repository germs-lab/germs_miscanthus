explore_phyloseq_list <- function(physeq_obj, obj_name) {
  cat("=== Analysis for:", obj_name, "===\n")

  # Basic dimensions
  cat("OTU table dimensions:", dim(otu_table(physeq_obj)), "\n")
  cat("Number of taxa:", ntaxa(physeq_obj), "\n")
  cat("Number of samples:", nsamples(physeq_obj), "\n")

  # Summary statistics
  cat("\nRaw count summary:\n")
  print(summary(as.vector(otu_table(physeq_obj))))

  # # Sample sums
  # cat("\nSample sums (total reads per sample):\n")
  # print(sample_sums(physeq_obj))

  # cat("\n", rep("=", 50), "\n\n")
}
