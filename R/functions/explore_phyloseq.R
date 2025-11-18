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


explore_nested_phyloseq <- function(nested_list) {
  results <- purrr::imap(nested_list, function(project_list, project_name) {
    cat("\n", rep("=", 40), "\n")
    cat("PROJECT:", toupper(project_name), "\n")
    cat(rep("=", 40), "\n")

    purrr::imap(project_list, function(physeq_obj, seq_type) {
      cat("\n", rep("-", 40), "\n")
      cat("Sequencing Type:", toupper(seq_type), "\n")
      cat(rep("-", 40), "\n")

      # Basic phyloseq summaries
      cat("Basic Summary:\n")
      print(metagMisc::phyloseq_summary(
        physeq_obj,
        more_stats = FALSE,
        long = FALSE
      ))

      cat("\nRead/Sequencing Summary:\n")
      print(microbiome::summarize_phyloseq(physeq_obj))

      # # Taxonomic distribution
      # if ("phylum" %in% colnames(tax_table(physeq_obj))) {
      #   cat("\nPhylum Distribution:\n")
      #   phyla_dist <- phyloseq_ntaxa_by_tax(
      #     physeq_obj,
      #     TaxRank = "phylum",
      #     relative = FALSE,
      #     add_meta_data = FALSE
      #   ) %>%
      #     as.data.frame() %>%
      #     mutate(sum = sum(N.OTU)) %>%
      #     group_by(phylum) %>%
      #     summarise(occurance_in_samples = n()) %>%
      #     arrange(desc(occurance_in_samples))

      #   print(phyla_dist)
      # }

      # Sample and taxa counts
      list(
        project = project_name,
        region = gsub(".*_([^_]+)_physeq$", "\\1", seq_type),
        n_samples = nsamples(physeq_obj),
        n_taxa = ntaxa(physeq_obj),
        sample_vars = sample_variables(physeq_obj),
        rank_names = rank_names(physeq_obj),
        total_reads = sum(sample_sums(physeq_obj)),
        min_reads_per_sample = min(sample_sums(physeq_obj)),
        max_reads_per_sample = max(sample_sums(physeq_obj)),
        mean_reads_per_sample = mean(sample_sums(physeq_obj)),
        median_reads_per_sample = median(sample_sums(physeq_obj))
      )
    })
  })

  return(results)
}
