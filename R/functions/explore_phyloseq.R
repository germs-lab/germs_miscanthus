explore_nested_phyloseq <- function(nested_list) {
  results <- purrr::imap(nested_list, function(project_list, project_name) {
    cat("\n", rep("=", 40), "\n")
    cat("PROJECT:", toupper(project_name), "\n")
    cat(rep("=", 40), "\n")

    purrr::imap(project_list, function(physeq_obj, seq_type) {
      cat("\n", rep("-", 40), "\n")
      cat("Sequencing Type:", toupper(seq_type), "\n")
      cat(rep("-", 40), "\n")

      cat("Basic Summary:\n")
      print(metagMisc::phyloseq_summary(
        physeq_obj,
        more_stats = FALSE,
        long = FALSE
      ))

      cat("\nRead/Sequencing Summary:\n")
      print(microbiome::summarize_phyloseq(physeq_obj))

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
