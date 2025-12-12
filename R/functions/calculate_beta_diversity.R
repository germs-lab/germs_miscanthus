beta_workflow <- function(physeq_obj, physeq_name, method, distance) {
  physeq_filtered <- physeq_obj %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)

  ord <- ordinate(physeq_filtered, method = method, distance = distance)

  # Return both ordination and phyloseq object for plotting
  return(list(
    ordination = ord,
    physeq = physeq_filtered,
    physeq_name = physeq_name
  ))
}

calculate_beta_diversity <- function(
  project_list,
  method = "NMDS",
  distance = "bray"
) {
  purrr::imap(project_list, function(physeq_obj, physeq_name) {
    # Calculate ordination
    result <- beta_workflow(
      physeq_obj = physeq_obj,
      physeq_name = physeq_name,
      method = method,
      distance = distance
    )
  })
}

# Calculate beta diversity (ordination) for nested phyloseq objects
calculate_beta_diversity_nested <- function(
  nested_list,
  method = "NMDS",
  distance = "bray"
) {
  purrr::imap(nested_list, function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      # Calculate ordination
      result <- beta_workflow(
        physeq_obj = physeq_obj,
        physeq_name = physeq_name,
        method = method,
        distance = distance
      )
    })
  })
}
