calculate_beta_diversity <- function(
  phyloseq_list,
  method = "NMDS",
  distance = "bray",
  nested = FALSE,
  join_df = FALSE
) {
  if (nested) {
    if (join_df) {
      # Join all phyloseq objects into a single data frame for analysis
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
      return(combined_df)
    } else {
      purrr::imap(phyloseq_list, function(project_list, project_name) {
        purrr::imap(project_list, function(physeq_obj, physeq_name) {
          result <- .beta_workflow(
            physeq_obj = physeq_obj,
            physeq_name = physeq_name,
            method = method,
            distance = distance
          )
        })
      })
    }
  } else {
    purrr::imap(phyloseq_list, function(physeq_obj, physeq_name) {
      result <- .beta_workflow(
        physeq_obj = physeq_obj,
        physeq_name = physeq_name,
        method = method,
        distance = distance
      )
    })
  }
}

.beta_workflow <- function(physeq_obj, physeq_name, method, distance) {
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


beta_plot_workflow <- function(
  beta_data,
  physeq_name,
  type = "samples",
  color = "crop"
) {
  plot_title <- gsub("_physeq$", " ", physeq_name) %>%
    str_to_upper(.)
  # Ordination plot colored by crop
  p_ordination <- plot_ordination(
    beta_data$physeq,
    beta_data$ordination,
    type = type,
    color = color
  ) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = ifelse(
        class(beta_data$ordination) == "nmds",
        paste("NMDS (Bray-Curtis) -", plot_title),
        paste("PCoA (Bray-Curtis) -", plot_title)
      ),
      color = str_to_title(color)
    ) +
    theme_bw() +
    theme(legend.position = "right")

  return(p_ordination)
}
