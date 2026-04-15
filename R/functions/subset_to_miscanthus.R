#' Subset phyloseq objects to Miscanthus crops with optional filtering
#'
#' @param physeq_input Either a single phyloseq object or a nested list of phyloseq objects
#' @param crop_patterns Character vector of crop patterns to match (default: c("MXG", "M", "Miscanthus"))
#' @param filter_16S_only Logical. If TRUE, only keep 16S objects in nested lists (default: FALSE)
#' @param filter_DNA_only Logical. If TRUE, subset to DNA nucleotide only (default: FALSE)
#' @param rename_suffix Character. Suffix to add when renaming (default: "_DNA")
#'
#' @return Filtered phyloseq object or nested list structure
subset_to_miscanthus <- function(
  physeq_input,
  crop_patterns = c("MXG", "M", "Miscanthus"),
  filter_16S_only = FALSE,
  filter_DNA_only = FALSE,
  rename_suffix = "_DNA"
) {
  # Helper function to subset a single phyloseq object
  .subset_single_physeq <- function(
    ps,
    check_nucleotide = TRUE,
    .crop_patterns = crop_patterns
  ) {
    keep <- sample_data(ps)$crop %in% .crop_patterns
    ps_subset <- prune_samples(keep, ps) |>
      filter_taxa(function(x) sum(x) > 0, TRUE)

    if (
      check_nucleotide &&
        filter_DNA_only &&
        "nucleotide" %in% colnames(sample_data(ps_subset))
    ) {
      keep_dna <- sample_data(ps_subset)$nucleotide == "DNA"
      ps_subset <- prune_samples(keep_dna, ps_subset) |>
        filter_taxa(function(x) sum(x > 0) > 0, TRUE)
    }

    return(ps_subset)
  }

  # Check if input is a single phyloseq object
  if (inherits(physeq_input, "phyloseq")) {
    return(.subset_single_physeq(physeq_input, check_nucleotide = TRUE))
  }

  # Handle nested list structure
  if (is.list(physeq_input)) {
    result <- purrr::imap(
      physeq_input,
      function(project_list, project_name) {
        if (inherits(project_list, "phyloseq")) {
          return(.subset_single_physeq(
            project_list,
            check_nucleotide = TRUE,
            .crop_patterns = crop_patterns
          ))
        }
        purrr::imap(project_list, function(physeq_obj, physeq_name) {
          # Skip non-16S objects if filtering
          if (
            filter_16S_only && !grepl("16S", physeq_name, ignore.case = TRUE)
          ) {
            return(NULL)
          }

          # Subset the phyloseq object
          .subset_single_physeq(
            physeq_obj,
            check_nucleotide = TRUE,
            .crop_patterns = crop_patterns
          )
        }) %>%
          # Keep only non-NULL objects and rename
          purrr::discard(is.null) %>%
          set_names(names(.) %>% str_replace("_physeq$", rename_suffix))
      }
    ) %>%
      set_names(names(.) %>% str_remove("_physeq_list$"))

    return(result)
  }

  stop("Input must be a phyloseq object or a nested list of phyloseq objects")
}
