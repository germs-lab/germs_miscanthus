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
  subset_single_physeq <- function(ps, check_nucleotide = TRUE) {
    # Subset to crops
    ps_subset <- subset_samples(ps, crop %in% crop_patterns) %>%
      filter_taxa(function(x) sum(x) > 0, TRUE)

    # Optionally subset to DNA nucleotide
    if (
      check_nucleotide &&
        filter_DNA_only &&
        "nucleotide" %in% colnames(sample_data(ps_subset))
    ) {
      ps_subset <- ps_subset %>%
        subset_samples(nucleotide == "DNA") %>%
        filter_taxa(function(x) sum(x > 0) > 0, TRUE)
    }

    return(ps_subset)
  }

  # Check if input is a single phyloseq object
  if (inherits(physeq_input, "phyloseq")) {
    return(subset_single_physeq(physeq_input, check_nucleotide = TRUE))
  }

  # Handle nested list structure
  if (is.list(physeq_input)) {
    result <- purrr::imap(
      physeq_input,
      function(project_list, project_name) {
        purrr::imap(project_list, function(physeq_obj, physeq_name) {
          # Skip non-16S objects if filtering
          if (
            filter_16S_only && !grepl("16S", physeq_name, ignore.case = TRUE)
          ) {
            return(NULL)
          }

          # Subset the phyloseq object
          subset_single_physeq(physeq_obj, check_nucleotide = TRUE)
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
