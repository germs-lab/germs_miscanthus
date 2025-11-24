#' Add prevalent and rare taxa subsets as alternative experiments to a TreeSummarizedExperiment
#'
#' Iterates over specified prevalence thresholds and adds both prevalent and rare taxa
#' subsets as alternative experiments to the provided TreeSummarizedExperiment object.
#'
#' @param ps A TreeSummarizedExperiment (from mia) object.
#' @param prev_name Optional string prefix for prevalent altExp names. If NULL, uses "prevalent".
#' @param rare_name Optional string prefix for rare altExp names. If NULL, uses "rare".
#' @param thresholds Numeric vector of prevalence thresholds (as percentages, e.g. c(90, 80, 70, 60)).
#' @param rank Taxonomic rank for prevalence calculation (default "phylum").
#' @param rank_name Optional string to append to altExp names, e.g. "Phylum" (default NULL).
#' @param detection Detection threshold for mia functions (default 1/100).
#' @param assay.type Assay type for mia functions (default "counts").
#'
#' @return A named list of phyloseq objects. Each element corresponds to a prevalent or rare subset for a given threshold, with names constructed as "<prev_name>_<threshold>" or "<rare_name>_<threshold>" (e.g., "prevalent_90", "rare_80"). The list may contain both prevalent and rare subsets for each threshold, depending on the data.

get_prevalent_rare <- function(
  ps,
  prev_name = NULL,
  rare_name = NULL,
  thresholds = c(90, 80, 70, 60),
  detection = 1 / 100,
  include.lowest = FALSE
) {
  result_list <- list()

  for (thr in thresholds) {
    # Build names
    prev_suffix <- paste0("_", thr)

    prev_exp_name <- if (!is.null(prev_name)) {
      paste0(prev_name, prev_suffix)
    } else {
      paste0("prevalent", prev_suffix)
    }
    rare_exp_name <- if (!is.null(rare_name)) {
      paste0(rare_name, prev_suffix)
    } else {
      paste0("rare", prev_suffix)
    }

    # Get core/prevalent taxa members
    core_taxa <- microbiome::core_members(
      ps,
      detection = detection,
      prevalence = thr / 100,
      include.lowest = include.lowest
    )

    # Get rare taxa members
    rare_taxa <- microbiome::rare_members(
      ps,
      detection = detection,
      prevalence = thr / 100,
      include.lowest = include.lowest
    )

    # Create phyloseq subsets
    if (length(core_taxa) > 0) {
      prevalent_ps <- prune_taxa(core_taxa, ps)
      result_list[[prev_exp_name]] <- prevalent_ps
    }

    if (length(rare_taxa) > 0) {
      rare_ps <- prune_taxa(rare_taxa, ps)
      result_list[[rare_exp_name]] <- rare_ps
    }
  }

  return(result_list)
}
