#' Export Reference Sequences from Phyloseq Objects to FASTA Files
#'
#' Extracts reference sequences from phyloseq objects and exports them to FASTA
#' format with explicit, standardized naming.
#'
#' @param physeq_obj A a phyloseq object or named list of phyloseq objects.
#' @param out_dir Output directory for FASTA files.
#' @param crop_subset Character string describing the crop subset
#'   (e.g., "all_crops", "mxg_only"). Default is NULL.
#' @param out_ext File extension. Default is ".fa".
#'
#' @return Invisibly returns a named list of output file paths.
#'
#' @details
#' Output naming convention: {project}_{target_region}_{nucleotide}_{crop_subset}.fa
#'
#' Examples:
#' - ef_16S_DNA_all_crops.fa (Energy Farm, 16S, DNA, all crops)
#' - lamps_2018_ITS_DNA_mxg_only.fa (LAMPS 2018, ITS, DNA, Miscanthus only)
#' - lamps_2022_AMF_DNA_all_crops.fa (LAMPS 2022, AMF, DNA, all crops)
#'
#' @export
refseq2fasta <- function(
  physeq_obj,
  out_dir,
  crop_subset = NULL,
  out_ext = ".fa"
) {
  # Validate inputs
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  .phyloseq_export_fasta <- function(ps, nm) {
    seqs <- phyloseq::refseq(ps)

    # Build explicit filename with crop subset
    fname <- if (!is.null(crop_subset)) {
      paste0(nm, "_", crop_subset, out_ext)
    } else {
      paste0(nm, out_ext)
    }

    # Sanitize filename (remove spaces, keep alphanumeric, underscore, hyphen, dot)
    fname <- gsub("[[:space:]]+", "_", fname)
    fname <- gsub("[^A-Za-z0-9_\\-\\.]", "", fname)

    outfile <- file.path(out_dir, fname)

    # Convert Biostrings to phylotools format
    new_seqs <- Biostrings::as.data.frame(seqs) %>%
      tibble::rownames_to_column(., var = "seq.name") %>%
      dplyr::rename(seq.text = x)

    tryCatch(
      {
        phylotools::dat2fasta(new_seqs, outfile = outfile)
        message("Exported: ", fname, " (", nrow(new_seqs), " sequences)")
      },
      error = function(e) {
        warning("Failed to write ", outfile, ": ", conditionMessage(e))
      }
    )

    outfile
  }

  # Check if list or single object
  if (is.list(physeq_obj) && all(sapply(physeq_obj, inherits, "phyloseq"))) {
    # Get base names from list (remove _physeq suffix)
    base_names <- names(physeq_obj)
    base_names <- gsub("_physeq$", "", base_names)

    # Process each phyloseq object and collect output paths
    output_files <- purrr::map2(
      physeq_obj,
      base_names,
      .phyloseq_export_fasta
    )

    names(output_files) <- base_names
  } else {
    # Single phyloseq: extract name after $ if it exists
    obj_expr <- deparse(substitute(physeq_obj))
    obj_name <- if (grepl("\\$", obj_expr)) {
      sub(".*\\$", "", obj_expr) # Extract after $
    } else {
      obj_expr
    }

    output_files <- .phyloseq_export_fasta(physeq_obj, obj_name)
    names(output_files) <- obj_name
  }

  invisible(output_files)
}
