refseq2fasta <- function(phyloseq_list, out_dir, out_ext = ".fa") {
  # Names for files
  nms <- names(phyloseq_list)
  nms <- gsub("_physeq", "", nms)

  # Refseqs to a FASTA named by the list name
  purrr::walk2(
    phyloseq_list,
    nms,
    ~ {
      seqs <- phyloseq::refseq(.x)

      # build a safe output path
      fname <- paste0(.y, out_ext)
      fname <- gsub("[[:space:]]+", "_", fname)
      fname <- gsub("[^A-Za-z0-9_\\-\\.]", "", fname)

      outfile <- file.path(out_dir, fname)

      tryCatch(
        Biostrings::writeXStringSet(
          seqs,
          filepath = outfile,
          format = "fasta",
          append = FALSE,
          compress = FALSE
        ),
        error = function(e) {
          message("Failed to write ", outfile, ": ", conditionMessage(e))
        }
      )
    }
  )
}
