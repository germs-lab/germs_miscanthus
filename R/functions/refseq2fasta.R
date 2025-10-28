refseq2fasta <- function(
  phyloseq_list,
  out_dir,
  extra_id = NULL,
  out_ext = ".fa"
) {
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
      fname <- paste0(.y, extra_id, out_ext)
      fname <- gsub("[[:space:]]+", "_", fname)
      fname <- gsub("[^A-Za-z0-9_\\-\\.]", "", fname)

      outfile <- file.path(out_dir, fname)

      # Getting Biostring to talk to phylotools
      new_seqs <- Biostrings::as.data.frame(seqs) %>%
        tibble::rownames_to_column(., var = "seq.name") %>%
        rename(seq.text = x)

      tryCatch(
        phylotools::dat2fasta(new_seqs, outfile = outfile),
        error = function(e) {
          message("Failed to write ", outfile, ": ", conditionMessage(e))
        }
      )
    }
  )
}
