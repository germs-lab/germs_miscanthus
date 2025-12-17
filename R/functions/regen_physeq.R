regen_physeq <- function(physeq, sample_metadata, rownames = "sample_id") {
  if (!is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))) {
    physeq_regen <- phyloseq(
      otu_table(physeq, taxa_are_rows = TRUE),
      tax_table(physeq),
      sample_data(
        sample_metadata %>% column_to_rownames(., var = rownames)
      ),
      phyloseq::refseq(physeq)
    )
  }

  if (is.null(phyloseq::refseq(physeq, errorIfNULL = FALSE))) {
    physeq_regen <- phyloseq(
      otu_table(physeq, taxa_are_rows = TRUE),
      tax_table(physeq),
      sample_data(
        sample_metadata %>% column_to_rownames(., var = rownames)
      )
    )

    physeq_regen <- add_refseq(physeq_regen, tag = NA, seq_to_name = TRUE)
  }

  return(physeq_regen)
}
