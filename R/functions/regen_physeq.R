regen_physeq <- function(physeq, sample_metadata, rownames = "sample_id") {
  phyloseq(
    otu_table(physeq),
    tax_table(physeq),
    sample_data(sample_metadata %>% column_to_rownames(., var = rownames))
  )
}
