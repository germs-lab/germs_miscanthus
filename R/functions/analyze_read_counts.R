# Analyze read counts across nested phyloseq objects

analyze_read_counts <- function(nested_list) {
  read_summaries <- purrr::imap(
    nested_list,
    function(project_list, project_name) {
      purrr::imap(project_list, function(physeq_obj, seq_type) {
        # Get read counts
        reads <- readcount(physeq_obj) %>%
          as.data.frame() %>%
          rownames_to_column(var = "sample_id") %>%
          rename(n_seqs = ".") %>%
          group_by(sample_id) %>%
          mutate(
            n_singletons = sum(n_seqs == 1),
            goods = 1 - (n_singletons / n_seqs),
            project = gsub("_physeq", "", seq_type),
            region = gsub(".*_([^_]+)_physeq$", "\\1", seq_type)
          )

        metadata <- physeq_obj %>%
          physeq2df() %>%
          select(sample_id, crop)

        new_df <- dplyr::left_join(reads, metadata, by = "sample_id")

        return(new_df)
      })
    }
  ) %>%
    purrr::list_flatten(name_spec = "{outer}_{inner}") %>%
    purrr::map_dfr(~.x)

  return(read_summaries)
}
