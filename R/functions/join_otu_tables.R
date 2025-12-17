join_otu_tables <- function(table_a, table_b, table_c) {
  combined_seqtab_DT <- suppressMessages({
    full_join(
      {
        full_join(
          table_a,
          table_b,
          by = NULL
        )
      },
      table_c,
      by = NULL
    )
  }) %>%
    data.table::as.data.table(.) # We need some speed with this big data

  # Replace resulting NAs
  for (j in seq_len(ncol(combined_seqtab_DT))) {
    data.table::set(
      combined_seqtab_DT,
      which(is.na(combined_seqtab_DT[[j]])),
      j,
      0
    )
  }

  # Ensuring types
  combined_seqtab_DT[,
    (names(combined_seqtab_DT)) := lapply(.SD, function(x) {
      if (is.numeric(x)) as.integer(x) else x
    })
  ]

  new_seqtab <- combined_seqtab_DT %>%
    column_to_rownames(., var = "sample_id") %>%
    as.matrix()

  return(new_seqtab)
}
