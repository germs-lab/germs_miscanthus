create_binary_df_from_flat <- function(flat_list) {
  all_items <- unique(unlist(flat_list))
  set_names <- names(flat_list)
  binary_df <- map_dfc(set_names, function(set_name) {
    tibble(!!set_name := as.integer(all_items %in% flat_list[[set_name]]))
  })
  binary_df$item <- all_items
  binary_df |> relocate(item)
}

list_to_binary_df <- function(presence_list) {
  flattened <- purrr::flatten(presence_list)
  create_binary_df_from_flat(flattened)
}

get_mxg_sets <- function(presence_list) {
  flattened <- purrr::flatten(presence_list)
  mxg_patterns <- c("MXG", "Miscanthus", "_M$")
  mxg_idx <- grep(paste(mxg_patterns, collapse = "|"), names(flattened))
  flattened[mxg_idx]
}
