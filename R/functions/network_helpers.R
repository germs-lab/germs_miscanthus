#' Prevalence-filter a phyloseq object
#' @param physeq      phyloseq object (taxa_are_rows = TRUE expected)
#' @param prev_thresh minimum fraction of samples an ASV must appear in
#' @param min_count   minimum raw count to consider a sample "present"
prev_filter <- function(physeq, prev_thresh = 0.10, min_count = 2) {
  mat <- as.matrix(otu_table(physeq))
  if (!taxa_are_rows(physeq)) {
    mat <- t(mat)
  }
  prev <- rowMeans(mat >= min_count)
  prune_taxa(prev >= prev_thresh, physeq)
}

#' Retain only shared samples between two phyloseq objects
align_taxa <- function(ps1, ps2) {
  shared <- intersect(taxa_names(ps1), taxa_names(ps2))
  if (length(shared) == 0) {
    stop("No shared taxa between the two physeq objects.")
  }
  list(
    ps1 = prune_taxa(shared, ps1),
    ps2 = prune_taxa(shared, ps2)
  )
}
align_samples <- function(ps1, ps2) {
  shared <- intersect(sample_names(ps1), sample_names(ps2))
  if (length(shared) == 0) {
    stop("No shared samples between the two physeq objects.")
  }
  list(
    ps1 = prune_samples(shared, ps1),
    ps2 = prune_samples(shared, ps2)
  )
}


#' Network-level summary statistics for one igraph object.
network_summary <- function(graph, site, kind) {
  n_v <- vcount(graph) #nodes
  n_e <- ecount(graph)
  n_pos <- sum(E(graph)$link_sign == "positive", na.rm = TRUE)
  n_neg <- sum(E(graph)$link_sign == "negative", na.rm = TRUE)

  # Edge count logic for cross-kingdom networks: count only edges between bacteria and fungi
  bf_edges <- if ("Fungi" %in% V(graph)$kingdom) {
    e_ends <- igraph::as_edgelist(graph)
    bact_v <- V(graph)$name[V(graph)$kingdom == "Bacteria"]
    fung_v <- V(graph)$name[V(graph)$kingdom == "Fungi"]
    sum(
      (e_ends[, 1] %in% bact_v & e_ends[, 2] %in% fung_v) |
        (e_ends[, 2] %in% bact_v & e_ends[, 1] %in% fung_v)
    )
  } else {
    NA_integer_
  }

  # Density and connectance calculations
  max_edges <- n_v * (n_v - 1) / 2
  multi_level_comm <- cluster_louvain(graph)
  grdy <- cluster_fast_greedy(graph)

  # Nestedness calculation
  full_nodf <- NA_real_
  full_nodf <- tryCatch(
    {
      adj <- as_adjacency_matrix(graph, sparse = FALSE)
      vegan::nestednodf(adj)$statistic[["NODF"]]
    },
    error = function(e) NA_real_
  )

  sub_matx_nodf <- NA_real_
  if (kind == "cross" && length(b) > 1 && length(f) > 1) {
    adj_b <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
    bpt <- adj_b[b, f]
    nodf <- tryCatch(
      bipartite::networklevel(bpt, index = "NODF")[["NODF"]],
      error = function(e) NA_real_
    )
  }

  summary_results <- tibble::tibble(
    site = site,
    network_type = kind,
    n_taxa = n_v,
    n_edges = n_e,
    max_edges = max_edges,
    positive_links = n_pos,
    negative_links = n_neg,
    density_connectance = edge_density(graph),
    connectance = if (max_edges > 0) n_e / max_edges else NA_real_,
    bf_edge_fraction = if (!is.na(bf_edges)) bf_edges / n_e else NA_real_,
    pos_neg_link_ratio = if (n_neg > 0) n_pos / (n_pos + n_neg) else NA_real_,
    transitivity_avgCC = transitivity(
      graph,
      type = "average",
      isolates = "zero"
    ),
    cl_modularity = modularity(multi_level_comm),
    grdy_modularity = modularity(grdy),
    full_nestedness_NODF = full_nodf,
    subnetwork_nestedness_NODF = sub_matx_nodf
  )
  return(summary_results)
}


#' Per-ASV centrality shift: bacteria-only graph vs. bacterial subgraph of the
#' cross-kingdom network.
centrality_shift <- function(g_bact, g_cross, site) {
  bact_v <- V(g_cross)$name[V(g_cross)$kingdom == "Bacteria"]
  g_sub <- induced_subgraph(g_cross, vids = bact_v)

  deg_b <- degree(g_bact, normalized = FALSE)
  btw_b <- betweenness(g_bact, normalized = TRUE)
  eig_b <- eigen_centrality(g_bact)$vector

  deg_c <- degree(g_sub, normalized = FALSE)
  btw_c <- betweenness(g_sub, normalized = TRUE)
  eig_c <- eigen_centrality(g_sub)$vector

  shared <- intersect(names(deg_b), names(deg_c))
  tibble::tibble(
    site = site,
    asv = shared,
    deg_bact_only = deg_b[shared],
    deg_cross = deg_c[shared],
    btw_bact_only = btw_b[shared],
    btw_cross = btw_c[shared],
    eig_bact_only = eig_b[shared],
    eig_cross = eig_c[shared]
  )
}

# centrality_shift <- function(g_bact, g_cross, bact_ids) {
#   shared <- intersect(V(g_bact)$name, V(g_cross)$name)
#   shared <- intersect(shared, bact_ids)

#   deg_b <- igraph::degree(g_bact)[shared]
#   deg_c <- igraph::degree(g_cross)[shared]
#   btw_b <- igraph::betweenness(g_bact, normalized = TRUE)[shared]
#   btw_c <- igraph::betweenness(g_cross, normalized = TRUE)[shared]
#   eig_b <- igraph::eigen_centrality(g_bact)$vector[shared]
#   eig_c <- igraph::eigen_centrality(g_cross)$vector[shared]

#   tibble(
#     node = shared,
#     delta_degree = deg_c - deg_b,
#     delta_betweenness = btw_c - btw_b,
#     delta_eigenvector = eig_c - eig_b
#   )
# }
