# Functions ----
#
# Workflow overview:
#
#   OTU tables
#       │
#       ▼
# [prev_filter()]       - prevalence-filter phyloseq objects
# [align_taxa()/align_samples()] - align phyloseq objects to shared taxa/samples
#       │
#       ▼
#   [build_network()]      - .cor2tran():correlation matrix -> list(graph, trans_mat) (single cutoff)
#                          - .add_attributes(): adds node attributes: degree, module, zi/pi roles. link attributes: kingdom (bacteria/fungi)

#   [build_asym_network()]- correlation matrix -> list(graph, trans_mat) (separate within/cross cutoffs)
#       │
#       ▼
# [network_properties()]   - calculates network-level properties: node/edge counts, density, transitivity, modularity, power law fit, etc.
# [network_summary()]
#   [centrality_shift()]   - compare node centrality: e.g.,bacteria-only vs cross-kingdom
#   [cyto_gephi_output()]  - export node/edge tables for Cytoscape / Gephi

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


# build_network() and build_asym_network() discussion

# They are not the same, they differ in two concrete ways:

# 1. What gets thresholded

# build_network(bxf_cor_matrices[[s]], ..., kind = "cross") applies a single uniform cutoff to the already-zeroed B×F matrix. Within-kingdom edges are structurally absent (zeroed upstream), not filtered by cutoff.

# build_asym_network(full_cor_matrices[[s]], ..., within_cutoff, bf_cutoff) applies two separate cutoffs to the full matrix, one for B×B and F×F blocks, another for the B×F block. Both within- and cross-kingdom edges are present and competing.

# 2. Node retention

# With bxf_cor_matrices, bacteria and fungi with no cross-kingdom correlations above cutoff are dropped as isolates. With build_asym_network(), a bacterium can be retained if it has strong within-kingdom edges even if all its B×F correlations fall below bf_cutoff.

# Equivalence depends on intent:

# Goal-> Function to use
# Pure B×F bipartite topology, single cutoff ->	build_network(bxf_cor_matrices, kind = "cross")
# Mixed topology — within + cross, asymmetric cutoffs ->	build_asym_network(full_cor_matrices, ...)
# Within-kingdom only ->	build_network(full_cor_matrices, kind = "bacteria"/"fungi")

# build_asym_network() would only differ if called on full_cor_matrices with within_cutoff != bf_cutoff.

# Three distinct matrix types, each producing a different network topology:

# Matrix types and network topologies:

# * full_cor_matrices — B×B ✓, F×F ✓, B×F ✓ -> Mixed: all within + cross edges
# * bxf_cor_matrices — B×B 0, F×F 0, B×F ✓ -> Strict bipartite: B×F edges only
# * build_asym_network(full_cor_matrices) — B×B ✓ within_cutoff, F×F ✓ within_cutoff, B×F ✓ bxf_cutoff -> Mixed with asymmetric cutoffs

# Practical implication for your analysis:

# ef_bxf_net -> use for bipartite ecological metrics (networklevel(), NODF)
# ef_bact_net / ef_fungi_net -> use for within-kingdom topology and RMT null model
# build_asym_network(full_cor_matrices) -> use when you want to ask how cross-kingdom edges restructure a network that also has within-kingdom edges

## Graph construction ----

# build_network: correlation matrix -> list(graph, trans_mat)
# Subsets to bacteria, fungi, or both kingdoms, applies a single |r| >= cutoff,
# drops isolated nodes, assigns edge sign and vertex kingdom.
# Optionally subsets to a focal neighbourhood (keep) or excludes nodes (rmv).
# Returns list: $graph = igraph object, $trans_mat = filtered correlation matrix.
.cor2tran <- function(cor_mat, cutoff, keep = c(), rmv = c()) {
  diag(cor_mat) <- 0 # Diagonal thresholding is redundant but ensures no self-loops if the input matrix isn't pre-zeroed on the diagonal.
  cor_mat[abs(cor_mat) < cutoff] <- 0 # Zero out weak correlations

  # drop isolated nodes
  cor_mat <- cor_mat[rowSums(cor_mat) != 0, colSums(cor_mat) != 0]

  # keep / rmv subsetting
  if (length(keep) > 0 && length(rmv) > 0) {
    stop("Provide either `keep` or `rmv`, not both.")
  }

  if (length(keep) > 0) {
    keep_idx <- which(rownames(cor_mat) %in% keep)
    if (length(keep_idx) == 0) {
      stop("None of the `keep` nodes found.")
    }
    nb_idx <- unique(c(
      keep_idx,
      which(colSums(cor_mat[, keep_idx, drop = FALSE]) != 0)
    ))
    cor_mat <- cor_mat[nb_idx, nb_idx]
  }

  if (length(rmv) > 0) {
    rmv_idx <- which(rownames(cor_mat) %in% rmv)
    if (length(rmv_idx) > 0) cor_mat <- cor_mat[-rmv_idx, -rmv_idx]
  }

  # re-drop isolates after subsetting
  cor_mat <- cor_mat[rowSums(cor_mat) != 0, colSums(cor_mat) != 0]

  return(cor_mat)
}


.add_attributes <- function(graph, tran_mat, bact_ids, fungi_ids) {
  if (igraph::vcount(graph) == 0) {
    cli::cli_alert_warning("Returning the same graph object. No vertices found")
    return(graph)
  }

  # Ensure vertex names exist — graph_from_adjacency_matrix only sets $name
  # if the input matrix has dimnames. Fall back to rownames(tran_mat).
  if (is.null(igraph::V(graph)$name)) {
    if (!is.null(rownames(tran_mat))) {
      igraph::V(graph)$name <- rownames(tran_mat)
    } else {
      cli::cli_alert_warning(
        ".add_attributes: graph has no vertex names and tran_mat has no rownames. ",
        "kingdom assignment will be NA."
      )
      igraph::V(graph)$name <- as.character(seq_len(igraph::vcount(graph)))
    }
  }

  # Link sign assignment logic (same as in add_link_attribute())
  if (igraph::ecount(graph) > 0) {
    el <- igraph::as_edgelist(graph)
    igraph::E(graph)$link_sign <- ifelse(
      tran_mat[el] > 0,
      "positive",
      "negative"
    )
  }

  igraph::V(graph)$kingdom <- dplyr::case_when(
    igraph::V(graph)$name %in% bact_ids ~ "Bacteria",
    igraph::V(graph)$name %in% fungi_ids ~ "Fungi",
    .default = NA_character_
  )

  # Node attribute assignment logic (same as in add_node_attribute())
  igraph::V(graph)$biomarker <- igraph::V(graph)$name
  igraph::V(graph)$node_degree <- centr_degree(graph)$res

  module_separation <- tryCatch(
    cluster_fast_greedy(graph),
    error = function(e) NULL
  )

  if (
    !is.null(module_separation) &&
      length(unique(membership(module_separation))) > 1
  ) {
    module_membership <- membership(module_separation)
    pi <- part_coeff(g = graph, memb = module_membership)
    zi <- within_module_deg_z_score(g = graph, memb = module_membership)
    role <- rep("peripherals", igraph::vcount(graph))
    role[which(pi >= 6.2 & zi < 2.5)] <- "connector"
    role[which(pi >= 6.2 & zi >= 2.5)] <- "network_hub"
    role[which(pi < 6.2 & zi >= 2.5)] <- "module_hub"
  } else {
    module_membership <- rep(NA_integer_, igraph::vcount(graph))
    pi <- rep(NA_real_, igraph::vcount(graph))
    zi <- rep(NA_real_, igraph::vcount(graph))
    role <- rep(NA_character_, igraph::vcount(graph))
  }

  igraph::V(graph)$module_membership <- as.integer(module_membership)
  igraph::V(graph)$pi <- pi
  igraph::V(graph)$zi <- zi
  igraph::V(graph)$vertex_role <- role

  return(graph)
}


build_network <- function(
  cor_mat,
  bact_ids,
  fungi_ids,
  cutoff = my_cutoff,
  kind = "cross",
  keep = c(), # focal node names - subset to their neighbourhood
  rmv = c() # node names to exclude
) {
  nodes <- rownames(cor_mat)
  b_ids <- intersect(bact_ids, nodes)
  f_ids <- intersect(fungi_ids, nodes)

  # For kind = "cross", cor_mat may be a rectangular B×F matrix.
  # Index rows (bacteria) and columns (fungi) separately, then symmetrise
  # into a square matrix so igraph can build an undirected graph.
  if (kind == "cross") {
    col_nodes <- colnames(cor_mat)
    f_ids_col <- intersect(fungi_ids, col_nodes)
    b_ids_row <- intersect(bact_ids, nodes)

    rect_mat <- cor_mat[b_ids_row, f_ids_col, drop = FALSE]

    # Symmetrise: build a square (B+F) × (B+F) matrix with B×F and F×B blocks
    all_ids <- c(b_ids_row, f_ids_col)
    n_all <- length(all_ids)
    sq_mat <- matrix(
      0,
      nrow = n_all,
      ncol = n_all,
      dimnames = list(all_ids, all_ids)
    )
    sq_mat[b_ids_row, f_ids_col] <- rect_mat
    sq_mat[f_ids_col, b_ids_row] <- t(rect_mat)

    sub_mat <- sq_mat
    b_ids <- b_ids_row
    f_ids <- f_ids_col
  } else {
    idx <- switch(kind, bacteria = b_ids, fungi = f_ids)
    sub_mat <- cor_mat[idx, idx]
  }

  # .cor2tran logic ---- applies cutoff, keeps/rmvs nodes, drops isolates
  sub_mat <- .cor2tran(
    cor_mat = sub_mat,
    cutoff = cutoff,
    keep = keep,
    rmv = rmv
  )

  # Make the graph (0/1 for absence/presence of edges, regardless of sign)
  adj <- (sub_mat != 0) * 1L
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    weighted = NULL
  )

  # .add_attributes logic ---- applies node and link attributes based on the filtered graph
  g <- .add_attributes(
    g,
    tran_mat = sub_mat,
    bact_ids = bact_ids,
    fungi_ids = fungi_ids
  )

  return(
    list(
      graph = g,
      trans_mat = sub_mat
    )
  )
}

# build_asym_network: correlation matrix -> list(graph, trans_mat) with asymmetric cutoffs
# Like build_network() but applies separate thresholds for within-kingdom edges
# (within_bact_cutoff and within_fungi_cutoff) and cross-kingdom bacteria–fungi edges (bxf_cutoff).
# Use this when B–F correlations are structurally weaker than within-kingdom ones.
build_asym_network <- function(
  cor_mat,
  bact_ids,
  fungi_ids,
  within_bact_cutoff,
  within_fungi_cutoff,
  bxf_cutoff
) {
  nodes <- rownames(cor_mat)
  b_ids <- intersect(bact_ids, nodes)
  f_ids <- intersect(fungi_ids, nodes)
  idx <- c(b_ids, f_ids)

  sub_mat <- cor_mat[idx, idx]
  diag(sub_mat) <- 0

  # apply asymmetric cutoffs via a logical mask
  mask <- matrix(
    FALSE,
    nrow = length(idx),
    ncol = length(idx),
    dimnames = list(idx, idx)
  )
  mask[b_ids, b_ids] <- abs(sub_mat[b_ids, b_ids]) >= within_bact_cutoff
  mask[f_ids, f_ids] <- abs(sub_mat[f_ids, f_ids]) >= within_fungi_cutoff
  mask[b_ids, f_ids] <- abs(sub_mat[b_ids, f_ids]) >= bxf_cutoff
  mask[f_ids, b_ids] <- t(mask[b_ids, f_ids])
  sub_mat[!mask] <- 0

  # # .cor2tran logic ---- applies cutoff, keeps/rmvs nodes, drops isolates
  # sub_mat <- .cor2tran(
  #   cor_mat = sub_mat,
  #   cutoff = cutoff,
  #   keep = keep,
  #   rmv = rmv
  # )

  # Make the graph (0/1 for absence/presence of edges, regardless of sign)
  adj <- (sub_mat != 0) * 1L
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    weighted = NULL
  )

  # .add_attributes logic ---- applies node and link attributes based on the filtered graph
  g <- .add_attributes(
    g,
    tran_mat = sub_mat,
    bact_ids = bact_ids,
    fungi_ids = fungi_ids
  )

  return(
    list(
      graph = g,
      trans_mat = sub_mat
    )
  )
}

## Network metrics ----
# network_summary: returns a one-row tibble of network-level metrics for a graph.
# Metrics: node/edge counts, density, transitivity, fraction of cross-kingdom
# edges, positive:negative edge ratio, modularity, and NODF nestedness
# (cross-kingdom networks only).

network_summary <- function(graph, site, kind) {
  n_v <- vcount(graph) #nodes
  n_e <- ecount(graph)
  n_pos <- sum(E(graph)$link_sign == "positive", na.rm = TRUE)
  n_neg <- sum(E(graph)$link_sign == "negative", na.rm = TRUE)

  # Edge count logic for cross-kingdom networks: count only edges between bacteria and fungi
  bf_edges <- if ("Fungi" %in% V(graph)$kingdom) {
    e_ends <- igraph::as_edgelist(graph)
    bact_v <- V(graph)$name[V(graph)$kingdom == "Bacteria"]
    fungi_v <- V(graph)$name[V(graph)$kingdom == "Fungi"]
    sum(
      (e_ends[, 1] %in% bact_v & e_ends[, 2] %in% fungi_v) |
        (e_ends[, 2] %in% bact_v & e_ends[, 1] %in% fungi_v)
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
  if (kind == "cross" && length(bact_v) > 1 && length(fungi_v) > 1) {
    adj_b <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
    bpt <- adj_b[bact_v, fungi_v, drop = FALSE]
    sub_matx_nodf <- tryCatch(
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
    # connectance = if (max_edges > 0) n_e / max_edges else NA_real_,
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
#'
centrality_shift <- function(g_bact, g_fungi, g_cross, site) {
  cross_b_v <- V(g_cross)$name[V(g_cross)$kingdom == "Bacteria"]
  cross_f_v <- V(g_cross)$name[V(g_cross)$kingdom == "Fungi"]

  # Subgraphs of the cross-kingdom graph for bacteria and fungi
  b_sub_cross <- induced_subgraph(g_cross, vids = cross_b_v)
  f_sub_cross <- induced_subgraph(g_cross, vids = cross_f_v)

  # Bacteria-only metrics
  deg_b <- degree(g_bact, normalized = FALSE)
  btw_b <- betweenness(g_bact, normalized = TRUE)
  eig_b <- eigen_centrality(g_bact)$vector

  deg_b_c <- degree(b_sub_cross, normalized = FALSE)
  btw_b_c <- betweenness(b_sub_cross, normalized = TRUE)
  eig_b_c <- eigen_centrality(b_sub_cross)$vector

  # Fungi-only metrics
  deg_f <- degree(g_fungi, normalized = FALSE)
  btw_f <- betweenness(g_fungi, normalized = TRUE)
  eig_f <- eigen_centrality(g_fungi)$vector

  deg_f_c <- degree(f_sub_cross, normalized = FALSE)
  btw_f_c <- betweenness(f_sub_cross, normalized = TRUE)
  eig_f_c <- eigen_centrality(f_sub_cross)$vector

  # Cross-kingdom centrality metrics
  deg_cr <- degree(g_cross, normalized = FALSE)
  btw_cr <- betweenness(g_cross, normalized = TRUE)
  eig_cr <- eigen_centrality(g_cross)$vector

  shared_b <- intersect(names(deg_b), names(deg_b_c))
  shared_f <- intersect(names(deg_f), names(deg_f_c))
  cross_names <- names(deg_cr)

  list(
    bacteria_metrics = tibble::tibble(
      site = site,
      asv = shared_b,
      deg_bact_only = deg_b[shared_b],
      deg_cross = deg_b_c[shared_b],
      btw_bact_only = btw_b[shared_b],
      btw_cross = btw_b_c[shared_b],
      eig_bact_only = eig_b[shared_b],
      eig_cross = eig_b_c[shared_b],
      delta_degree = deg_b_c[shared_b] - deg_b[shared_b],
      delta_betweenness = btw_b_c[shared_b] - btw_b[shared_b],
      delta_eigenvector = eig_b_c[shared_b] - eig_b[shared_b]
    ),
    fungi_metrics = tibble::tibble(
      site = site,
      asv = shared_f,
      deg_fungi_only = deg_f[shared_f],
      deg_cross = deg_f_c[shared_f],
      btw_fungi_only = btw_f[shared_f],
      btw_cross = btw_f_c[shared_f],
      eig_fungi_only = eig_f[shared_f],
      eig_cross = eig_f_c[shared_f],
      delta_degree = deg_f_c[shared_f] - deg_f[shared_f],
      delta_betweenness = btw_f_c[shared_f] - btw_f[shared_f],
      delta_eigenvector = eig_f_c[shared_f] - eig_f[shared_f]
    ),
    cross_network_centrality = tibble::tibble(
      site = site,
      asv = cross_names,
      deg_cross = deg_cr,
      btw_cross = btw_cr,
      eig_cross = eig_cr
    )
  )
}

## Export ----

# cyto_gephi_output: writes node and edge attribute tables for Cytoscape and
# Gephi. Requires a graph with all vertex and edge attributes already attached
# (via add_node_attribute() and add_link_attribute()).
cyto_gephi_output <- function(
  graph_with_all_attributes,
  file_path = "data/output/networks/",
  output_prefix = NULL
) {
  node_table <- as.data.frame(vertex.attributes(graph_with_all_attributes))
  edge_table <- as.data.frame(edge.attributes(graph_with_all_attributes))

  vnames_attr <- base::attr(igraph::E(graph_with_all_attributes), "vnames")

  if (!is.null(vnames_attr)) {
    edge_names <- unlist(strsplit(vnames_attr, "|", fixed = TRUE))
    edge_table$node1 <- edge_names[seq(
      from = 1,
      to = length(edge_names) - 1,
      by = 2
    )]
    edge_table$node2 <- edge_names[seq(
      from = 2,
      to = length(edge_names),
      by = 2
    )]
  } else {
    el <- igraph::as_edgelist(graph_with_all_attributes, names = TRUE)
    edge_table$node1 <- el[, 1]
    edge_table$node2 <- el[, 2]
  }

  gephi_node_table <- node_table
  if (ncol(gephi_node_table) > 0) {
    names(gephi_node_table)[1] <- "Id"
  }

  gephi_edge_table <- edge_table
  if (ncol(gephi_edge_table) >= 3) {
    names(gephi_edge_table)[2:3] <- c("Source", "Target")
  }
  gephi_edge_table$Type <- rep("Undirected", nrow(gephi_edge_table))

  if (is.null(file_path) || is.null(output_prefix)) {
    stop("file_path and output_prefix must be provided.")
  }
  if (!dir.exists(file_path)) {
    dir.create(file_path, recursive = TRUE)
  }

  write.table(
    node_table,
    paste0(
      file_path,
      output_prefix,
      "_cytoscape_node_attribute.txt"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    edge_table,
    paste0(
      file_path,
      output_prefix,
      "_cytoscape_edge_attribute.txt"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    gephi_node_table,
    paste0(file_path, output_prefix, "_gephi_node_attribute.csv"),
    sep = ",",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    gephi_edge_table,
    paste0(file_path, output_prefix, "_gephi_edge_attribute.csv"),
    sep = ",",
    row.names = FALSE,
    quote = FALSE
  )
}

## RMT cutoff ----

#' find_rmt_cutoff: Random Matrix Theory-based correlation threshold selection.
#'
#' Algorithm (Zhou et al. 2010 / Luo et al. 2006 — the MENA approach):
#'   For each candidate cutoff in `cutoff_seq`:
#'     1. Threshold `cor_mat`: zero entries where |r| < cutoff.
#'     2. Compute eigenvalues of the symmetric thresholded matrix.
#'     3. Unfold eigenvalues via a natural cubic spline fit to the empirical
#'        staircase function N(E).
#'        This maps the eigenvalues so that spacings have unit mean (removes
#'        the bulk shape).
#'     4. Compute Nearest Neighbor Spacing Distribution (NNSD):
#'        s_i = lambda_{i+1} - lambda_i of the unfolded eigenvalues.
#'     5. Bin NNSD into `n_bins` bins over s in [0, 3] and chi-square test
#'        against:
#'          - Poisson:       p(s) = exp(-s)           (random, no structure)
#'          - Wigner (Gaussian orthogonal ensemble -GOE):  p(s) = (pi/2)*s*exp
#'            (-pi/4*s^2)  (structured)
#'   Optimal cutoff: highest threshold at which Poisson is still *rejected*
#'   (p_poisson < alpha) — i.e., the network retains biological signal.
#'
#' @param cor_mat      Symmetric correlation matrix (output of full_cor_matrices
#'                     or bxf_cor_matrices). Diagonal is ignored.
#' @param bact_ids     Character vector of bacterial OTU IDs.
#' @param fungi_ids    Character vector of fungal OTU IDs.
#' @param kind         "bacteria", "fungi", or "cross" - which submatrix to analyze.
#'
#' @param cutoff_seq   Numeric vector of candidate thresholds to sweep,
#'                     e.g. seq(0.50, 0.10, by = -0.01).
#' @param n_bins       Number of histogram bins for NNSD (default 15).
#' @param alpha        Significance level for Poisson rejection (default 0.05).
#' @param verbose      Print progress for each cutoff step (default FALSE).
#'
#' @return A list with:
#'   $results        tibble — one row per cutoff: cutoff, n_eigenvalues,
#'                   chi2_poisson, p_poisson, chi2_wigner, p_wigner,
#'                   poisson_rejected, wigner_rejected.
#'   $optimal_cutoff scalar — highest cutoff where poisson_rejected is TRUE.
#'   $plot           ggplot — chi2 statistics vs. cutoff with transition marked.
#'
#' @examples
#' rmt <- find_rmt_cutoff(full_cor_matrices$ef)
#' rmt$optimal_cutoff
#' rmt$plot
#' rmt$results
#'

find_rmt_cutoff <- function(
  cor_mat,
  bact_ids = NULL,
  fungi_ids = NULL,
  kind = "cross",
  cutoff_seq = seq(0.50, 0.10, by = -0.01),
  n_bins = 15,
  alpha = 0.05,
  verbose = FALSE
) {
  mat_name <- deparse(substitute(cor_mat)) # capture name before any evaluation

  # 0. Helpers ----

  # Wigner-Dyson (GOE) density over bin midpoints
  .wigner_phf <- function(s) (pi / 2) * s * exp(-(pi / 4) * s^2)

  # Poisson density over bin midpoints
  .poisson_phf <- function(s) exp(-s)

  # Unfold eigenvalues via a natural cubic spline fit to the empirical staircase function N(E).
  .unfold <- function(eigs) {
    eigs <- sort(unique(eigs)) # Drop duplicates to avoid issues with spline interpolation
    n <- length(eigs)
    staircase <- seq_len(n) # N(E_i) = i

    # Natural cubic spline interpolating the staircase N_av(E)
    spline_fn <- stats::splinefun(x = eigs, y = staircase, method = "natural")

    # Unfolded spectrum: e_i = N_av(E_i)
    spline_fn(eigs)
  }

  # Chi-square GOF: observed vs. expected bin counts
  .chisq_gof <- function(spacings, phf_fn, breaks) {
    obs <- as.numeric(
      graphics::hist(spacings, breaks = breaks, plot = FALSE)$counts
    )
    mids <- (breaks[-length(breaks)] + breaks[-1]) / 2
    expected_probs <- phf_fn(mids)
    expected_probs <- expected_probs / sum(expected_probs)
    exp_counts <- expected_probs * length(spacings)

    # Merge sparse bins (Cochran's rule)
    merge_mask <- exp_counts >= 1
    if (sum(merge_mask) < 2) {
      return(list(chi2 = NA_real_, p = NA_real_))
    }

    obs <- obs[merge_mask]
    exp_counts <- exp_counts[merge_mask]

    chi2 <- sum((obs - exp_counts)^2 / exp_counts)
    df <- sum(merge_mask) - 1
    p <- stats::pchisq(chi2, df = df, lower.tail = FALSE)
    list(chi2 = chi2, p = p)
  }

  # 1. Prepare matrix ----
  cli::cli_h1("Finding RMT cutoff for kind='{kind}'")

  nodes <- rownames(cor_mat)

  if (!is.null(bact_ids)) {
    cli::cli_alert_info(
      "Filtering to {length(intersect(bact_ids, nodes))} bacterial IDs for kind='{kind}'."
    )
    b_ids <- intersect(bact_ids, nodes)
  } else {
    b_ids <- character(0)
  }
  if (!is.null(fungi_ids)) {
    cli::cli_alert_info(
      "Filtering to {length(intersect(fungi_ids, nodes))} fungal IDs for kind='{kind}'."
    )
    f_ids <- intersect(fungi_ids, nodes)
  } else {
    f_ids <- character(0)
  }

  idx <- switch(kind, bacteria = b_ids, fungi = f_ids, cross = c(b_ids, f_ids))

  mat0 <- cor_mat[idx, idx]

  if (is.null(bact_ids) && is.null(fungi_ids)) {
    cli::cli_alert_info(
      "No bacterial or fungal IDs provided. Analyzing the full matrix."
    )
    mat0 <- cor_mat
  }

  cli::cli_alert_info(
    "Analyzing a {nrow(mat0)}×{ncol(mat0)} matrix for kind='{kind}'."
  )
  diag(mat0) <- 0 # Avoid self-loops

  breaks <- seq(0, 3, length.out = n_bins + 1)

  # 2. Sweep cutoffs ----
  results_list <- vector("list", length(cutoff_seq))

  progressbar_calc_rmt <- cli::cli_progress_bar(
    name = "Calculating RMT cutoff",
    total = length(cutoff_seq),
    format = "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}",
    .auto_close = TRUE
  )
  for (i in seq_along(cutoff_seq)) {
    thr <- cutoff_seq[i]
    mat_thr <- mat0
    mat_thr[abs(mat_thr) < thr] <- 0

    # Drop isolated nodes — they add spurious zero eigenvalues
    keep <- rowSums(mat_thr != 0) > 0
    mat_thr <- mat_thr[keep, keep]
    n_keep <- nrow(mat_thr)

    if (n_keep < 10) {
      if (verbose) {
        cli::cli_inform(sprintf(
          "cutoff %.3f -> %d nodes, skipping",
          thr,
          n_keep
        ))
      }
      results_list[[i]] <- tibble::tibble(
        cutoff = thr,
        n_eigenvalues = n_keep,
        chi2_poisson = NA_real_,
        p_poisson = NA_real_,
        chi2_wigner = NA_real_,
        p_wigner = NA_real_,
        poisson_rejected = NA,
        wigner_rejected = NA
      )
      next
    }

    # Eigenvalues (all real for symmetric matrix)
    eigs <- eigen(mat_thr, symmetric = TRUE, only.values = TRUE)$values

    # Remove the dominant Perron–Frobenius eigenvalue before unfolding
    eigs <- eigs[-which.max(eigs)]

    # Unfold
    unfolded <- tryCatch(.unfold(eigs), error = function(e) NULL)
    if (is.null(unfolded) || length(unfolded) < 5) {
      results_list[[i]] <- tibble::tibble(
        cutoff = thr,
        n_eigenvalues = length(eigs),
        chi2_poisson = NA_real_,
        p_poisson = NA_real_,
        chi2_wigner = NA_real_,
        p_wigner = NA_real_,
        poisson_rejected = NA,
        wigner_rejected = NA
      )
      next
    }

    # Nearest-neighbour spacings, clipped to [0, 3]
    spacings <- diff(unfolded)
    spacings <- spacings[spacings > 0 & spacings <= 3]

    res_p <- .chisq_gof(spacings, .poisson_phf, breaks)
    res_w <- .chisq_gof(spacings, .wigner_phf, breaks)

    if (verbose) {
      cli::cli_inform(sprintf(
        "cutoff %.3f | n=%d | chi2_Poisson=%.2f (p=%.3f) | chi2_Wigner=%.2f (p=%.3f)",
        thr,
        length(eigs),
        res_p$chi2 %||% NA,
        res_p$p %||% NA,
        res_w$chi2 %||% NA,
        res_w$p %||% NA
      ))
    }

    results_list[[i]] <- tibble::tibble(
      cutoff = thr,
      n_eigenvalues = length(eigs),
      chi2_poisson = res_p$chi2,
      p_poisson = res_p$p,
      chi2_wigner = res_w$chi2,
      p_wigner = res_w$p,
      poisson_rejected = !is.na(res_p$p) & res_p$p < alpha,
      wigner_rejected = !is.na(res_w$p) & res_w$p < alpha
    )

    cli::cli_progress_update(id = progressbar_calc_rmt)
  }

  # Process alerts
  cli::cli_alert_warning(
    "Duplicate eigenvalues dropped to prevent issues with spline interpolation, \navoiding zero-spacing entries that can distort the NNSD."
  )

  results <- dplyr::bind_rows(results_list)

  # 3. Identify optimal cutoff ----
  # Highest cutoff where the NNSD still departs from Poisson = last cutoff
  # with biological signal before the network collapses to random.
  candidates <- results |> dplyr::filter(poisson_rejected == TRUE)

  if (nrow(candidates) == 0) {
    cli::cli_alert_warning(
      "find_rmt_cutoff: Poisson never rejected across the cutoff range. ",
      "The matrix may lack signal or cutoff_seq range is inappropriate. ",
      "Returning NA for optimal_cutoff."
    )
    optimal <- NA_real_
  } else {
    optimal <- max(candidates$cutoff)
  }

  cli::cli_alert_success(
    "RMT cutoff analysis complete. Optimal cutoff: {round(optimal, 3)}"
  )

  # 4. Diagnostic chi-squared plot ----
  cli::cli_h2("RMT cutoff diagnostic plots")
  cli::cli_alert_info(
    "Generating RMT cutoff diagnostic plot: chi-squared vs. cutoff with optimal marked."
  )

  plot_df <- results |>
    dplyr::filter(!is.na(chi2_poisson)) |>
    tidyr::pivot_longer(
      cols = c(chi2_poisson, chi2_wigner),
      names_to = "distribution",
      values_to = "chi2"
    ) |>
    dplyr::mutate(
      distribution = dplyr::recode(
        distribution,
        chi2_poisson = "Poisson",
        chi2_wigner = "Wigner-Dyson"
      )
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = cutoff, y = chi2, color = distribution)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(
      values = c("Poisson" = "#E05C5C", "Wigner-Dyson" = "#4C8CBF"),
      name = "Reference\ndistribution"
    ) +
    ggplot2::geom_vline(
      xintercept = optimal,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.7
    ) +
    ggplot2::annotate(
      "text",
      x = optimal,
      y = max(plot_df$chi2, na.rm = TRUE) * 0.95,
      label = paste0("optimal\n", round(optimal, 3)),
      hjust = -0.1,
      size = 3,
      color = "grey20"
    ) +
    ggplot2::labs(
      title = paste0("RMT cutoff selection: ", mat_name, " (", kind, ")"),
      subtitle = paste0(
        "Optimal = highest cutoff where Poisson is rejected (p < ",
        alpha,
        ")"
      ),
      x = "Correlation threshold",
      y = expression(chi^2 ~ statistic)
    ) +
    ggplot2::theme_bw()

  # 5. NNSD plot: mirroring Luo (2006) Fig. 1 ----

  cli::cli_alert_info(
    "Generating NNSD plot for optimal cutoff: observed vs. Poisson and Wigner-Dyson."
  )

  s_grid <- seq(0, 3, length.out = 300)

  theory_df <- dplyr::bind_rows(
    tibble::tibble(
      s = s_grid,
      density = exp(-s_grid),
      distribution = "Poisson"
    ),
    tibble::tibble(
      s = s_grid,
      density = (pi / 2) * s_grid * exp(-(pi / 4) * s_grid^2),
      distribution = "Wigner-Dyson"
    )
  )
  # Rescale theoretical curves to relative frequency over the same bin grid
  # by integrating each pdf over bin width so they sum to ~1.

  theory_rf <- theory_df |>
    dplyr::mutate(density = density * 1)

  # Build tidy DF: one row per bin midpoint per cutoff.
  # Luo (2006) Fig. 1.

  h <- graphics::hist(spacings, breaks = breaks, plot = FALSE)
  nnsd_df <- tibble::tibble(
    s = h$mids,
    density = h$counts / sum(h$counts)
  )

  nnsd_p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = nnsd_df,
      ggplot2::aes(x = s, y = density),
      color = "red",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = nnsd_df,
      ggplot2::aes(x = s, y = density),
      color = "red",
      size = 2
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(theory_rf, distribution == "Wigner-Dyson"),
      ggplot2::aes(x = s, y = density),
      color = "black",
      linetype = "solid",
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(theory_rf, distribution == "Poisson"),
      ggplot2::aes(x = s, y = density),
      color = "black",
      linetype = "dashed",
      linewidth = 0.9
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
    ggplot2::coord_cartesian(xlim = c(0, 3), ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0(
        "NNSD: ",
        mat_name,
        " (",
        kind,
        ") - cutoff = ",
        round(optimal, 3)
      ),
      subtitle = "Black solid = Wigner-Dyson (GOE); black dashed = Poisson",
      x = "Level spacing  s",
      y = "P(s)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "right")

  list(
    results = results,
    optimal_cutoff = optimal,
    plot = p,
    nnsd_plot = nnsd_p
  )
}
