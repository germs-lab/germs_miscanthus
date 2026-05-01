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

#' Run SpiecEasi on one (bacteria-only / fungi-only) or two (cross-kingdom)
#' phyloseq objects.
#' @param ps1     primary phyloseq (bacteria)
#' @param ps2     optional second phyloseq (fungi); triggers multi.spiec.easi
#' @param method  "mb" (neighbourhood selection) or "glasso"
run_spiec <- function(
  ps1,
  ps2 = NULL,
  method = "mb",
  nlambda = 20,
  lambda.min.ratio = 1e-2,
  rep.num = 20,
  ncores = 4
) {
  pulsar_params <- list(rep.num = rep.num, ncores = ncores)
  if (is.null(ps2)) {
    spiec.easi(
      data = ps1,
      method = method,
      nlambda = nlambda,
      lambda.min.ratio = lambda.min.ratio,
      pulsar.params = pulsar_params
    )
  } else {
    multi.spiec.easi(
      datalist = list(ps1, ps2),
      method = method,
      nlambda = nlambda,
      lambda.min.ratio = lambda.min.ratio,
      pulsar.params = pulsar_params
    )
  }
}

#' Convert a SpiecEasi fit to a signed, weighted igraph object.
#' Vertex attribute `kingdom` is set to "Bacteria" or "Fungi".
#' Edge attribute `sign` encodes positive (+1) / negative (-1) association.
#'

spiec_to_igraph <- function(
  fit,
  aligned_ps1,
  method = "mb",
  mode = "maxabs",
  aligned_ps2 = NULL
) {
  stars_refit <- getRefit(fit)
  if (method == "mb") {
    opt_mat <- symBeta(getOptBeta(fit), mode = mode)
  } else if (method == "glasso") {
    opt_mat <- symBeta(getOptCov(fit), mode = mode)
  } else {
    stop("Unsupported method: ", method)
  }

  diag(opt_mat) <- 0

  g <- graph_from_adjacency_matrix(
    opt_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  bact_ids <- taxa_names(aligned_ps1)
  if (!is.null(aligned_ps2)) {
    fung_ids <- taxa_names(aligned_ps2)
    V(g)$name <- c(bact_ids, fung_ids)
    V(g)$kingdom <- c(
      rep("Bacteria", length(bact_ids)),
      rep("Fungi", length(fung_ids))
    )
  } else {
    V(g)$name <- bact_ids
    V(g)$kingdom <- rep("Bacteria", length(bact_ids))
  }
  # E(g)$sign <- sign(E(g)$weight)
  E(g)$sign <- sign(as.numeric(E(g)$weight))
  g
}


plot_igraph <- function(
  fit,
  aligned_ps1,
  source_matrix,
  type,
  color,
  method = "mb",
  mode = "maxabs",
  label = NULL,
  kingdom = "Bacteria"
) {
  #TODO: this function is a bit redundant with spiec_to_igraph; consider merging them and just returning the igraph object, then plotting separately. The main reason for keeping them separate for now is that this function can also return the ggplot object from plot_network, which is more flexible for downstream customization.
  vnames <- taxa_names(aligned_ps1)
  stars_refit <- getRefit(fit)
  if (method == "mb") {
    opt_mat <- symBeta(getOptBeta(fit), mode = mode)
  } else if (method == "glasso") {
    opt_mat <- symBeta(getOptCov(fit), mode = mode)
  } else {
    stop("Unsupported method: ", method)
  }

  ig_obj <- adj2igraph(
    getRefit(fit),
    vertex.attr = list(
      name = vnames,
      kingdom = rep(kingdom, length(vnames))
    )
  )

  if (igraph::ecount(ig_obj) == 0) {
    warning("Network has no edges; returning NULL.")
    return(NULL)
  }

  # pull edge signs from the weighted beta matrix
  edge_idx <- igraph::get.edgelist(ig_obj, names = FALSE)
  E(ig_obj)$sign <- sign(as.numeric(opt_mat[edge_idx]))

  p <- phyloseq::plot_network(
    ig_obj,
    aligned_ps1,
    type = type,
    color = color,
    label = label
  )
  return(p)
}

#' Network-level summary statistics for one igraph object.
net_summary <- function(g, site, type) {
  n_v <- vcount(g)
  n_e <- ecount(g)
  n_pos <- sum(E(g)$sign == 1, na.rm = TRUE)
  n_neg <- sum(E(g)$sign == -1, na.rm = TRUE)

  bf_edges <- if ("Fungi" %in% V(g)$kingdom) {
    e_ends <- ends(g, E(g))
    bact_v <- V(g)$name[V(g)$kingdom == "Bacteria"]
    fung_v <- V(g)$name[V(g)$kingdom == "Fungi"]
    sum(
      (e_ends[, 1] %in% bact_v & e_ends[, 2] %in% fung_v) |
        (e_ends[, 2] %in% bact_v & e_ends[, 1] %in% fung_v)
    )
  } else {
    NA_integer_
  }

  max_edges <- n_v * (n_v - 1) / 2
  wc <- cluster_louvain(g)

  nestedness_val <- tryCatch(
    {
      adj <- as_adjacency_matrix(g, sparse = FALSE)
      vegan::nestednodf(adj)$statistic[["NODF"]]
    },
    error = function(e) NA_real_
  )

  tibble::tibble(
    site = site,
    network_type = type,
    n_taxa = n_v,
    n_edges = n_e,
    density = graph.density(g),
    connectance = if (max_edges > 0) n_e / max_edges else NA_real_,
    bf_edge_fraction = if (!is.na(bf_edges)) bf_edges / n_e else NA_real_,
    pos_neg_ratio = if (n_neg > 0) n_pos / n_neg else NA_real_,
    transitivity = transitivity(g, type = "global"),
    modularity = modularity(wc),
    nestedness_NODF = nestedness_val
  )
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
