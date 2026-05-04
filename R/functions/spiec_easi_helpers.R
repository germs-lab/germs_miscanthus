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
  pulsar_params <- list(
    rep.num = rep.num,
    ncores = ncores,
    seed = 10010
  )
  if (is.null(ps2)) {
    spiec.easi(
      data = ps1,
      method = method,
      nlambda = nlambda,
      verbose = TRUE,
      lambda.min.ratio = lambda.min.ratio,
      pulsar.params = pulsar_params
    )
  } else {
    multi.spiec.easi(
      datalist = list(ps1, ps2),
      method = method,
      nlambda = nlambda,
      verbose = TRUE,
      lambda.min.ratio = lambda.min.ratio,
      pulsar.params = pulsar_params
    )
  }
}


plot_igraph <- function(
  fit,
  aligned_matrix1,
  aligned_matrix2 = NULL,
  adj_matrix = "StARS-refit",
  type,
  color,
  method = "mb",
  mode = "maxabs",
  label = NULL,
  main_title = NULL,
  kingdom = c("Bacteria", "Fungi")
) {
  # Vertex names
  vnames_1 <- colnames(aligned_matrix1) #samples must be columns
  vnames_2 <- if (!is.null(aligned_matrix2)) {
    colnames(aligned_matrix2)
  } else {
    character(0)
  }
  if (!is.null(aligned_matrix2)) {
    vnames <- c(vnames_1, vnames_2)
  } else {
    vnames <- vnames_1
  }

  # Kingdom assignment
  kingdom_attr <- if (!is.null(aligned_matrix2)) {
    c(rep("Bacteria", length(vnames_1)), rep("Fungi", length(vnames_2)))
  } else {
    rep(kingdom[1], length(vnames_1))
  }

  # Stars refit matrix (binary adjacency)
  stars_refit <- getRefit(fit)

  if (method == "mb") {
    opt_mat <- symBeta(getOptBeta(fit), mode = mode)
    diag(opt_mat) <- 0
  } else if (method == "glasso") {
    opt_mat <- symBeta(getOptCov(fit), mode = mode)
    diag(opt_mat) <- 0
  } else {
    stop("Unsupported method: ", method)
  }

  # Adjacency matrix for igraph construction
  if (adj_matrix == "StARS-refit") {
    adj <- getRefit(fit)
  } else if (adj_matrix == "optMat") {
    adj <- opt_mat
  } else {
    stop("Unsupported adj_matrix option: ", adj_matrix)
  }

  # Build igraph object
  ig_obj <- adj2igraph(
    adj,
    vertex.attr = list(
      name = vnames,
      kingdom = kingdom_attr
    )
  )

  if (igraph::ecount(ig_obj) == 0) {
    warning("Network has no edges; returning NULL.")
    return(NULL)
  }

  # pull edge signs from the weighted beta matrix
  edge_idx <- igraph::as_edgelist(ig_obj, names = FALSE)
  E(ig_obj)$sign <- sign(as.numeric(opt_mat[edge_idx]))
  E(ig_obj)$weight <- abs(as.numeric(opt_mat[edge_idx]))

  # Plot with phyloseq's plot_network
  vsize <- rowMeans(clr(aligned_matrix1, 1)) + 6
  am.coord <- layout_with_fr(ig_obj)

  p1 <- ggraph(ig_obj, layout = am.coord) +
    geom_edge_link(aes(color = sign), alpha = 0.5) +
    geom_node_point(aes(color = kingdom), size = vsize) +
    labs(title = main_title) +
    theme_graph()

  if (!is.null(aligned_matrix2)) {
    p2 <- plot(
      ig_obj,
      layout = am.coord,
      vertex.size = vsize,
      vertex.label = NA,
      main = main_title
    )
    p <- cowplot::plot_grid(p1, p2, ncol = 2)
  } else {
    p <- p1
  }

  return(list(
    igraph_obj = ig_obj,
    plot = p
  ))
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
