##########################################################################

# Exploring networks with

# "Fungal-Bacterial Cooccurrence Patterns Differ between Arbuscular Mycorrhizal Fungi and Nonmycorrhizal Fungi across Soil Niches". Citation: Yuan MM, Kakouridis A, Starr E, Nguyen N, Shi S, Pett-Ridge J, Nuccio E, Zhou J, Firestone M. 2021. Fungal-bacterial cooccurrence patterns differ between arbuscularmycorrhizal fungi and nonmycorrhizal fungi across soil niches. mBio 12:e03509-20. https://doi.org/10.1128/mBio.03509-20.

# This code calculates the full correlation matrix for constructing the network from OTU tables.

##########################################################################

#  Load data ----

source("R/utils/00_setup.R")
library("igraph")
library("ggraph")
library(brainGraph)
library(bipartite)

source("R/functions/rmt_approach/fit_power_law.R")
source("R/functions/rmt_approach/global.R")
source("R/functions/rmt_approach/bpt.R")
source("R/functions/rmt_approach/rand_adj_gen.R")

# Functions ----
#
# Workflow overview:
#
#   OTU tables
#       │
#       ▼
#   [build_network()]      — correlation matrix → list(graph, trans_mat) (single cutoff)
#   [build_cross_network()]— correlation matrix → list(graph, trans_mat) (separate within/cross cutoffs)
#       │
#       ▼
#   [add_node_attribute()] — attach degree, module, zi/pi roles to vertices
#       │
#       ▼
#   [net_summary()]        — tidy tibble of network-level metrics
#   [centrality_shift()]   — compare node centrality: bacteria-only vs cross-kingdom
#   [cyto_gephi_output()]  — export node/edge tables for Cytoscape / Gephi

## Graph construction ----

# build_network: correlation matrix → list(graph, trans_mat)
# Subsets to bacteria, fungi, or both kingdoms, applies a single |r| >= cutoff,
# drops isolated nodes, assigns edge sign and vertex kingdom.
# Optionally subsets to a focal neighbourhood (keep) or excludes nodes (rmv).
# Returns list: $graph = igraph object, $trans_mat = filtered correlation matrix.
build_network <- function(
  cor_mat,
  bact_ids,
  fungi_ids,
  cutoff = my_cutoff,
  kind = "cross",
  keep = character(0), # focal node names — subset to their neighbourhood
  rmv = character(0) # node names to exclude
) {
  nodes <- rownames(cor_mat)
  b_ids <- intersect(bact_ids, nodes)
  f_ids <- intersect(fungi_ids, nodes)

  idx <- switch(kind, bacteria = b_ids, fungi = f_ids, cross = c(b_ids, f_ids))

  sub_mat <- cor_mat[idx, idx]

  # cor2tran logic ----
  diag(sub_mat) <- 0
  sub_mat[abs(sub_mat) < cutoff] <- 0

  # drop isolated nodes
  sub_mat <- sub_mat[rowSums(sub_mat) != 0, colSums(sub_mat) != 0]

  # keep / rmv subsetting
  if (length(keep) > 0 && length(rmv) > 0) {
    stop("Provide either `keep` or `rmv`, not both.")
  }

  if (length(keep) > 0) {
    keep_idx <- which(rownames(sub_mat) %in% keep)
    if (length(keep_idx) == 0) {
      stop("None of the `keep` nodes found.")
    }
    nb_idx <- unique(c(
      keep_idx,
      which(colSums(sub_mat[, keep_idx, drop = FALSE]) != 0)
    ))
    sub_mat <- sub_mat[nb_idx, nb_idx]
  }

  if (length(rmv) > 0) {
    rmv_idx <- which(rownames(sub_mat) %in% rmv)
    if (length(rmv_idx) > 0) sub_mat <- sub_mat[-rmv_idx, -rmv_idx]
  }

  # re-drop isolates after subsetting
  sub_mat <- sub_mat[rowSums(sub_mat) != 0, colSums(sub_mat) != 0]
  # ----

  adj <- (sub_mat != 0) * 1L
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    weighted = NULL
  )

  if (igraph::ecount(g) > 0) {
    el <- igraph::as_edgelist(g)
    igraph::E(g)$link_sign <- ifelse(sub_mat[el] > 0, "positive", "negative")
  }

  igraph::V(g)$kingdom <- ifelse(
    igraph::V(g)$name %in% bact_ids,
    "Bacteria",
    "Fungi"
  )

  list(graph = g, trans_mat = sub_mat)
}

# build_cross_network: correlation matrix → list(graph, trans_mat) with asymmetric cutoffs
# Like build_network() but applies separate thresholds for within-kingdom edges
# (within_cutoff) and cross-kingdom bacteria–fungi edges (bf_cutoff).
# Use this when B–F correlations are structurally weaker than within-kingdom ones.
build_cross_network <- function(
  cor_mat,
  bact_ids,
  fungi_ids,
  within_cutoff,
  bf_cutoff
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
  mask[b_ids, b_ids] <- abs(sub_mat[b_ids, b_ids]) >= within_cutoff
  mask[f_ids, f_ids] <- abs(sub_mat[f_ids, f_ids]) >= within_cutoff
  mask[b_ids, f_ids] <- abs(sub_mat[b_ids, f_ids]) >= bf_cutoff
  mask[f_ids, b_ids] <- t(mask[b_ids, f_ids])
  sub_mat[!mask] <- 0

  # drop isolated nodes
  sub_mat <- sub_mat[rowSums(sub_mat) != 0, colSums(sub_mat) != 0]

  adj <- (sub_mat != 0) * 1L
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    weighted = NULL
  )

  # Graph annotation logic (same as in add_node_attribute() and add_link_attribute())

  # Link sign assignment logic (same as in add_link_attribute())
  if (igraph::ecount(g) > 0) {
    el <- igraph::as_edgelist(g)
    igraph::E(g)$link_sign <- ifelse(sub_mat[el] > 0, "positive", "negative")
  }

  igraph::V(g)$kingdom <- ifelse(
    igraph::V(g)$name %in% bact_ids,
    "Bacteria",
    "Fungi"
  )

  # Node attribute assignment logic (same as in add_node_attribute())
  igraph::V(g)$biomarker <- igraph::V(g)$name
  igraph::V(g)$node_degree <- centr_degree(g)$res

  module_separation <- cluster_fast_greedy(g)
  module_membership <- membership(module_separation)
  pi <- part_coeff(g = g, memb = module_membership)
  zi <- within_module_deg_z_score(g = g, memb = module_membership)
  role <- rep("peripherals", igraph::vcount(g))
  role[which(pi >= 6.2 & zi < 2.5)] <- "connector"
  role[which(pi >= 6.2 & zi >= 2.5)] <- "network_hub"
  role[which(pi < 6.2 & zi >= 2.5)] <- "module_hub"

  igraph::V(g)$module_membership <- module_membership
  igraph::V(g)$pi <- pi
  igraph::V(g)$zi <- zi
  igraph::V(g)$vertex_role <- role

  return(
    list(
      graph = g,
      trans_mat = sub_mat
    )
  )
}

## Graph annotation ----

# # add_node_attribute: attaches computed node-level metrics as vertex attributes.
# # Adds: biomarker label, degree, module membership (fast-greedy), participation
# # coefficient (pi), within-module degree z-score (zi), and topological role
# # (peripheral / connector / module hub / network hub).
# add_node_attribute <- function(graph) {
#   biomarker <- gsub("OTU_.*_", "", vertex_attr(graph)$name)

#   OTU_name <- vertex_attr(graph)$name
#   node_degree <- centr_degree(graph)$res

#   module_separation <- cluster_fast_greedy(graph)
#   module_membership <- membership(module_separation)
#   pi <- part_coeff(g = graph, memb = module_membership)
#   zi <- within_module_deg_z_score(g = graph, memb = module_membership)
#   role <- rep("peripherals", gorder(graph))
#   role[which(pi >= 6.2 & zi < 2.5)] <- "connector"
#   role[which(pi >= 6.2 & zi >= 2.5)] <- "network_hub"
#   role[which(pi < 6.2 & zi >= 2.5)] <- "module_hub"

#   graph <- set_vertex_attr(graph, "biomarker", index = V(graph), biomarker)
#   graph <- set_vertex_attr(graph, "node_degree", index = V(graph), node_degree)
#   graph <- set_vertex_attr(
#     graph,
#     "module_membership",
#     index = V(graph),
#     module_membership
#   )
#   graph <- set_vertex_attr(graph, "pi", index = V(graph), pi)
#   graph <- set_vertex_attr(graph, "zi", index = V(graph), zi)
#   graph <- set_vertex_attr(graph, "vertex_role", index = V(graph), role)

#   return(graph)
# }

# # add_link_attribute: attaches positive/negative sign to each edge.
# # Reads sign from a transition matrix (values > 0 → "positive", < 0 →
# # "negative") and stores it as the link_sign edge attribute.
# # Note: build_network() and build_cross_network() assign link_sign internally,
# # so this function is only needed when working with graphs built via cor2tran()
# # + graph_from_adjacency_matrix().
# add_link_attribute <- function(graph, tranmatx) {
#   mat_p_n <- tranmatx
#   mat_p_n[mat_p_n > 0] <- "positive"
#   mat_p_n[mat_p_n < 0] <- "negative"

#   edge_names <- unlist(strsplit(
#     base::attr(E(graph), "vnames"),
#     "|",
#     fixed = T
#   ))
#   edge_rownames <- edge_names[seq(1, (length(edge_names) - 1), by = 2)]
#   edge_colnames <- edge_names[seq(2, length(edge_names), by = 2)]
#   edge_dir <- unlist(lapply(
#     c(1:length(edge_rownames)),
#     FUN = function(id, np_matrix, rownames, colnames) {
#       dir <- np_matrix[
#         which(rownames(np_matrix) == edge_rownames[id]),
#         which(colnames(np_matrix) == edge_colnames[id])
#       ]
#     },
#     np_matrix = mat_p_n,
#     rownames = edge_rownames,
#     colnames = edge_colnames
#   ))

#   graph <- set_edge_attr(graph, "link_sign", index = E(graph), edge_dir)

#   return(graph)
# }

## Network metrics ----

# net_summary: returns a one-row tibble of network-level metrics for a graph.
# Metrics: node/edge counts, density, transitivity, fraction of cross-kingdom
# edges, positive:negative edge ratio, modularity, and NODF nestedness
# (cross-kingdom networks only).
net_summary <- function(g, bact_ids, fungi_ids, kind) {
  n <- igraph::vcount(g)
  ec <- igraph::ecount(g)
  el <- igraph::as_edgelist(g)
  b <- intersect(V(g)$name, bact_ids)
  f <- intersect(V(g)$name, fungi_ids)
  bf_e <- sum(
    (el[, 1] %in% b & el[, 2] %in% f) | (el[, 1] %in% f & el[, 2] %in% b)
  )
  pos <- sum(E(g)$link_sign == "positive")
  cl <- igraph::cluster_louvain(g)

  nodf <- NA_real_
  if (kind == "cross" && length(b) > 1 && length(f) > 1) {
    adj_b <- igraph::as_adjacency_matrix(g, sparse = FALSE)
    bpt <- adj_b[b, f]
    nodf <- tryCatch(
      bipartite::networklevel(bpt, index = "NODF")[["NODF"]],
      error = function(e) NA_real_
    )
  }

  tibble(
    kind = kind,
    nodes = n,
    edges = ec,
    density = round(igraph::edge_density(g), 4),
    transitivity = round(igraph::transitivity(g, type = "global"), 3),
    frac_BF_edges = round(bf_e / max(ec, 1), 3),
    pos_neg_ratio = round(pos / max(ec - pos, 1), 2),
    modularity = round(igraph::modularity(cl), 3),
    NODF_nestedness = round(nodf, 2)
  )
}

# centrality_shift: compares node centrality for bacteria shared between a
# bacteria-only graph and a cross-kingdom graph.
# Returns a tibble with delta_degree, delta_betweenness, delta_eigenvector
# for each shared bacterium. Positive deltas indicate centrality gained when
# fungi are included in the network.
centrality_shift <- function(g_bact, g_cross, bact_ids) {
  shared <- intersect(V(g_bact)$name, V(g_cross)$name)
  shared <- intersect(shared, bact_ids)

  deg_b <- igraph::degree(g_bact)[shared]
  deg_c <- igraph::degree(g_cross)[shared]
  btw_b <- igraph::betweenness(g_bact, normalized = TRUE)[shared]
  btw_c <- igraph::betweenness(g_cross, normalized = TRUE)[shared]
  eig_b <- igraph::eigen_centrality(g_bact)$vector[shared]
  eig_c <- igraph::eigen_centrality(g_cross)$vector[shared]

  tibble(
    node = shared,
    delta_degree = deg_c - deg_b,
    delta_betweenness = btw_c - btw_b,
    delta_eigenvector = eig_c - eig_b
  )
}

## Export ----

# cyto_gephi_output: writes node and edge attribute tables for Cytoscape and
# Gephi. Requires a graph with all vertex and edge attributes already attached
# (via add_node_attribute() and add_link_attribute()).
cyto_gephi_output <- function(graph_with_all_attributes) {
  node_table <- as.data.frame(vertex.attributes(graph_with_all_attributes))

  edge_table <- as.data.frame(edge.attributes(graph_with_all_attributes))
  edge_names <- unlist(strsplit(
    base::attr(E(graph_with_all_attributes), "vnames"),
    "|",
    fixed = T
  ))
  edge_table$node1 <- edge_names[seq(
    from = 1,
    to = (length(edge_names) - 1),
    by = 2
  )]
  edge_table$node2 <- edge_names[seq(from = 2, to = length(edge_names), by = 2)]

  gephi_node_table <- node_table
  names(gephi_node_table)[1] <- "Id"

  gephi_edge_table <- edge_table
  names(gephi_edge_table)[2:3] <- c("Source", "Target")
  gephi_edge_table$Type = rep("Undirected", nrow(gephi_edge_table))

  write.table(
    node_table,
    "data/output/networks/cytoscape_node_attribute.txt",
    sep = "\t",
    row.names = F,
    quote = F
  )
  write.table(
    edge_table,
    "data/output/networks/cytoscape_edge_attribute.txt",
    sep = "\t",
    row.names = F,
    quote = F
  )
  write.table(
    gephi_node_table,
    "data/output/networks/gephi_node_attribute.csv",
    sep = ",",
    row.names = F,
    quote = F
  )
  write.table(
    gephi_edge_table,
    "data/output/networks/gephi_edge_attribute.csv",
    sep = ",",
    row.names = F,
    quote = F
  )
}

# Pre-process ----

main_mxg_physeq_list <- subset_to_miscanthus(
  main_physeq_list,
  crop_patterns = c("MXG", "M", "Miscanthus"),
  filter_DNA_only = TRUE
)


sites <- names(main_mxg_physeq_list)

# Map each site name to its bacteria and fungi phyloseq objects
site_bact_list <- list(
  ef = main_mxg_physeq_list$ef$ef_16S_DNA,
  lamps_2018 = main_mxg_physeq_list$lamps_2018$lamps_2018_16S_DNA,
  lamps_2022 = main_mxg_physeq_list$lamps_2022$lamps_2022_16S_DNA
)

site_fungi_list <- list(
  ef = main_mxg_physeq_list$ef$ef_AMF_DNA,
  lamps_2018 = main_mxg_physeq_list$lamps_2018$lamps_2018_ITS_DNA,
  lamps_2022 = main_mxg_physeq_list$lamps_2022$lamps_2022_AMF_DNA
)

# Prevalence-filter each site independently
bact_filtered <- purrr::map(
  site_bact_list,
  ~ prev_filter(.x, prev_thresh = 0.2)
)
fungi_filtered <- purrr::map(
  site_fungi_list,
  ~ prev_filter(.x, prev_thresh = 0.2)
)

# Joining
joint_filtered <- purrr::map2(
  bact_filtered,
  fungi_filtered,
  ~ merge_phyloseq(.x, .y)
)

# Per site: pair bacteria and fungi phyloseq objects without requiring shared
# sample names.
aligned_matrices <- purrr::map(sites, function(s) {
  extract_matrices <- function(ps) {
    mat <- as.matrix(as.data.frame(otu_table(ps)))
    if (!taxa_are_rows(ps)) {
      #  with samples on columns and features/OTUs in rows
      mat <- t(mat)
    }
    mat
  }

  ps_b <- bact_filtered[[s]]
  ps_f <- fungi_filtered[[s]]
  ps_j <- joint_filtered[[s]]

  # Sort samples independently — SpiecEasi will use positional correspondence
  ps_b <- prune_samples(sort(sample_names(ps_b)), ps_b)
  ps_f <- prune_samples(sort(sample_names(ps_f)), ps_f)
  ps_j <- prune_samples(sort(sample_names(ps_j)), ps_j)

  ps_b_mat <- extract_matrices(ps_b)
  ps_f_mat <- extract_matrices(ps_f)
  ps_j_mat <- extract_matrices(ps_j)

  message(
    s,
    " – bacteria samples: ",
    nsamples(ps_b),
    " | fungi samples: ",
    nsamples(ps_f),
    " | joint samples: ",
    nsamples(ps_j),
    " | bacteria ASVs: ",
    ntaxa(ps_b),
    " | fungi ASVs: ",
    ntaxa(ps_f),
    " | joint ASVs: ",
    ntaxa(ps_j)
  )

  list(matx_b = ps_b_mat, matx_f = ps_f_mat, matx_j = ps_j_mat)
}) |>
  set_names(sites)

# Majority_corematrix ----

joint_cor_matrices <- purrr::imap(
  aligned_matrices,
  function(site_matx, site_name) {
    otu_matx <- aligned_matrices[[site_name]][["matx_j"]] #joint matrix with bacteria and fungi
    otu_matx[is.na(otu_matx)] <- 0

    ## choose correlation methods from: 1) Spearman; 2) Pearson; 3) central log-ratio Pearson
    # 1) Spearman
    # otu_matx_cor = cor(t(otu_matx), method = "spearman")

    # # 2) Pearson, used in this study
    # otu_matx_cor = cor(t(otu_matx), method = "pearson")

    # 3) central log-ratio Pearson
    clr <- apply((otu_matx + 1), 2, function(xc) {
      log(xc, 2) - sum(log(xc, 2)) / length(xc)
    })
    cor(t(clr), method = "pearson")
  }
)


save(joint_cor_matrices, file = "data/output/networks/correlation_matrices.rda")


# Cormatrix for bipartite network ----
load("data/output/networks/correlation_matrices.rda")

## Bacteria and Fungi OTUs IDs in the correlation matrix ----
# Since we know where they came from, we can directly get the id of 16S and ITS OTUs using the rownames of the aligned matrices.

bipartite_cor_matrices <- purrr::imap(
  joint_cor_matrices,
  function(cor_matrix, site_name) {
    bact_id <- rownames(aligned_matrices[[site_name]]$matx_b)
    fungi_id <- rownames(aligned_matrices[[site_name]]$matx_f)

    sub_cor_matrix <- cor_matrix
    sub_cor_matrix[bact_id, bact_id] <- 0
    sub_cor_matrix[fungi_id, fungi_id] <- 0

    sub_cor_matrix
  }
)

# Verify: check that within-kingdom blocks are zeroed
purrr::imap(bipartite_cor_matrices, function(mat, site) {
  bact_id <- rownames(aligned_matrices[[site]]$matx_b)
  fungi_id <- rownames(aligned_matrices[[site]]$matx_f)
  list(
    site = site,
    length_bact_id = length(bact_id),
    length_fungi_id = length(fungi_id),
    bact_bact_nonzero = sum(mat[bact_id, bact_id] != 0),
    fungi_fungi_nonzero = sum(mat[fungi_id, fungi_id] != 0),
    cross_nonzero = sum(mat[bact_id, fungi_id] != 0)
  )
})


save(
  bipartite_cor_matrices,
  file = "data/output/networks/bipartite_cor_matrices.rda"
)

# Subnetworks ----

# amf_node_list = read.table("example_data/OTU_list_AMF-in-network_22.txt")$V1

# my_cor_matrix <- read.table("example_data/correlation_matrix_bipartite.txt")
# my_cutoff <- 0.757 # the correlation cutoff was determined in MENAP (http://ieg4.rccc.ou.edu/MENA/)

ef_bpt_net <- build_network(
  cor_mat = bipartite_cor_matrices[[1]], # ef site
  bact_ids = rownames(aligned_matrices$ef$matx_b),
  fungi_ids = rownames(aligned_matrices$ef$matx_f),
  cutoff = 0.25, # Testing with this cutoff. I need to determine with MENAP.
  kind = "cross"
) # list(graph, trans_mat) for the network containing all fungal-bacterial links


my_graph <- ef_bpt_net$graph

# tran_matrix_nonAMF = build_network(
#   cor_mat   = my_cor_matrix, bact_ids = ..., fungi_ids = ...,
#   cutoff    = my_cutoff, kind = "cross", rmv = amf_node_list
# )
# write.table(tran_matrix_nonAMF$trans_mat, "example_data/transition-matrix_nonAMF-bacteria.txt", sep = "\t", quote = F)

# tran_matrix_AMF = build_network(
#   cor_mat   = my_cor_matrix, bact_ids = ..., fungi_ids = ...,
#   cutoff    = my_cutoff, kind = "cross", keep = amf_node_list
# )
# write.table(tran_matrix_AMF$trans_mat, "example_data/transition-matrix_AMF-bacteria.txt", sep = "\t", quote = F)

write.table(
  ef_bpt_net$trans_mat,
  "data/output/networks/ef_joint_transition-matrix2.txt",
  sep = "\t",
  quote = F
)

# Analyze networks ----
# Calculates network topological features, node and link attributes, and
# generates input files for visualization using Cytoscape and Gephi.

########################

my_tranmatx <- ef_bpt_net$trans_mat

p_link <- sum(my_tranmatx > 0) / 2 # number of positive links
n_link <- sum(my_tranmatx < 0) / 2 # number of negative links
prop_p_link <- p_link / (p_link + n_link) # proportion of positive links

adj_mat <- my_tranmatx
adj_mat[adj_mat != 0] <- 1 # get adjacency matrix for the network. All links are indicated by 1 in adjacency matrix (for both + and - links).
write.table(
  adj_mat,
  "data/output/networks/ef_joint_adjacency-matrix.txt",
  sep = "\t"
)

# graph properties
conn_nodes <- gorder(my_graph) # node number
links <- ecount(my_graph) # link number
r2 <- fit_power_law(my_graph) # r2 of power law fit
avgK <- mean(centr_degree(my_graph)$res) # average degree
avgCC <- transitivity(my_graph, type = "average", isolates = "zero") # average clustering coefficient. Zero for bipartite network (no triangular subnetwork).
GD <- mean_distance(my_graph, directed = F, unconnected = T) # geodesic distance
gd <- cluster_fast_greedy(my_graph) # greedy module separation
modules <- length(gd) # number of modules
M <- modularity(gd) # modularity
largest_connected <- round(
  max(component_distribution(my_graph) * conn_nodes),
  0
) # node number in the largest connected component

# bipartite properties from the adjacency matrix
rows_16s <- rownames(aligned_matrices$ef$matx_b)[
  rownames(aligned_matrices$ef$matx_b) %in% rownames(adj_mat)
]
cols_its <- rownames(aligned_matrices$ef$matx_f)[
  rownames(aligned_matrices$ef$matx_f) %in% rownames(adj_mat)
]

bpt_matx <- adj_mat[rows_16s, cols_its]

bipartite_result = networklevel(
  index = c(
    "connectance",
    "nestedness",
    "NODF",
    "weighted connectance",
    "number of species",
    "cluster coefficient",
    "niche overlap",
    "partner diversity"
  ),
  bpt_matx
)

# add graph attributes — link_sign already set by build_network()
my_graph <- add_node_attribute(graph = my_graph)

#  output cytoscape and gephi input files for visualisation
cyto_gephi_output(graph_with_all_attributes = my_graph)

# Random bipartite network generation and analysis ----
# this code generates random bipartite networks that preserve the link and node numbers but rewire the links among nodes.
# then calculates the mean and standard deviation of network properties for multiple random networks.

# input: adjacency matrix
my_adj_mat <- read.table(
  "data/output/networks/ef_joint_adjacency-matrix.txt",
  sep = "\t",
  row.names = 1,
  header = T
)

# Identify bacteria and fungi node IDs present in the adjacency matrix
ef_bact_in_adj <- rownames(aligned_matrices$ef$matx_b)[
  rownames(aligned_matrices$ef$matx_b) %in% rownames(my_adj_mat)
]
ef_fungi_in_adj <- rownames(aligned_matrices$ef$matx_f)[
  rownames(aligned_matrices$ef$matx_f) %in% rownames(my_adj_mat)
]

# generate random networks
times <- 3 # number of random networks to generate, depending on the size of the empirical network. 100 times randomization would be good for ~several hundred nodes.
rand_adj_mat_list <- lapply(
  as.list(c(1:times)),
  FUN = rand_adj_gen,
  adj = my_adj_mat,
  bact_ids = ef_bact_in_adj,
  fungi_ids = ef_fungi_in_adj
)
names(rand_adj_mat_list) <- c(1:times)
length(rand_adj_mat_list) # check the number of random networks
dim(rand_adj_mat_list[[1]]) # check the dimension of the first random network

## Analyze random networks ----
random_global <- as.data.frame(lapply(rand_adj_mat_list, FUN = global))
names(random_global) <- c(1:times)
random_bpt <- as.data.frame(lapply(
  rand_adj_mat_list,
  FUN = bpt,
  bact_ids = ef_bact_in_adj,
  fungi_ids = ef_fungi_in_adj
))
names(random_bpt) <- c(1:times)
random_properties <- as.data.frame(rbind(random_global, random_bpt))

# get mean and sd
means <- rowMeans(random_properties)
sds <- apply(random_properties, 1, FUN = sd)


# Plotting ----
# Pull empirical values from the already-computed objects
empirical_global <- global(adj_mat)
empirical_bpt <- bpt(
  adj_mat,
  bact_ids = ef_bact_in_adj,
  fungi_ids = ef_fungi_in_adj
)
empirical_vals <- c(empirical_global, empirical_bpt)

# Build a tidy data frame; drop properties with zero SD (uninformative)
plot_df <- tibble(
  property = names(means),
  random_mean = means,
  random_sd = sds,
  empirical = empirical_vals[names(means)]
) |>
  filter(random_sd > 0) |>
  mutate(
    property = factor(property, levels = property),
    z_score = (empirical - random_mean) / random_sd
  )

# Panel A: empirical vs. random mean ± 2SD
p1 <- ggplot(plot_df, aes(x = property)) +
  geom_errorbar(
    aes(ymin = random_mean - 2 * random_sd, ymax = random_mean + 2 * random_sd),
    width = 0.3,
    colour = "grey60"
  ) +
  geom_point(aes(y = random_mean), colour = "grey40", size = 2) +
  geom_point(aes(y = empirical), colour = "#E05C5C", size = 3, shape = 17) +
  facet_wrap(~property, scales = "free", ncol = 4) +
  labs(
    title = "Empirical vs. random network properties",
    subtitle = "Grey: random mean ± 2SD  |  Red triangle: empirical",
    x = NULL,
    y = "Value"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 7)
  )

# Panel B: z-scores (how many SDs the empirical value departs from random)
p2 <- ggplot(plot_df, aes(x = property, y = z_score, fill = z_score > 0)) +
  geom_col() +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "grey40") +
  scale_fill_manual(
    values = c("TRUE" = "#4C8CBF", "FALSE" = "#E05C5C"),
    labels = c("TRUE" = "above random", "FALSE" = "below random"),
    name = NULL
  ) +
  labs(title = "Z-scores: empirical vs. random", x = NULL, y = "Z-score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

p1 / p2 # requires patchwork

ggsave(
  "data/output/networks/empirical_vs_random_network_properties.png",
  width = 10,
  height = 8,
  dpi = 300
)


## igraph ----

# Assign node type: bacteria vs fungi
V(my_graph)$kingdom <- ifelse(
  V(my_graph)$name %in% ef_bact_in_adj,
  "Bacteria",
  "Fungi"
)

V(my_graph)$type <- V(my_graph)$kingdom == "Fungi"

# Node aesthetics
node_colour <- ifelse(V(my_graph)$kingdom == "Bacteria", "#4C8CBF", "#E8A838")
node_size <- sqrt(igraph::degree(my_graph)) * 3 # scale by degree

# Edge colour by sign (positive/negative correlations)
edge_colour <- ifelse(
  E(my_graph)$link_sign == "positive",
  "#2ECC71AA",
  "#E05C5CAA"
)


# Layout: Fruchterman-Reingold works well for bipartite-ish networks
set.seed(42)
layout <- igraph::layout_with_fr(my_graph)

layout_bi <- igraph::layout_as_bipartite(my_graph)

ef_bipartite_network <- ggraph(my_graph, layout = layout) +
  geom_edge_arc(
    aes(colour = link_sign),
    width = 0.4,
    alpha = 0.7,
    circular = FALSE,
    strength = 0.1
  ) +
  geom_node_point(aes(color = "white", size = node_size * 1.03), shape = 21) +
  geom_node_point(aes(colour = kingdom, size = node_size)) +
  scale_edge_colour_manual(
    values = c("positive" = "#2ECC71AA", "negative" = "#E05C5CAA"),
    name = "Link sign"
  ) +
  scale_colour_manual(
    values = c("Bacteria" = "#4C8CBF", "Fungi" = "#E8A838"),
    name = "Kingdom"
  ) +
  scale_size_continuous(range = c(1, 8), guide = "none") +
  theme_graph() +
  labs(title = "EF bipartite network")


ef_bipartite_network
# plot(
#   my_graph,
#   layout = layout,
#   vertex.color = node_colour,
#   vertex.size = node_size,
#   vertex.label = NA, # hide long ASV names
#   vertex.frame.color = "white",
#   edge.color = edge_colour,
#   edge.width = 0.8,
#   edge.curved = 0.2,
#   margin = 0
# )

# legend(
#   "bottomleft",
#   legend = c("Bacteria", "Fungi", "Positive", "Negative"),
#   pch = c(21, 21, NA, NA),
#   lty = c(NA, NA, 1, 1),
#   col = c("white", "white", "#2ECC71", "#E05C5C"),
#   pt.bg = c("#4C8CBF", "#E8A838", NA, NA),
#   pt.cex = 2,
#   bty = "n",
#   text.col = "grey20"
# )

ggsave(
  filename = "data/output/networks/ef_bipartite_network.png",
  plot = ef_bipartite_network,
  width = 10,
  height = 8,
  dpi = 300
)

#TODO: clean  build_network() and centrality_shift() sections. Needs and data driven prev_filter() cutoff and correlation cutoff to be determined first.
# Insights

# EF cross-kingdom network analysis ----

# Arbitrary cutoff of 0.25 ----
# for exploratory analysis; will determine data-driven cutoffs later.

# Centrality shifts: bacteria-only → cross-kingdom for EF

ef_bact_id <- rownames(aligned_matrices$ef$matx_b)
ef_fungi_id <- rownames(aligned_matrices$ef$matx_f)
my_cutoff <- 0.25

ef_bact_g <- build_network(
  joint_cor_matrices$ef,
  ef_bact_id,
  ef_fungi_id,
  kind = "bacteria"
)
ef_fungi_g <- build_network(
  joint_cor_matrices$ef,
  ef_bact_id,
  ef_fungi_id,
  kind = "fungi"
)
ef_cross_g <- build_network(
  joint_cor_matrices$ef,
  ef_bact_id,
  ef_fungi_id,
  kind = "cross"
)

shifts <- centrality_shift(ef_bact_g$graph, ef_cross_g$graph, ef_bact_id)

# Summary of shifts
shift_summary1 <- shifts |>
  summarise(
    across(
      starts_with("delta"),
      list(mean = mean, sd = sd),
      .names = "{.col}_{.fn}"
    )
  )

shift_summary1

# ── Metrics table ─────────────────────────────────────────────────────────────

ef_metrics <- bind_rows(
  net_summary(ef_bact_g$graph, ef_bact_id, ef_fungi_id, "bacteria"),
  net_summary(ef_fungi_g$graph, ef_bact_id, ef_fungi_id, "fungi"),
  net_summary(ef_cross_g$graph, ef_bact_id, ef_fungi_id, "cross")
)
ef_metrics

# New cutoff for EF bacteria network ----
# EF bacteria network has 0 edges — cutoff is too high for this site.
# Check what cutoff gives a reasonable network for EF bacteria
cor_b <- joint_cor_matrices$ef[ef_bact_id, ef_bact_id]
quantile(abs(cor_b[upper.tri(cor_b)]), probs = c(0.90, 0.95, 0.97, 0.99, 0.999))


# my_cutoff = 0.25 is too lax for EF (only 468 bacteria, small dataset).
# Use the 99th percentile as a data-driven cutoff for EF specifically
ef_cutoff <- unname(quantile(abs(cor_b[upper.tri(cor_b)]), 0.99))

ef_bact_g2 <- build_network(
  joint_cor_matrices$ef,
  ef_bact_id,
  ef_fungi_id,
  cutoff = ef_cutoff,
  kind = "bacteria"
)
ef_fungi_g2 <- build_network(
  joint_cor_matrices$ef,
  ef_bact_id,
  ef_fungi_id,
  cutoff = ef_cutoff,
  kind = "fungi"
)
ef_cross_g2 <- build_network(
  joint_cor_matrices$ef,
  ef_bact_id,
  ef_fungi_id,
  cutoff = ef_cutoff,
  kind = "cross"
)

ef_metrics <- bind_rows(
  net_summary(ef_bact_g2$graph, ef_bact_id, ef_fungi_id, "bacteria"),
  net_summary(ef_fungi_g2$graph, ef_bact_id, ef_fungi_id, "fungi"),
  net_summary(ef_cross_g2$graph, ef_bact_id, ef_fungi_id, "cross")
)

cat("EF cutoff used:", round(ef_cutoff, 3), "\n")
ef_metrics

# Another cutoff ----

# cross network has 0 BF edges — the cross correlation matrix still has 0 bacteria-fungi entries
# because bipartite_cor_matrices zeroes those out. We need joint_cor_matrices (unzeroed) for cross.
# Check cross-kingdom correlations in joint vs bipartite
cor_cross <- joint_cor_matrices$ef[ef_bact_id, ef_fungi_id]
range(cor_cross)
quantile(abs(cor_cross), probs = c(0.90, 0.95, 0.99))


# Max B-F correlation is 0.39 — far below ef_cutoff (0.43).
# Use a lower cross-kingdom cutoff based on the 95th percentile of B-F correlations
ef_bf_cutoff <- unname(quantile(abs(cor_cross), 0.95))
cat("B-F cutoff:", round(ef_bf_cutoff, 3), "\n")

# Rebuild cross network with the B-F specific cutoff
ef_cross_g2 <- build_cross_network(
  joint_cor_matrices$ef,
  ef_bact_id,
  ef_fungi_id,
  within_cutoff = ef_cutoff,
  bf_cutoff = ef_bf_cutoff
)

ef_metrics <- bind_rows(
  net_summary(ef_bact_g2$graph, ef_bact_id, ef_fungi_id, "bacteria"),
  net_summary(ef_fungi_g2$graph, ef_bact_id, ef_fungi_id, "fungi"),
  net_summary(ef_cross_g2$graph, ef_bact_id, ef_fungi_id, "cross")
)
