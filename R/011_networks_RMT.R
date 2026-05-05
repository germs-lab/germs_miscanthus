# Exploring networks with

# "Fungal-Bacterial Cooccurrence Patterns Differ between Arbuscular Mycorrhizal Fungi and Nonmycorrhizal Fungi across Soil Niches". Citation: Yuan MM, Kakouridis A, Starr E, Nguyen N, Shi S, Pett-Ridge J, Nuccio E, Zhou J, Firestone M. 2021. Fungal-bacterial cooccurrence patterns differ between arbuscularmycorrhizal fungi and nonmycorrhizal fungi across soil niches. mBio 12:e03509-20. https://doi.org/10.1128/mBio.03509-20.

# This code calculates the full correlation matrix for constructing the network from OTU tables.

#  Load data ----

source("R/utils/00_setup.R")
library("igraph")
library("ggraph")
library(brainGraph)
library(bipartite)


# Assemble fungi physeqs per site (keyed to match main_16S_physeq_list names)
# EF         → AMF
# LAMPS 2018 → ITS (DNA only)
# LAMPS 2022 → AMF

main_mxg_physeq_list <- subset_to_miscanthus(
  main_physeq_list,
  crop_patterns = c("MXG", "M", "Miscanthus")
)


# Pre-process ----

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
write.table(
  joint_cor_matrices[[1]],
  file = "data/output/networks/correlation_matrix_2.txt",
  sep = "\t",
  quote = F
)


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

# this code is to get subnetworks by removing certain nodes or keeping certain nodes.

# input: cormatx - correlation matrix
# input: cutoff - correlation cutoff
# input: keep - list of nodes to keep
# input: rmv - list of nodes to remove

# output: transition matrix,
# where positive correlation above cutoff is 1, and negative correlation above cutoff is -1.
# All correlations weaker than cutoff are 0.
# OTUs with no link to other OTUs are removed from the matrix.

cor2tran <- function(cormatx, cutoff, keep, rmv) {
  tranmat <- cormatx
  diag(tranmat) <- 0
  tranmat[abs(tranmat) < cutoff] <- 0
  nod_id_tranmat <- which(rowSums(tranmat) != 0)
  tranmat_noiso <- tranmat[nod_id_tranmat, nod_id_tranmat]

  keep_id <- which(names(tranmat_noiso) %in% keep)
  rmv_id <- which(names(tranmat_noiso) %in% rmv)

  if (length(rmv) > 0 & length(keep) > 0) {
    print("Do not accept keep and remove at the same time")
  } else {
    if (length(keep) > 0 & length(keep_id) == 0) {
      print("Cannot find node to keep.")
    } else {
      if (length(keep_id) > 0) {
        nb_id = c()
        for (i in 1:length(keep_id)) {
          check_col = tranmat_noiso[, keep_id[i]]
          nb_id = c(nb_id, which(check_col != 0))
        }
        keep_id = unique(c(keep_id, nb_id))
        tranmat_sub = tranmat_noiso[keep_id, keep_id]
      }

      if (length(rmv_id) > 0) {
        tranmat_sub = tranmat_noiso[-rmv_id, -rmv_id]
      }

      if (length(rmv) > 0 & length(rmv_id) == 0) {
        tranmat_sub = tranmat_noiso
      }

      if (length(rmv) == 0 & length(keep) == 0) {
        tranmat_sub = tranmat_noiso
      }

      nod_id = which(rowSums(tranmat_sub) != 0)
      tranmat_sub_noiso = tranmat_sub[nod_id, nod_id]
    }
  }
  return(tranmat_sub_noiso)
}

amf_node_list = read.table("example_data/OTU_list_AMF-in-network_22.txt")$V1

my_cor_matrix <- read.table("example_data/correlation_matrix_bipartite.txt")
my_cutoff <- 0.757 # the correlation cutoff was determined in MENAP (http://ieg4.rccc.ou.edu/MENA/)

ef_tran_joint_matx <- cor2tran(
  cormatx = bipartite_cor_matrices[[1]], #ef site,
  cutoff = 0.25, # Testing with this cutoff. I need to determine with MENAP.
  rmv = c(),
  keep = c()
) # transition matrix for the network containing all fungal-bacterial links


write.table(
  ef_tran_joint_matx,
  "example_data/transition-matrix_fungi-bacteria.txt",
  sep = "\t",
  quote = F
)

# tran_matrix_nonAMF = cor2tran(
#   cormatx = my_cor_matrix,
#   cutoff = my_cutoff,
#   rmv = amf_node_list,
#   keep = c()
# ) # transition matrix for the network containing all nonAMF-bacterial links
# write.table(
#   tran_matrix_nonAMF,
#   "example_data/transition-matrix_nonAMF-bacteria.txt",
#   sep = "\t",
#   quote = F
# )

# tran_matrix_AMF = cor2tran(
#   cormatx = my_cor_matrix,
#   cutoff = my_cutoff,
#   rmv = c(),
#   keep = amf_node_list
# ) # transition matrix for the network containing all AMF-bacterial links
# write.table(
#   tran_matrix_AMF,
#   "example_data/transition-matrix_AMF-bacteria.txt",
#   sep = "\t",
#   quote = F
# )

# Analyze networks ----
# this code calculates network topological features, node and link attributes, and generates input files for visualization using Cytoscape and Gephi.

# Load
# library(brainGraph)
# library(bipartite)

## Helper functions ----
source("R/functions/rmt_approach/fit_power_law.R")

add_node_attribute = function(graph) {
  biomarker = gsub("OTU_.*_", "", vertex_attr(graph)$name)

  OTU_name = vertex_attr(graph)$name # OTU_name[1:10]
  node_degree = centr_degree(graph)$res

  module_separation = cluster_fast_greedy(graph)
  module_membership = membership(module_separation)
  pi = part_coeff(g = graph, memb = module_membership)
  zi = within_module_deg_z_score(g = graph, memb = module_membership)
  role = rep("peripherals", gorder(graph))
  role[which(pi >= 6.2 & zi < 2.5)] = "connector"
  role[which(pi >= 6.2 & zi >= 2.5)] = "network_hub"
  role[which(pi < 6.2 & zi >= 2.5)] = "module_hub"

  graph = set_vertex_attr(graph, "biomarker", index = V(graph), biomarker)
  graph = set_vertex_attr(graph, "node_degree", index = V(graph), node_degree)
  graph = set_vertex_attr(
    graph,
    "module_membership",
    index = V(graph),
    module_membership
  )
  graph = set_vertex_attr(graph, "pi", index = V(graph), pi)
  graph = set_vertex_attr(graph, "zi", index = V(graph), zi)
  graph = set_vertex_attr(graph, "vertex_role", index = V(graph), role)

  return(graph)
}

add_link_attribute = function(graph, tranmatx) {
  # get positive vs negative links from transition matrix
  mat_p_n = tranmatx
  mat_p_n[mat_p_n > 0] <- "positive"
  mat_p_n[mat_p_n < 0] <- "negative"

  # get edge attribute of negative/positive
  edge_names = unlist(strsplit(
    base::attr(E(graph), "vnames"),
    "|",
    fixed = T
  ))
  edge_rownames = edge_names[seq(1, (length(edge_names) - 1), by = 2)]
  edge_colnames = edge_names[seq(2, length(edge_names), by = 2)]
  edge_dir = unlist(lapply(
    c(1:length(edge_rownames)),
    FUN = function(id, np_matrix, rownames, colnames) {
      dir = np_matrix[
        which(rownames(np_matrix) == edge_rownames[id]),
        which(colnames(np_matrix) == edge_colnames[id])
      ]
    },
    np_matrix = mat_p_n,
    rownames = edge_rownames,
    colnames = edge_colnames
  ))

  graph = set_edge_attr(graph, "link_sign", index = E(graph), edge_dir)

  return(graph)
}

cyto_gephi_output = function(graph_with_all_attributes) {
  node_table = as.data.frame(vertex.attributes(graph_with_all_attributes))

  edge_table = as.data.frame(edge.attributes(graph_with_all_attributes))
  edge_names = unlist(strsplit(
    base::attr(E(graph_with_all_attributes), "vnames"),
    "|",
    fixed = T
  ))
  edge_table$node1 = edge_names[seq(
    from = 1,
    to = (length(edge_names) - 1),
    by = 2
  )]
  edge_table$node2 = edge_names[seq(from = 2, to = length(edge_names), by = 2)]

  gephi_node_table = node_table
  names(gephi_node_table)[1] <- "Id"

  gephi_edge_table = edge_table
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


########################

# input: transition matrix
my_tranmatx <- ef_tran_joint_matx

p_link = sum(my_tranmatx > 0) / 2 # number of positive links
n_link = sum(my_tranmatx < 0) / 2 # number of negative links
prop_p_link = p_link / (p_link + n_link) # proportion of positive links

adj_mat <- my_tranmatx
adj_mat[adj_mat != 0] <- 1 # get adjacency matrix for the network. All links are indicated by 1 in adjacency matrix (for both + and - links).
write.table(
  adj_mat,
  "data/output/networks/ef_joint_adjacency-matrix.txt",
  sep = "\t"
)

# construct graph
my_graph = graph_from_adjacency_matrix(
  as.matrix(adj_mat),
  mode = "undirected",
  weighted = NULL,
  diag = FALSE,
  add.colnames = NULL
)

# graph properties
conn_nodes = gorder(my_graph) # node number
links = ecount(my_graph) # link number
r2 = fit_power_law(my_graph) # r2 of power law fit
avgK = mean(centr_degree(my_graph)$res) # average degree
avgCC = transitivity(my_graph, type = "average", isolates = "zero") # average clustering coefficient. Zero for bipartite network (no triangular subnetwork).
GD = mean_distance(my_graph, directed = F, unconnected = T) # geodesic distance
gd = cluster_fast_greedy(my_graph) # greedy module separation
modules = length(gd) # number of modules
M = modularity(gd) # modularity
largest_connected = round(max(component_distribution(my_graph) * conn_nodes), 0) # node number in the largest connected component

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

# add graph attributes
my_graph = add_node_attribute(graph = my_graph)
my_graph = add_link_attribute(graph = my_graph, tranmatx = my_tranmatx)

#  output cytoscape and gephi input files for visulization
cyto_gephi_output(graph_with_all_attributes = my_graph)

# Random bipartite network generation and analysis ----
# this code generates random bipartite networks that preserve the link and node numbers but rewire the links among nodes.
# then calculates the mean and standard deviation of network properties for multiple random networks.

source("R/functions/rmt_approach/global.R")
source("R/functions/rmt_approach/bpt.R")
source("R/functions/rmt_approach/rand_adj_gen.R")


# input: adjacency matrix
my_adj_mat = read.table(
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
times = 3 # number of random networks to generate, depending on the size of the empirical network. 100 times randomization would be good for ~several hundred nodes.
rand_adj_mat_list = lapply(
  as.list(c(1:times)),
  FUN = rand_adj_gen,
  adj = my_adj_mat,
  bact_ids = ef_bact_in_adj,
  fungi_ids = ef_fungi_in_adj
)
names(rand_adj_mat_list) = c(1:times)
length(rand_adj_mat_list) # check the number of random networks
dim(rand_adj_mat_list[[1]]) # check the dimension of the first random network

## analyze random networks ----
random_global = as.data.frame(lapply(rand_adj_mat_list, FUN = global))
names(random_global) = c(1:times)
random_bpt = as.data.frame(lapply(
  rand_adj_mat_list,
  FUN = bpt,
  bact_ids = ef_bact_in_adj,
  fungi_ids = ef_fungi_in_adj
))
names(random_bpt) = c(1:times)
random_properties = as.data.frame(rbind(random_global, random_bpt))

# get mean and sd
means = rowMeans(random_properties)
sds = apply(random_properties, 1, FUN = sd)


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
library(igraph)

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
    values = c("positive" = "#2ECC71", "negative" = "#E05C5C"),
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
