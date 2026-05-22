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
source("R/functions/rmt_approach/network_properties.R")


# SECTION 1: Pre-process ----

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
# These should be informed by abundance occupancy curves and the distribution of prevalence values for each dataset. For now, I'm using 10% on all datasets as a starting point for exploratory analysis, but this should be data-driven and may differ for bacteria vs fungi and across sites. See "reports/mxg_abun_occu_report.pdf" for more details on the curves.

bact_filtered <- purrr::map(
  site_bact_list,
  ~ prev_filter(.x, prev_thresh = 0.1, min_count = 2)
)

fungi_filtered <- purrr::map(
  site_fungi_list,
  ~ prev_filter(.x, prev_thresh = 0.1, min_count = 2)
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

  # Sort samples independently - SpiecEasi will use positional correspondence
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
    " | bipartite ASVs: ",
    ntaxa(ps_j)
  )

  list(matx_b = ps_b_mat, matx_f = ps_f_mat, matx_j = ps_j_mat)
}) |>
  set_names(sites)

# SECTION 2: Correlation matrices ----
## Majority_corematrix ----
# Example with "ef" site.

# full_cor_matrices (505×505 for ef):     bxf_cor_matrices result:
# ┌─────────┬─────────┐                   ┌─────────┬─────────┐
# │  B×B    │  B×F    │        ->         │    0    │  B×F    │
# ├─────────┼─────────┤                   ├─────────┼─────────┤
# │  F×B    │  F×F    │                   │  F×B    │    0    │
# └─────────┴─────────┘                   └─────────┴─────────┘

full_cor_matrices <- purrr::imap(
  aligned_matrices,
  function(site_matx, site_name) {
    otu_matx <- aligned_matrices[[site_name]][["matx_j"]]
    otu_matx[is.na(otu_matx)] <- 0

    ## choose correlation methods from: 1) Spearman; 2) Pearson; 3) CLR-Pearson
    # 1) Spearman
    otu_matx_cor <- cor(t(otu_matx), method = "spearman")

    # # 2) Pearson, used in this study
    # otu_matx_cor = cor(t(otu_matx), method = "pearson")

    # 3) central log-ratio Pearson
    # clr <- apply((otu_matx + 1), 2, function(xc) {
    #   log(xc, 2) - sum(log(xc, 2)) / length(xc)
    # })
    # cor(t(clr), method = "pearson")

    # Remove ASVs involved in perfect off-diagonal correlations (|r| = 1).
    # These arise from proportionally identical abundance profiles (Spearman
    # ranks ties), not just exact duplicates. Iteratively drop one of each
    # perfectly-correlated pair until no off-diagonal |r| = 1 remains.
    diag(otu_matx_cor) <- NA
    repeat {
      perfect <- which(abs(otu_matx_cor) == 1, arr.ind = TRUE)
      if (nrow(perfect) == 0) {
        break
      }
      # Drop the taxon with the higher index in the first offending pair
      drop <- max(perfect[1, ])
      otu_matx_cor <- otu_matx_cor[-drop, -drop]
    }
    diag(otu_matx_cor) <- 1

    cli::cli_alert_info(
      "{site_name}: {nrow(otu_matx_cor)} ASVs retained after removing perfect correlations \\
      (started with {nrow(otu_matx)} ASVs)"
    )

    otu_matx_cor
  }
)

purrr::imap(full_cor_matrices, function(mat, site) {
  list(
    site = site,
    dim = dim(mat),
    abs_range = range(abs(mat), na.rm = TRUE)
  )
})


save(
  full_cor_matrices,
  file = "data/output/networks/full_correlation_matrices.rda"
)


## Cormatrix for bipartite network ----
# Stored as a rectangular B×F matrix (not the full square with zeroed blocks).
# This avoids inflating the eigenspectrum with structural zeros during RMT.

bxf_cor_matrices <- purrr::imap(
  full_cor_matrices,
  function(cor_matrix, site_name) {
    bact_id <- rownames(aligned_matrices[[site_name]]$matx_b)
    fungi_id <- rownames(aligned_matrices[[site_name]]$matx_f)

    # Keep only existing IDs after duplicate removal in full_cor_matrices
    bact_id <- intersect(bact_id, rownames(cor_matrix))
    fungi_id <- intersect(fungi_id, rownames(cor_matrix))

    # Rectangular B×F submatrix: bacteria rows × fungi columns
    cor_matrix[bact_id, fungi_id, drop = FALSE]
  }
)

# Verify: dimensions and range of cross-domain correlations
purrr::imap(bxf_cor_matrices, function(mat, site) {
  list(
    site = site,
    dim = dim(mat),
    abs_range = range(abs(mat), na.rm = TRUE)
  )
})


save(
  bxf_cor_matrices,
  file = "data/output/networks/bxf_correlation_matrices.rda"
)

# SECTION 3: Networks and Subnetworks ----

## SECTION 3.1: RMT cutoff for each site ----
ef_bact_rmt <- find_rmt_cutoff(
  full_cor_matrices$ef,
  bact_ids = rownames(aligned_matrices$ef$matx_b),
  kind = "bacteria",
  cutoff_seq = seq(1, 0.10, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)

ef_fungi_rmt <- find_rmt_cutoff(
  full_cor_matrices$ef,
  fungi_ids = rownames(aligned_matrices$ef$matx_f),
  kind = "fungi",
  cutoff_seq = seq(1, 0.10, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)


lamps2018_bact_rmt <- find_rmt_cutoff(
  full_cor_matrices$lamps_2018,
  bact_ids = rownames(aligned_matrices$lamps_2018$matx_b),
  kind = "bacteria",
  cutoff_seq = seq(1, 0.10, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)

lamps2018_fungi_rmt <- find_rmt_cutoff(
  full_cor_matrices$lamps_2018,
  fungi_ids = rownames(aligned_matrices$lamps_2018$matx_f),
  kind = "fungi",
  cutoff_seq = seq(1, 0.10, by = -0.01),
  n_bins = 15,
  alpha = 0.05,
  verbose = FALSE
)


lamps2022_bact_rmt <- find_rmt_cutoff(
  full_cor_matrices$lamps_2022,
  bact_ids = rownames(aligned_matrices$lamps_2022$matx_b),
  kind = "bacteria",
  cutoff_seq = seq(1, 0.10, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)

lamps2022_fungi_rmt <- find_rmt_cutoff(
  full_cor_matrices$lamps_2022,
  fungi_ids = rownames(aligned_matrices$lamps_2022$matx_f),
  kind = "fungi",
  cutoff_seq = seq(1, 0.10, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)

### BxF RMT ranges ----
# these are slightly different because the correlation values are different when we zero out the within-kingdom blocks. So we need to run RMT separately on the bxf matrices to get the appropriate cutoff for cross-kingdom edges.
# Check range first
bxf_rmt_ranges <- purrr::map(bxf_cor_matrices, ~ range(abs(.x), na.rm = TRUE))

ef_bxf_rmt <- find_rmt_cutoff(
  full_cor_matrices$ef,
  bact_ids = rownames(aligned_matrices$ef$matx_b),
  fungi_ids = rownames(aligned_matrices$ef$matx_f),
  kind = "cross",
  cutoff_seq = seq(max(bxf_rmt_ranges$ef), 0, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)

lamps2018_bxf_rmt <- find_rmt_cutoff(
  full_cor_matrices$lamps_2018,
  bact_ids = rownames(aligned_matrices$lamps_2018$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2018$matx_f),
  kind = "cross",
  cutoff_seq = seq(max(bxf_rmt_ranges$lamps_2018), 0, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)

lamps2022_bxf_rmt <- find_rmt_cutoff(
  full_cor_matrices$lamps_2022,
  bact_ids = rownames(aligned_matrices$lamps_2022$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2022$matx_f),
  kind = "cross",
  cutoff_seq = seq(max(bxf_rmt_ranges$lamps_2022), 0, by = -0.01),
  n_bins = 30,
  alpha = 0.05,
  verbose = FALSE
)

### Summary ----
rmt_cutoff_summary <- tibble(
  site = c("ef", "lamps_2018", "lamps_2022"),
  bacteria_cutoff = c(
    ef_bact_rmt$optimal_cutoff,
    lamps2018_bact_rmt$optimal_cutoff,
    lamps2022_bact_rmt$optimal_cutoff
  ),
  fungi_cutoff = c(
    ef_fungi_rmt$optimal_cutoff,
    lamps2018_fungi_rmt$optimal_cutoff,
    lamps2022_fungi_rmt$optimal_cutoff
  ),
  bxf_cutoff = c(
    ef_bxf_rmt$optimal_cutoff,
    lamps2018_bxf_rmt$optimal_cutoff,
    lamps2022_bxf_rmt$optimal_cutoff
  )
)

rmt_cutoff_summary

## SECTION 3.2: Within-kingdom networks ----

## EF site
ef_bact_net <- build_network(
  cor_mat = full_cor_matrices[[1]], # ef site
  bact_ids = rownames(aligned_matrices$ef$matx_b),
  fungi_ids = rownames(aligned_matrices$ef$matx_f),
  cutoff = ef_bact_rmt$optimal_cutoff, # data-driven cutoff from RMT approach
  kind = "bacteria"
)

ef_fungi_net <- build_network(
  cor_mat = full_cor_matrices[[1]], # ef site
  bact_ids = rownames(aligned_matrices$ef$matx_b),
  fungi_ids = rownames(aligned_matrices$ef$matx_f),
  cutoff = ef_fungi_rmt$optimal_cutoff,
  kind = "fungi"
)


## LAMPS 2018 site
lamps_2018_bact_net <- build_network(
  cor_mat = full_cor_matrices[[2]], # lamps_2018 site
  bact_ids = rownames(aligned_matrices$lamps_2018$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2018$matx_f),
  cutoff = lamps2018_bact_rmt$optimal_cutoff,
  kind = "bacteria"
)

lamps_2018_fungi_net <- build_network(
  cor_mat = full_cor_matrices[[2]], # lamps_2018 site
  bact_ids = rownames(aligned_matrices$lamps_2018$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2018$matx_f),
  cutoff = lamps2018_fungi_rmt$optimal_cutoff,
  kind = "fungi"
)


## LAMPS 2022 site
lamps_2022_bact_net <- build_network(
  cor_mat = full_cor_matrices[[3]], # lamps_2022 site
  bact_ids = rownames(aligned_matrices$lamps_2022$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2022$matx_f),
  cutoff = lamps2022_bact_rmt$optimal_cutoff,
  kind = "bacteria"
)

lamps_2022_fungi_net <- build_network(
  cor_mat = full_cor_matrices[[3]], # lamps_2022 site
  bact_ids = rownames(aligned_matrices$lamps_2022$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2022$matx_f),
  cutoff = lamps2022_fungi_rmt$optimal_cutoff,
  kind = "fungi"
)

## SECTION 3.3: Cross-kingdom networks ----

# EF site
ef_bxf_net <- build_network(
  cor_mat = bxf_cor_matrices[[1]], # ef site
  bact_ids = rownames(aligned_matrices$ef$matx_b),
  fungi_ids = rownames(aligned_matrices$ef$matx_f),
  cutoff = ef_bxf_rmt$optimal_cutoff,
  kind = "cross"
)

#  how cross-kingdom edges restructure a network that also has within-kingdom edges
ef_asym_net <- build_asym_network(
  cor_mat = full_cor_matrices[[1]], # ef site
  bact_ids = rownames(aligned_matrices$ef$matx_b),
  fungi_ids = rownames(aligned_matrices$ef$matx_f),
  within_bact_cutoff = ef_bact_rmt$optimal_cutoff,
  within_fungi_cutoff = ef_fungi_rmt$optimal_cutoff,
  bxf_cutoff = ef_bxf_rmt$optimal_cutoff
)

# LAMPS 2018
lamps_2018_bxf_net <- build_network(
  cor_mat = bxf_cor_matrices[[2]], # lamps_2018 site
  bact_ids = rownames(aligned_matrices$lamps_2018$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2018$matx_f),
  cutoff = lamps2018_bxf_rmt$optimal_cutoff,
  kind = "cross"
)

lamps_2018_asym_net <- build_asym_network(
  cor_mat = full_cor_matrices[[2]], # lamps_2018 site
  bact_ids = rownames(aligned_matrices$lamps_2018$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2018$matx_f),
  within_bact_cutoff = lamps2018_bact_rmt$optimal_cutoff,
  within_fungi_cutoff = lamps2018_fungi_rmt$optimal_cutoff,
  bxf_cutoff = lamps2018_bxf_rmt$optimal_cutoff
)

# LAMPS 2022

lamps_2022_bxf_net <- build_network(
  cor_mat = bxf_cor_matrices[[3]], # lamps_2022 site
  bact_ids = rownames(aligned_matrices$lamps_2022$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2022$matx_f),
  cutoff = lamps2022_bxf_rmt$optimal_cutoff,
  kind = "cross"
)

lamps_2022_asym_net <- build_asym_network(
  cor_mat = full_cor_matrices[[3]], # lamps_2022 site
  bact_ids = rownames(aligned_matrices$lamps_2022$matx_b),
  fungi_ids = rownames(aligned_matrices$lamps_2022$matx_f),
  within_bact_cutoff = lamps2022_bact_rmt$optimal_cutoff,
  within_fungi_cutoff = lamps2022_fungi_rmt$optimal_cutoff,
  bxf_cutoff = lamps2022_bxf_rmt$optimal_cutoff
)


### Save transition matrix ----
all_networks <- list(
  ef_bact_net = ef_bact_net,
  ef_fungi_net = ef_fungi_net,
  ef_bxf_net = ef_bxf_net,
  ef_asym_net = ef_asym_net,
  lamps_2018_bact_net = lamps_2018_bact_net,
  lamps_2018_fungi_net = lamps_2018_fungi_net,
  lamps_2018_bxf_net = lamps_2018_bxf_net,
  lamps_2018_asym_net = lamps_2018_asym_net,
  lamps_2022_bact_net = lamps_2022_bact_net,
  lamps_2022_fungi_net = lamps_2022_fungi_net,
  lamps_2022_bxf_net = lamps_2022_bxf_net,
  lamps_2022_asym_net = lamps_2022_asym_net
)


purrr::iwalk(
  all_networks,
  function(net, net_name) {
    write.table(
      net$trans_mat,
      paste0("data/output/networks/", net_name, "_transition-matrix.txt"),
      sep = "\t",
      quote = FALSE
    )
  }
)

### Save all networks ----
save(
  all_networks,
  file = "data/output/networks/all_networks.rda"
)


# SECTION 4: Analyze networks ----
# Calculates network topological features, node and link attributes, and
# generates input files for visualization using Cytoscape and Gephi.

# Graph properties
load("data/output/networks/all_networks.rda")

networks_properties_summary <- purrr::imap(
  all_networks,
  function(net, net_name) {
    site <- strsplit(net_name, "_")[[1]][1]
    kind <- strsplit(net_name, "_")[[1]][2]

    if (kind == "asym" | kind == "bxf") {
      kind <- "cross"
    }

    if (igraph::vcount(net$graph) == 0 || igraph::ecount(net$graph) == 0) {
      return(list(summary = NULL, properties = NULL))
    }

    list(
      summary = network_summary(net$graph, site = site, kind = kind),
      properties = network_properties(net$graph, power_law_engine = "OLS")
    )
  }
)

# centrality_shift = centrality_shift(
#   g_bact,
#   g_fungi,
#   g_cross,
#   site = site,
# )

## SECTION 4.1: Centrality shifts ----
# Compares node centrality (degree, betweenness, eigenvector) for bacteria and
# fungi in their within-kingdom networks vs. their position in the cross-kingdom
# (bxf) network.

centrality_shifts <- purrr::map(
  c("ef", "lamps_2018", "lamps_2022"),
  function(site) {
    g_bact <- all_networks[[paste0(site, "_bact_net")]]$graph
    g_fungi <- all_networks[[paste0(site, "_fungi_net")]]$graph
    g_cross <- all_networks[[paste0(site, "_bxf_net")]]$graph

    if (
      igraph::vcount(g_bact) == 0 ||
        igraph::vcount(g_fungi) == 0 ||
        igraph::vcount(g_cross) == 0
    ) {
      message(
        "Skipping centrality_shift for site '",
        site,
        "': one or more graphs are empty."
      )
      return(NULL)
    }

    centrality_shift(g_bact, g_fungi, g_cross, site = site)
  }
) |>
  set_names(c("ef", "lamps_2018", "lamps_2022"))

# Combined bacteria and fungi shift tables across sites
all_bact_shifts <- purrr::map(centrality_shifts, "bacteria_metrics") |>
  purrr::compact() |>
  dplyr::bind_rows()

all_fungi_shifts <- purrr::map(centrality_shifts, "fungi_metrics") |>
  purrr::compact() |>
  dplyr::bind_rows()

all_cross_centrality <- purrr::map(
  centrality_shifts,
  "cross_network_centrality"
) |>
  purrr::compact() |>
  dplyr::bind_rows()

## SECTION 4.2: Bipartite properties from the adjacency matrix ----

# Adjacency matrix makes binary networks for calculating bipartite network properties. The transition matrix has the correlation values, but we need to convert it to binary (0/1) for calculating bipartite properties.
adj_matrices <- purrr::imap(
  all_networks,
  function(net, net_name) {
    adj_mat <- net$trans_mat
    adj_mat[adj_mat != 0] <- 1 # get adjacency matrix for the network. All links are indicated by 1 in adjacency matrix (for both + and - links).
    list(
      adj_mat = adj_mat
    )
  }
)


purrr::iwalk(
  adj_matrices,
  function(net, net_name) {
    write.table(
      net$adj_mat,
      paste0("data/output/networks/", net_name, "_adjacency-matrix.txt"),
      sep = "\t",
      quote = FALSE
    )
  }
)


# Calculate bipartite network properties using the bipartite package.
# IDs of bacteria and fungi nodes in the adjacency matrix

# Inline: bpt
# Calculates bipartite network-level properties from an adjacency matrix.
bpt <- function(adj, bact_ids, fungi_ids) {
  adj <- as.matrix(adj)
  id_bact <- which(rownames(adj) %in% bact_ids)
  id_fungi <- which(rownames(adj) %in% fungi_ids)
  bpt_matx <- adj[id_bact, id_fungi, drop = FALSE]

  has_dims <- length(dim(bpt_matx)) > 0 &&
    nrow(bpt_matx) > 0 &&
    ncol(bpt_matx) > 0
  has_links <- has_dims && any(bpt_matx != 0)

  if (has_links) {
    bipartite_result <- bipartite::networklevel(
      bpt_matx,
      index = c(
        "connectance",
        "nestedness",
        "NODF",
        "weighted connectance",
        "number of species",
        "cluster coefficient",
        "niche overlap",
        "partner diversity"
      )
    )
  } else {
    bipartite_result <- rep(NA_real_, 13)
  }

  return(bipartite_result)
}

bpt_properties <- purrr::imap(
  adj_matrices[stringr::str_detect(names(adj_matrices), "_bxf_net")],
  function(net, net_name) {
    site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
    bact_ids <- rownames(aligned_matrices[[site]][["matx_b"]])
    fungi_ids <- rownames(aligned_matrices[[site]][["matx_f"]])

    list(
      adj_mat = net$adj_mat,
      bpt_result = bpt(net$adj_mat, bact_ids = bact_ids, fungi_ids = fungi_ids)
    )
  }
)

bpt_properties$ef_bxf_net$bpt_result
bpt_properties$lamps_2018_bxf_net$bpt_result


# Since no cross-kingdom edges passed the RMT cutoff for the EF and LAMPS 2018 site, the bxf network is empty and we cannot calculate bipartite properties for that site. Only LAMPS 2022 has a non-empty bxf network, so we can only calculate bipartite properties for that site.

bpt_properties$lamps_2022_bxf_net$bpt_result

# SECTION 4.3: Cytoscape and Gephi Visualization ----
#  Output cytoscape and gephi input files for visualisation

purrr::iwalk(
  all_networks,
  function(net, net_name) {
    cyto_gephi_output(
      graph_with_all_attributes = net$graph,
      file_path = "data/output/networks/cyto_gephi_inputs/",
      output_prefix = net_name
    )
  }
)


# SECTION 5: Random bipartite network generation and analysis ----
# Generates random bipartite networks that preserve node and link counts but
# rewire links among nodes. Only LAMPS 2022 has a valid cross-kingdom network
# (ef and lamps_2018 bxf networks are empty after RMT thresholding).
# 100 randomizations is recommended for ~several hundred nodes; 3 is used here
# for fast prototyping.

# Inline: rand_adj_gen
# Generates adjacency matrix for random networks based on the empirical
# adjacency matrix.
rand_adj_gen <- function(ID, adj, bact_ids, fungi_ids) {
  adj <- as.matrix(adj)
  id_bact <- which(rownames(adj) %in% bact_ids)
  id_fungi <- which(rownames(adj) %in% fungi_ids)
  network_L <- sum(adj[id_bact, id_fungi])

  rand_adj <- adj
  rand_adj[rand_adj != 0] <- 0

  if (length(id_bact) >= length(id_fungi)) {
    X <- id_bact
    Y <- id_fungi
  } else {
    X <- id_fungi
    Y <- id_bact
  }

  rand_order <- sample(c(
    Y,
    sample(Y, size = length(X) - length(Y), replace = TRUE)
  ))
  p <- 1
  for (j in X) {
    rand_adj[j, rand_order[p]] <- 1
    p <- p + 1
  }

  expanded <- as.numeric(rand_adj[X, Y])
  rand_remain_id <- sample(which(expanded == 0), size = (network_L - length(X)))
  expanded[rand_remain_id] <- 1
  rand_adj[X, Y] <- matrix(expanded, nrow = length(X), ncol = length(Y))
  rand_adj[Y, X] <- t(rand_adj[X, Y])

  if (sum(rand_adj) == sum(adj)) {
    return(rand_adj)
  } else {
    return("Error")
  }
}

# Use lamps_2022 bxf adjacency matrix — the only site with cross-kingdom edges
lamps2022_adj_mat <- adj_matrices[["lamps_2022_bxf_net"]]$adj_mat

lamps2022_bact_in_adj <- rownames(aligned_matrices$lamps_2022$matx_b)[
  rownames(aligned_matrices$lamps_2022$matx_b) %in% rownames(lamps2022_adj_mat)
]
lamps2022_fungi_in_adj <- rownames(aligned_matrices$lamps_2022$matx_f)[
  rownames(aligned_matrices$lamps_2022$matx_f) %in% rownames(lamps2022_adj_mat)
]

## SECTION 5.1: Generate random networks ----
times <- 3 # increase to ~100 for publication-quality null model comparison
rand_adj_mat_list <- lapply(
  as.list(seq_len(times)),
  FUN = rand_adj_gen,
  adj = lamps2022_adj_mat,
  bact_ids = lamps2022_bact_in_adj,
  fungi_ids = lamps2022_fungi_in_adj
)
names(rand_adj_mat_list) <- seq_len(times)

# Sanity checks
length(rand_adj_mat_list)
dim(rand_adj_mat_list[[1]])

## SECTION 5.2: Analyze random networks ----
# Network-level topological properties for each randomization
random_properties_list <- lapply(rand_adj_mat_list, function(rand_adj) {
  g_rand <- igraph::graph_from_adjacency_matrix(
    rand_adj,
    mode = "undirected",
    weighted = NULL,
    diag = FALSE
  )
  network_properties(g_rand, power_law_engine = "OLS")
})

random_properties <- dplyr::bind_rows(random_properties_list, .id = "replicate")

# Bipartite properties for each randomization
random_bpt_list <- lapply(rand_adj_mat_list, function(rand_adj) {
  result <- bpt(
    rand_adj,
    bact_ids = lamps2022_bact_in_adj,
    fungi_ids = lamps2022_fungi_in_adj
  )
  as.list(result)
})

random_bpt <- dplyr::bind_rows(random_bpt_list, .id = "replicate")

## SECTION 5.3: Null model summary (mean +/- sd) ----
# Compare empirical lamps_2022 bxf network against the null distribution
empirical_bpt <- bpt_properties$lamps_2022_bxf_net$bpt_result

null_summary <- tibble::tibble(
  metric = names(empirical_bpt),
  empirical = as.numeric(empirical_bpt),
  null_mean = colMeans(random_bpt[, -1], na.rm = TRUE),
  null_sd = apply(random_bpt[, -1], 2, sd, na.rm = TRUE)
) |>
  dplyr::mutate(
    z_score = (empirical - null_mean) / null_sd
  )

null_summary
