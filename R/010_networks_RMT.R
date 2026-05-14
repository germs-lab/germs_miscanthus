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
source("R/functions/rmt_approach/bpt.R")
source("R/functions/rmt_approach/rand_adj_gen.R")


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

## SECTION 4.1: Bipartite properties from the adjacency matrix ----

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

# IDs of bacteria and fungi nodes in the adjacency matrix
rows_16s <- rownames(aligned_matrices$ef$matx_b)[
  rownames(aligned_matrices$ef$matx_b) %in% rownames(adj_mat)
]
cols_its <- rownames(aligned_matrices$ef$matx_f)[
  rownames(aligned_matrices$ef$matx_f) %in% rownames(adj_mat)
]

bpt_matx <- adj_mat[rows_16s, cols_its]

# Calculate bipartite network properties using the bipartite package.
bipartite_result <- networklevel(
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

# SECTION 4.2: Cytoscape and Gephi Visualization ----
#  Output cytoscape and gephi input files for visualisation
cyto_gephi_output(graph_with_all_attributes = my_graph)

# SECTION 5: Random bipartite network generation and analysis ----
# this code generates random bipartite networks that preserve the link and node numbers but rewire the links among nodes.
# then calculates the mean and standard deviation of network properties for multiple random networks.

# input: adjacency matrix from SECTION 4.1, which is the binary version of the transition matrix (i.e., all non-zero entries are set to 1).
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

## SECTION 5.1: Generate random networks ----
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

## SECTION 5.2: Analyze random networks ----
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

# ############################################
# # Plotting ----
# # Pull empirical values from the already-computed objects
# empirical_global <- global(adj_mat)
# empirical_bpt <- bpt(
#   adj_mat,
#   bact_ids = ef_bact_in_adj,
#   fungi_ids = ef_fungi_in_adj
# )
# empirical_vals <- c(empirical_global, empirical_bpt)

# # Build a tidy data frame; drop properties with zero SD (uninformative)
# plot_df <- tibble(
#   property = names(means),
#   random_mean = means,
#   random_sd = sds,
#   empirical = empirical_vals[names(means)]
# ) |>
#   filter(random_sd > 0) |>
#   mutate(
#     property = factor(property, levels = property),
#     z_score = (empirical - random_mean) / random_sd
#   )

# # Panel A: empirical vs. random mean ± 2SD
# p1 <- ggplot(plot_df, aes(x = property)) +
#   geom_errorbar(
#     aes(ymin = random_mean - 2 * random_sd, ymax = random_mean + 2 * random_sd),
#     width = 0.3,
#     colour = "grey60"
#   ) +
#   geom_point(aes(y = random_mean), colour = "grey40", size = 2) +
#   geom_point(aes(y = empirical), colour = "#E05C5C", size = 3, shape = 17) +
#   facet_wrap(~property, scales = "free", ncol = 4) +
#   labs(
#     title = "Empirical vs. random network properties",
#     subtitle = "Grey: random mean ± 2SD  |  Red triangle: empirical",
#     x = NULL,
#     y = "Value"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     strip.text = element_text(size = 7)
#   )

# # Panel B: z-scores (how many SDs the empirical value departs from random)
# p2 <- ggplot(plot_df, aes(x = property, y = z_score, fill = z_score > 0)) +
#   geom_col() +
#   geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "grey40") +
#   scale_fill_manual(
#     values = c("TRUE" = "#4C8CBF", "FALSE" = "#E05C5C"),
#     labels = c("TRUE" = "above random", "FALSE" = "below random"),
#     name = NULL
#   ) +
#   labs(title = "Z-scores: empirical vs. random", x = NULL, y = "Z-score") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# p1 / p2 # requires patchwork

# ggsave(
#   "data/output/networks/empirical_vs_random_network_properties.png",
#   width = 10,
#   height = 8,
#   dpi = 300
# )

# ## igraph ----

# # Assign node type: bacteria vs fungi
# V(my_graph)$kingdom <- ifelse(
#   V(my_graph)$name %in% ef_bact_in_adj,
#   "Bacteria",
#   "Fungi"
# )

# V(my_graph)$type <- V(my_graph)$kingdom == "Fungi"

# # Node aesthetics
# node_colour <- ifelse(V(my_graph)$kingdom == "Bacteria", "#4C8CBF", "#E8A838")
# node_size <- sqrt(igraph::degree(my_graph)) * 3 # scale by degree

# # Edge colour by sign (positive/negative correlations)
# edge_colour <- ifelse(
#   E(my_graph)$link_sign == "positive",
#   "#2ECC71AA",
#   "#E05C5CAA"
# )

# # Layout: Fruchterman-Reingold works well for bipartite-ish networks
# set.seed(42)
# layout <- igraph::layout_with_fr(my_graph)

# layout_bi <- igraph::layout_as_bipartite(my_graph)

# ef_bipartite_network <- ggraph(my_graph, layout = layout) +
#   geom_edge_arc(
#     aes(colour = link_sign),
#     width = 0.4,
#     alpha = 0.7,
#     circular = FALSE,
#     strength = 0.1
#   ) +
#   geom_node_point(aes(color = "white", size = node_size * 1.03), shape = 21) +
#   geom_node_point(aes(colour = kingdom, size = node_size)) +
#   scale_edge_colour_manual(
#     values = c("positive" = "#2ECC71AA", "negative" = "#E05C5CAA"),
#     name = "Link sign"
#   ) +
#   scale_colour_manual(
#     values = c("Bacteria" = "#4C8CBF", "Fungi" = "#E8A838"),
#     name = "Kingdom"
#   ) +
#   scale_size_continuous(range = c(1, 8), guide = "none") +
#   theme_graph() +
#   labs(title = "EF bipartite network")

# ef_bipartite_network
# # plot(
# #   my_graph,
# #   layout = layout,
# #   vertex.color = node_colour,
# #   vertex.size = node_size,
# #   vertex.label = NA, # hide long ASV names
# #   vertex.frame.color = "white",
# #   edge.color = edge_colour,
# #   edge.width = 0.8,
# #   edge.curved = 0.2,
# #   margin = 0
# # )

# # legend(
# #   "bottomleft",
# #   legend = c("Bacteria", "Fungi", "Positive", "Negative"),
# #   pch = c(21, 21, NA, NA),
# #   lty = c(NA, NA, 1, 1),
# #   col = c("white", "white", "#2ECC71", "#E05C5C"),
# #   pt.bg = c("#4C8CBF", "#E8A838", NA, NA),
# #   pt.cex = 2,
# #   bty = "n",
# #   text.col = "grey20"
# # )

# ggsave(
#   filename = "data/output/networks/ef_bipartite_network.png",
#   plot = ef_bipartite_network,
#   width = 10,
#   height = 8,
#   dpi = 300
# )

# #TODO: clean  build_network() and centrality_shift() sections. Needs and data driven prev_filter() cutoff and correlation cutoff to be determined first.
# # Insights

# # EF cross-kingdom network analysis ----

# # Arbitrary cutoff of 0.25 ----
# # for exploratory analysis; will determine data-driven cutoffs later.

# # Centrality shifts: bacteria-only → cross-kingdom for EF

# ef_bact_id <- rownames(aligned_matrices$ef$matx_b)
# ef_fungi_id <- rownames(aligned_matrices$ef$matx_f)
# my_cutoff <- 0.25

# ef_bact_g <- build_network(
#   joint_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   kind = "bacteria"
# )
# ef_fungi_g <- build_network(
#   joint_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   kind = "fungi"
# )
# ef_cross_g <- build_network(
#   joint_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   kind = "cross"
# )

# shifts <- centrality_shift(ef_bact_g$graph, ef_cross_g$graph, ef_bact_id)

# # Summary of shifts
# shift_summary1 <- shifts |>
#   summarise(
#     across(
#       starts_with("delta"),
#       list(mean = mean, sd = sd),
#       .names = "{.col}_{.fn}"
#     )
#   )

# shift_summary1

# # ── Metrics table ─────────────────────────────────────────────────────────────

# ef_metrics <- bind_rows(
#   net_summary(ef_bact_g$graph, ef_bact_id, ef_fungi_id, "bacteria"),
#   net_summary(ef_fungi_g$graph, ef_bact_id, ef_fungi_id, "fungi"),
#   net_summary(ef_cross_g$graph, ef_bact_id, ef_fungi_id, "cross")
# )
# ef_metrics

# # New cutoff for EF bacteria network ----
# # EF bacteria network has 0 edges - cutoff is too high for this site.
# # Check what cutoff gives a reasonable network for EF bacteria
# cor_b <- joint_cor_matrices$ef[ef_bact_id, ef_bact_id]
# quantile(abs(cor_b[upper.tri(cor_b)]), probs = c(0.90, 0.95, 0.97, 0.99, 0.999))

# # my_cutoff = 0.25 is too lax for EF (only 468 bacteria, small dataset).
# # Use the 99th percentile as a data-driven cutoff for EF specifically
# ef_cutoff <- unname(quantile(abs(cor_b[upper.tri(cor_b)]), 0.99))

# ef_bact_g2 <- build_network(
#   joint_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   cutoff = ef_cutoff,
#   kind = "bacteria"
# )
# ef_fungi_g2 <- build_network(
#   joint_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   cutoff = ef_cutoff,
#   kind = "fungi"
# )
# ef_cross_g2 <- build_network(
#   joint_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   cutoff = ef_cutoff,
#   kind = "cross"
# )

# ef_metrics <- bind_rows(
#   net_summary(ef_bact_g2$graph, ef_bact_id, ef_fungi_id, "bacteria"),
#   net_summary(ef_fungi_g2$graph, ef_bact_id, ef_fungi_id, "fungi"),
#   net_summary(ef_cross_g2$graph, ef_bact_id, ef_fungi_id, "cross")
# )

# cat("EF cutoff used:", round(ef_cutoff, 3), "\n")
# ef_metrics

# # Another cutoff ----

# # cross network has 0 BF edges - the cross correlation matrix still has 0 bacteria-fungi entries
# # because bipartite_cor_matrices zeroes those out. We need joint_cor_matrices (unzeroed) for cross.
# # Check cross-kingdom correlations in joint vs bipartite
# cor_cross <- joint_cor_matrices$ef[ef_bact_id, ef_fungi_id]
# range(cor_cross)
# quantile(abs(cor_cross), probs = c(0.90, 0.95, 0.99))

# # Max B-F correlation is 0.39 - far below ef_cutoff (0.43).
# # Use a lower cross-kingdom cutoff based on the 95th percentile of B-F correlations
# ef_bf_cutoff <- unname(quantile(abs(cor_cross), 0.95))
# cat("B-F cutoff:", round(ef_bf_cutoff, 3), "\n")

# # Rebuild cross network with the B-F specific cutoff
# ef_cross_g2 <- build_asym_network(
#   joint_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   within_cutoff = ef_cutoff,
#   bf_cutoff = ef_bf_cutoff
# )

# ef_metrics <- bind_rows(
#   net_summary(ef_bact_g2$graph, ef_bact_id, ef_fungi_id, "bacteria"),
#   net_summary(ef_fungi_g2$graph, ef_bact_id, ef_fungi_id, "fungi"),
#   net_summary(ef_cross_g2$graph, ef_bact_id, ef_fungi_id, "cross")
# )

# # Centrality shifts: bacteria-only → cross-kingdom ----
# # centrality_shift(g_bact, g_cross, site) compares degree/betweenness/eigenvector
# # centrality for bacteria in a bacteria-only graph vs. the bacterial subgraph of
# # the cross-kingdom network.

# # ── EF ──────────────────────────────────────────────────────────────────────

# ef_bact_id <- rownames(aligned_matrices$ef$matx_b)
# ef_fungi_id <- rownames(aligned_matrices$ef$matx_f)
# my_cutoff <- 0.25

# ef_bact_g <- build_network(
#   full_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   kind = "bacteria"
# )
# ef_fungi_g <- build_network(
#   full_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   kind = "fungi"
# )
# ef_bxf_g <- build_network(
#   bxf_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   kind = "cross"
# )
# ef_asym_g <- build_asym_network(
#   full_cor_matrices$ef,
#   ef_bact_id,
#   ef_fungi_id,
#   within_cutoff = my_cutoff,
#   bxf_cutoff = my_cutoff * 0.75
# )

# ef_shifts <- centrality_shift(ef_bact_g$graph, ef_bxf_g$graph, site = "ef")

# ef_metrics <- bind_rows(
#   network_summary(ef_bact_g$graph, site = "ef", kind = "bacteria"),
#   network_summary(ef_fungi_g$graph, site = "ef", kind = "fungi"),
#   network_summary(ef_bxf_g$graph, site = "ef", kind = "bxf"),
#   network_summary(ef_asym_g$graph, site = "ef", kind = "asym")
# )

# # ── LAMPS 2018 ───────────────────────────────────────────────────────────────

# lamps18_bact_id <- rownames(aligned_matrices$lamps_2018$matx_b)
# lamps18_fungi_id <- rownames(aligned_matrices$lamps_2018$matx_f)

# lamps18_bact_g <- build_network(
#   full_cor_matrices$lamps_2018,
#   lamps18_bact_id,
#   lamps18_fungi_id,
#   kind = "bacteria"
# )
# lamps18_fungi_g <- build_network(
#   full_cor_matrices$lamps_2018,
#   lamps18_bact_id,
#   lamps18_fungi_id,
#   kind = "fungi"
# )
# lamps18_bxf_g <- build_network(
#   bxf_cor_matrices$lamps_2018,
#   lamps18_bact_id,
#   lamps18_fungi_id,
#   kind = "cross"
# )
# lamps18_asym_g <- build_asym_network(
#   full_cor_matrices$lamps_2018,
#   lamps18_bact_id,
#   lamps18_fungi_id,
#   within_cutoff = my_cutoff,
#   bxf_cutoff = my_cutoff * 0.75
# )

# lamps18_shifts <- centrality_shift(
#   lamps18_bact_g$graph,
#   lamps18_bxf_g$graph,
#   site = "lamps_2018"
# )

# lamps18_metrics <- bind_rows(
#   network_summary(lamps18_bact_g$graph, site = "lamps_2018", kind = "bacteria"),
#   network_summary(lamps18_fungi_g$graph, site = "lamps_2018", kind = "fungi"),
#   network_summary(lamps18_bxf_g$graph, site = "lamps_2018", kind = "bxf"),
#   network_summary(lamps18_asym_g$graph, site = "lamps_2018", kind = "asym")
# )

# # ── LAMPS 2022 ───────────────────────────────────────────────────────────────

# lamps22_bact_id <- rownames(aligned_matrices$lamps_2022$matx_b)
# lamps22_fungi_id <- rownames(aligned_matrices$lamps_2022$matx_f)

# lamps22_bact_g <- build_network(
#   full_cor_matrices$lamps_2022,
#   lamps22_bact_id,
#   lamps22_fungi_id,
#   kind = "bacteria"
# )
# lamps22_fungi_g <- build_network(
#   full_cor_matrices$lamps_2022,
#   lamps22_bact_id,
#   lamps22_fungi_id,
#   kind = "fungi"
# )
# lamps22_bxf_g <- build_network(
#   bxf_cor_matrices$lamps_2022,
#   lamps22_bact_id,
#   lamps22_fungi_id,
#   kind = "cross"
# )
# lamps22_asym_g <- build_asym_network(
#   full_cor_matrices$lamps_2022,
#   lamps22_bact_id,
#   lamps22_fungi_id,
#   within_cutoff = my_cutoff,
#   bxf_cutoff = my_cutoff * 0.75
# )

# lamps22_shifts <- centrality_shift(
#   lamps22_bact_g$graph,
#   lamps22_bxf_g$graph,
#   site = "lamps_2022"
# )

# lamps22_metrics <- bind_rows(
#   network_summary(lamps22_bact_g$graph, site = "lamps_2022", kind = "bacteria"),
#   network_summary(lamps22_fungi_g$graph, site = "lamps_2022", kind = "fungi"),
#   network_summary(lamps22_bxf_g$graph, site = "lamps_2022", kind = "bxf"),
#   network_summary(lamps22_asym_g$graph, site = "lamps_2022", kind = "asym")
# )

# # ── Combined tables ───────────────────────────────────────────────────────────

# all_shifts <- bind_rows(ef_shifts, lamps18_shifts, lamps22_shifts)
# all_metrics <- bind_rows(ef_metrics, lamps18_metrics, lamps22_metrics)

# all_shifts
# all_metrics
