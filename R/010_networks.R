################################################################################
# Cross-Kingdom Co-occurrence Networks for Miscanthus Microbiome
#
# Builds co-occurrence networks per site using SpiecEasi (MB or glasso method).
# Network types per site:
#   - Bacteria-only (16S DNA; all sequences)
#   - Fungi-only: ITS (LAMPS 2018) or AMF (EF, LAMPS 2022)
#   - Cross-kingdom: Bacteria + Fungi merged
#
# Outputs (per site):
#   - SpiecEasi fit objects
#   - igraph network objects
#   - Summary statistics table:
#       edge count, density, B-F edge fraction, pos:neg ratio,
#       degree/betweenness/eigenvector centrality (bacteria-only vs cross-kingdom),
#       connectance, transitivity, modularity, nestedness
#
# Author: Bolívar Aponte Rolón
# Date: 2026-04-24
################################################################################

source("R/utils/00_setup.R")
library("SpiecEasi") # v1.1.3+
library("igraph")
library("ggraph")

# ── 1. Load data ──────────────────────────────────────────────────────────────

# Assemble fungi physeqs per site (keyed to match main_16S_physeq_list names)
# EF         → AMF
# LAMPS 2018 → ITS (DNA only)
# LAMPS 2022 → AMF

main_mxg_physeq_list <- subset_to_miscanthus(
  main_physeq_list,
  crop_patterns = c("MXG", "M", "Miscanthus")
)


# ── 3. Pre-process ────────────────────────────────────────────────────────────

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
joint_filtered <- purrr::map2(
  bact_filtered,
  fungi_filtered,
  ~ merge_phyloseq(.x, .y)
)

# Per site: pair bacteria and fungi phyloseq objects without requiring shared
# sample names. SpiecEasi accepts two phyloseq objects directly for cross-kingdom
# inference; sample order must match, so we sort by sample name within each site.
aligned_matrices <- purrr::map(sites, function(s) {
  extract_matrices <- function(ps) {
    mat <- as.matrix(as.data.frame(otu_table(ps)))
    if (!taxa_are_rows(ps)) {
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

  list(ps1 = ps_b_mat, ps2 = ps_f_mat, ps3 = ps_j_mat)
}) |>
  set_names(sites)

# ── 4. SpiecEasi inference ────────────────────────────────────────────────────

cat("\nRunning SpiecEasi (this may take several minutes per site)...\n")

fits <- purrr::map(sites, function(s) {
  cat(" Site:", s, "\n")
  ncores <- 1
  list(
    mb = list(
      bact_only = run_spiec(
        aligned_matrices[[s]]$ps1,
        method = "mb",
        rep.num = 50,
        ncores = ncores
      ),
      fungi_only = run_spiec(
        aligned_matrices[[s]]$ps2,
        method = "mb",
        rep.num = 50,
        ncores = ncores
      ),
      joint_kingdom = run_spiec(
        aligned_matrices[[s]]$ps3,
        method = "mb",
        rep.num = 50,
        ncores = ncores
      )
    ),
    glasso = list(
      bact_only = run_spiec(
        aligned_matrices[[s]]$ps1,
        method = "glasso",
        rep.num = 50,
        ncores = ncores
      ),
      fungi_only = run_spiec(
        aligned_matrices[[s]]$ps2,
        method = "glasso",
        rep.num = 50,
        ncores = ncores
      ),
      joint_kingdom = run_spiec(
        aligned_matrices[[s]]$ps3,
        method = "glasso",
        rep.num = 50,
        ncores = ncores
      )
    )
  )
}) |>
  set_names(sites)


# ── 5. Build igraph objects ───────────────────────────────────────────────────

#

network_plots2 <- purrr::map(sites, function(s) {
  list(
    mb = list(
      bact_only = plot_igraph(
        fit = fits[[s]]$mb$bact_only, # remove $ef
        aligned_matrix1 = aligned_matrices[[s]]$ps1,
        adj_matrix = "StARS-refit",
        method = "mb",
        mode = "maxabs",
        label = NULL,
        kingdom = c("Bacteria")
      )
    )
  )
}) |>
  set_names(sites)

# ── 6. Summary statistics ─────────────────────────────────────────────────────

net_stats <- purrr::map_dfr(sites, function(s) {
  dplyr::bind_rows(
    net_summary(graphs[[s]]$bact_only, site = s, type = "bacteria_only"),
    net_summary(graphs[[s]]$fungi_only, site = s, type = "fungi_only"),
    net_summary(graphs[[s]]$cross, site = s, type = "cross_kingdom")
  )
})

cent_shifts <- purrr::map_dfr(sites, function(s) {
  centrality_shift(
    g_bact = graphs[[s]]$bact_only,
    g_cross = graphs[[s]]$cross,
    site = s
  )
})

net_stats

# ── 7. Save outputs ───────────────────────────────────────────────────────────

dir.create("data/output/networks", recursive = TRUE, showWarnings = FALSE)

save(fits, file = "data/output/networks/spiec_fits.rda")
save(graphs, file = "data/output/networks/igraph_objects.rda")

readr::write_csv(net_stats, "data/output/networks/network_summary_stats.csv")
readr::write_csv(cent_shifts, "data/output/networks/centrality_shifts.csv")

cat("\n010_networks.R complete.\n")
cat("Outputs saved to data/output/networks/\n")
