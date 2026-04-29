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

# ── 1. Load data ──────────────────────────────────────────────────────────────

load("data/output/rdata/main_16S_physeq_list_05.rda") # main_16S_physeq_list

load("data/output/rdata/phyloseq/ef_physeq_list.rda") # ef_physeq_list
load("data/output/rdata/phyloseq/lamps_2022_physeq_list.rda") # lamps_2022_physeq_list
load("data/output/rdata/phyloseq/lamps_2018_physeq_list.rda") # lamps_2018_physeq_list

# Assemble fungi physeqs per site (keyed to match main_16S_physeq_list names)
# EF         → AMF
# LAMPS 2018 → ITS (DNA only)
# LAMPS 2022 → AMF

main_mxg_physeq_list <- subset_to_miscanthus(
  main_physeq_list,
  crop_patterns = c("MXG", "M", "Miscanthus")
)
fungi_physeq_list <- list(
  ef_AMF = main_mxg_physeq_list$ef$ef_AMF_DNA,
  lamps_2018_16S_DNA = main_mxg_physeq_list$lamps_2018$lamps_2018_ITS_DNA,
  lamps_2022_16S_DNA = main_mxg_physeq_list$lamps_2022$lamps_2022_AMF_DNA
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
bact_filtered <- purrr::map(site_bact_list, prev_filter)
fungi_filtered <- purrr::map(site_fungi_list, prev_filter)

# Per site: pair bacteria and fungi phyloseq objects without requiring shared
# sample names. SpiecEasi accepts two phyloseq objects directly for cross-kingdom
# inference; sample order must match, so we sort by sample name within each site.
aligned <- purrr::map(sites, function(s) {
  ps_b <- bact_filtered[[s]]
  ps_f <- fungi_filtered[[s]]

  # Sort samples independently — SpiecEasi will use positional correspondence
  ps_b <- prune_samples(sort(sample_names(ps_b)), ps_b)
  ps_f <- prune_samples(sort(sample_names(ps_f)), ps_f)

  message(
    s,
    " – bacteria samples: ",
    nsamples(ps_b),
    " | fungi samples: ",
    nsamples(ps_f),
    " | bacteria ASVs: ",
    ntaxa(ps_b),
    " | fungi ASVs: ",
    ntaxa(ps_f)
  )

  list(ps1 = ps_b, ps2 = ps_f)
}) |>
  set_names(sites)

# ── 4. SpiecEasi inference ────────────────────────────────────────────────────

cat("\nRunning SpiecEasi (this may take several minutes per site)...\n")

fits <- purrr::map(sites, function(s) {
  cat(" Site:", s, "\n")
  list(
    bact_only = run_spiec(aligned[[s]]$ps1, ncores = 1),
    fungi_only = run_spiec(aligned[[s]]$ps2, ncores = 1),
    cross = run_spiec(aligned[[s]]$ps1, aligned[[s]]$ps2, ncores = 1)
  )
}) |>
  set_names(sites)

# ── 5. Build igraph objects ───────────────────────────────────────────────────

graphs <- purrr::map(sites, function(s) {
  list(
    bact_only = spiec_to_igraph(fits[[s]]$bact_only, aligned[[s]]$ps1),
    fungi_only = spiec_to_igraph(fits[[s]]$fungi_only, aligned[[s]]$ps2),
    cross = spiec_to_igraph(fits[[s]]$cross, aligned[[s]]$ps1, aligned[[s]]$ps2)
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
