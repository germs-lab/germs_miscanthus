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
aligned <- purrr::map(sites, function(s) {
  ps_b <- bact_filtered[[s]]
  ps_f <- fungi_filtered[[s]]
  ps_j <- joint_filtered[[s]]

  # Sort samples independently — SpiecEasi will use positional correspondence
  ps_b <- prune_samples(sort(sample_names(ps_b)), ps_b)
  ps_f <- prune_samples(sort(sample_names(ps_f)), ps_f)
  ps_j <- prune_samples(sort(sample_names(ps_j)), ps_j)

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

  list(ps1 = ps_b, ps2 = ps_f, ps3 = ps_j)
}) |>
  set_names(sites)

# ── 4. SpiecEasi inference ────────────────────────────────────────────────────

cat("\nRunning SpiecEasi (this may take several minutes per site)...\n")

fits <- purrr::map(sites, function(s) {
  cat(" Site:", s, "\n")
  list(
    mb = list(
      bact_only = run_spiec(
        aligned[[s]]$ps1,
        method = "mb",
        rep.num = 50,
        ncores = 1
      ),
      fungi_only = run_spiec(
        aligned[[s]]$ps2,
        method = "mb",
        rep.num = 50,
        ncores = 1
      ),
      joint_kingdom = run_spiec(
        aligned[[s]]$ps3,
        method = "mb",
        rep.num = 50,
        ncores = 1
      )
    ),
    glasso = list(
      bact_only = run_spiec(
        aligned[[s]]$ps1,
        method = "glasso",
        rep.num = 50,
        ncores = 1
      ),
      fungi_only = run_spiec(
        aligned[[s]]$ps2,
        method = "glasso",
        rep.num = 50,
        ncores = 1
      ),
      joint_kingdom = run_spiec(
        aligned[[s]]$ps3,
        method = "glasso",
        rep.num = 50,
        ncores = 1
      )
    )
  )
}) |>
  set_names(sites)


# ── 5. Build igraph objects ───────────────────────────────────────────────────

## getRefit
ig.mb <- adj2igraph(getRefit(fits$ef$bact_only))
vsize <- rowMeans(clr(bact_filtered$ef@otu_table, 1)) + 6
am.coord <- layout_with_fr(ig.mb)

plot(
  ig.mb,
  layout = am.coord,
  vertex.size = vsize,
  vertex.label = NA,
  main = "MB"
)


beta_mat <- symBeta(getOptBeta(fits$ef$bact_only), mode = "maxabs")
diag(beta_mat) <- 0

g <- graph_from_adjacency_matrix(
  beta_mat,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

bact_ids <- taxa_names(aligned$ef$ps1)
if (!is.null(ps2)) {
  fung_ids <- taxa_names(ps2)
  V(g)$name <- c(bact_ids, fung_ids)
  V(g)$kingdom <- c(
    rep("Bacteria", length(bact_ids)),
    rep("Fungi", length(fung_ids))
  )
} else {
  V(g)$name <- bact_ids
  V(g)$kingdom <- rep("Bacteria", length(bact_ids))
}
E(g)$sign <- sign(as.numeric(E(g)$weight))
g

graphs <- purrr::map(sites, function(s) {
  list(
    bact_only = spiec_to_igraph(
      fits[[s]]$bact_only,
      aligned[[s]]$ps1
    ),
    fungi_only = spiec_to_igraph(
      fits[[s]]$fungi_only,
      aligned[[s]]$ps2
    )
  )
}) |>
  set_names(sites)

graphs <- purrr::map(sites, function(s) {
  list(
    mb = list(
      bact_only = spiec_to_igraph(fits[[s]]$mb$bact_only, aligned[[s]]$ps1),
      fungi_only = spiec_to_igraph(fits[[s]]$mb$fungi_only, aligned[[s]]$ps2),
      joint_kingdom = spiec_to_igraph(
        fits[[s]]$mb$joint_kingdom,
        aligned[[s]]$ps1,
        aligned_ps2 = aligned[[s]]$ps2
      )
    ),
    glasso = list(
      bact_only = spiec_to_igraph(fits[[s]]$glasso$bact_only, aligned[[s]]$ps1),
      fungi_only = spiec_to_igraph(
        fits[[s]]$glasso$fungi_only,
        aligned[[s]]$ps2
      ),
      joint_kingdom = spiec_to_igraph(
        fits[[s]]$glasso$joint_kingdom,
        aligned[[s]]$ps1,
        aligned_ps2 = aligned[[s]]$ps2
      )
    )
  )
}) |>
  set_names(sites)

graphs$ef$bact_only
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
