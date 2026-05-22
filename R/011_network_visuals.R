##########################################################################

## Network visualizations for GERMS Miscanthus project
## Reference objects from: 010_networks_RMT.R

# By Bolívar Aponte Rolón
# Last updated: 2026-05-22
##########################################################################
source("R/utils/00_setup.R")
library(igraph)
library(ggraph)
library(patchwork)

source("R/functions/rmt_approach/network_properties.R")
source("R/functions/rmt_approach/fit_power_law.R")

load("data/output/networks/all_networks.rda")

# Shared aesthetics ----

outdir <- "data/output/figures/networks/"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

save_plot <- function(p, filename, width = 10, height = 7) {
  ggplot2::ggsave(
    file.path(outdir, filename),
    plot = p,
    width = width,
    height = height
  )
}

site_labels <- c(
  ef = "Energy Farm",
  lamps_2018 = "LAMPS 2018",
  lamps_2022 = "LAMPS 2022"
)

# Keyed on display labels so scale_color/fill_manual matches aes(color = site_label)
site_colors <- c(
  "Energy Farm" = "#E69F00",
  "LAMPS 2018" = "#56B4E9",
  "LAMPS 2022" = "#009E73"
)

kingdom_colors <- c(bacteria = "#4472C4", fungi = "#ED7D31")

net_type_labels <- c(
  bact = "Bacteria",
  fungi = "Fungi",
  bxf = "Cross-kingdom",
  asym = "Asymmetric"
)

theme_network <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right"
  )

# SECTION 1a: Network size (nodes & edges) ----

net_props <- purrr::imap_dfr(all_networks, function(net, net_name) {
  g <- net$graph
  site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
  type <- stringr::str_remove(net_name, paste0("^", site, "_"))
  type <- stringr::str_remove(type, "_net$")
  tibble::tibble(
    net_name = net_name,
    site = site,
    type = type,
    nodes = igraph::vcount(g),
    edges = igraph::ecount(g)
  )
}) |>
  dplyr::mutate(
    site_label = site_labels[site],
    type_label = net_type_labels[type]
  )

p_size <- net_props |>
  tidyr::pivot_longer(c(nodes, edges), names_to = "metric", values_to = "n") |>
  ggplot2::ggplot(ggplot2::aes(x = type_label, y = n, fill = site_label)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::scale_fill_manual(values = site_colors, name = "Site") +
  ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::labs(
    title = "Network size by type and site",
    x = NULL,
    y = "Count"
  ) +
  theme_network

save_plot(p_size, "01_network_size.pdf", width = 9, height = 5)

# SECTION 1b: Topological properties heatmap ----

prop_cols <- c(
  "degree_mean",
  "clustering_coef",
  "modularity",
  "avg_path_length",
  "connectance",
  "hub_score_mean"
)

networks_properties_summary <- purrr::imap(
  all_networks,
  function(net, net_name) {
    g <- net$graph
    site <- strsplit(net_name, "_")[[1]][1]
    kind <- strsplit(net_name, "_")[[1]][2]
    if (kind %in% c("asym", "bxf")) {
      kind <- "cross"
    }
    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
      return(list(summary = NULL, properties = NULL))
    }
    list(
      summary = network_summary(g, site = site, kind = kind),
      properties = network_properties(g, power_law_engine = "OLS")
    )
  }
)
networks_properties_summary$ef_asym_net$properties
props_df <- purrr::imap_dfr(
  networks_properties_summary,
  function(x, nm) {
    props <- x$properties # adjust if the key is different
    if (is.null(props)) {
      return(NULL)
    }
    nms <- names(props)
    nms[is.null(nms) | is.na(nms) | nms == ""] <- paste0(
      "V",
      which(is.null(nms) | is.na(nms) | nms == "")
    )
    names(props) <- make.unique(nms)
    tibble::as_tibble_row(as.list(props)) |>
      dplyr::mutate(net_name = nm)
  }
) |>
  dplyr::right_join(net_props, join_by(net_name))

props_scaled <- props_df |>
  dplyr::rename(
    degree_mean = avgK,
    clustering_coef = transitivity_avgCC,
    modularity = grdy_modularity,
    avg_path_length = GD
  ) |>
  dplyr::select(net_name, site_label, dplyr::any_of(prop_cols)) |>
  dplyr::mutate(dplyr::across(dplyr::any_of(prop_cols), scale)) |>
  tidyr::pivot_longer(
    dplyr::any_of(prop_cols),
    names_to = "metric",
    values_to = "z"
  )

p_heatmap <- props_scaled |>
  ggplot2::ggplot(ggplot2::aes(x = metric, y = net_name, fill = z)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#D6604D",
    midpoint = 0,
    name = "Z-score"
  ) +
  ggplot2::labs(
    title = "Network topological properties (z-scored)",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

save_plot(p_heatmap, "02_network_properties_heatmap.pdf", width = 9, height = 6)

# SECTION 2: RMT cutoffs ----

p_rmt <- rmt_cutoff_summary |>
  tidyr::pivot_longer(-site, names_to = "kingdom", values_to = "cutoff") |>
  dplyr::mutate(
    site_label = site_labels[site],
    kingdom = stringr::str_remove(kingdom, "_cutoff")
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = kingdom, y = cutoff, fill = site_label)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::scale_fill_manual(values = site_colors, name = "Site") +
  ggplot2::labs(
    title = "RMT correlation cutoffs by site and network type",
    x = "Network type",
    y = "Cutoff |r|"
  ) +
  theme_network

save_plot(p_rmt, "03_rmt_cutoffs.pdf", width = 7, height = 4)

# SECTION 3: Degree distributions ----

degree_df <- purrr::imap_dfr(all_networks, function(net, net_name) {
  g <- net$graph
  if (igraph::vcount(g) == 0) {
    return(NULL)
  }
  site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
  type <- stringr::str_remove(net_name, paste0("^", site, "_net$|^", site, "_"))
  type <- stringr::str_remove(type, "_net$")
  tibble::tibble(
    site = site,
    type = type,
    degree = igraph::degree(g)
  )
}) |>
  dplyr::mutate(
    site_label = site_labels[site],
    type_label = net_type_labels[type]
  )

p_degree <- degree_df |>
  ggplot2::ggplot(ggplot2::aes(x = degree, fill = site_label)) +
  ggplot2::geom_histogram(bins = 25, color = "white", alpha = 0.85) +
  ggplot2::scale_fill_manual(values = site_colors, name = "Site") +
  ggplot2::facet_grid(type_label ~ site_label, scales = "free") +
  ggplot2::labs(
    title = "Degree distributions by network type and site",
    x = "Degree",
    y = "Count"
  ) +
  theme_network +
  ggplot2::theme(legend.position = "none")

save_plot(p_degree, "04_degree_distributions.pdf", width = 10, height = 9)

# SECTION 4: Within-kingdom network graphs ----

plot_kingdom_network <- function(net_name) {
  net <- all_networks[[net_name]]
  g <- net$graph
  if (igraph::vcount(g) == 0) {
    return(NULL)
  }
  site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
  kind <- dplyr::if_else(grepl("bact", net_name), "bacteria", "fungi")
  label <- paste(site_labels[site], kind)

  set.seed(42)
  modules <- igraph::cluster_louvain(igraph::as_undirected(g))
  node_df <- tibble::tibble(
    degree = igraph::degree(g),
    module = factor(igraph::membership(modules))
  )

  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.15, color = "grey70") +
    ggraph::geom_node_point(
      ggplot2::aes(size = node_df$degree, color = node_df$module),
      alpha = 0.8
    ) +
    ggplot2::scale_size_continuous(range = c(1, 5), name = "Degree") +
    ggplot2::scale_color_discrete(name = "Module") +
    ggplot2::labs(title = label) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}

kingdom_nets <- names(all_networks)[grepl(
  "_bact_net|_fungi_net",
  names(all_networks)
)]
kingdom_plots <- purrr::map(kingdom_nets, plot_kingdom_network) |>
  purrr::compact()

p_kingdom_grid <- patchwork::wrap_plots(kingdom_plots, ncol = 3) +
  patchwork::plot_annotation(
    title = "Within-kingdom co-occurrence networks",
    subtitle = "Layout = Fruchterman-Reingold | Node color = Louvain module | Node size = degree"
  )

save_plot(p_kingdom_grid, "05_kingdom_networks.pdf", width = 12, height = 8)

# SECTION 5: Cross-kingdom network (LAMPS 2022 only) ----
# ef and lamps_2018 bxf networks are empty after RMT thresholding.

g_bxf <- all_networks[["lamps_2022_bxf_net"]]$graph

if (igraph::vcount(g_bxf) > 0) {
  bact_ids_l22 <- rownames(aligned_matrices[["lamps_2022"]][["matx_b"]])
  fungi_ids_l22 <- rownames(aligned_matrices[["lamps_2022"]][["matx_f"]])

  set.seed(42)
  layout_bxf <- igraph::layout_with_fr(g_bxf)
  node_names <- igraph::V(g_bxf)$name

  node_bxf <- tibble::tibble(
    name = node_names,
    x = layout_bxf[, 1],
    y = layout_bxf[, 2],
    kingdom = dplyr::case_when(
      name %in% bact_ids_l22 ~ "bacteria",
      name %in% fungi_ids_l22 ~ "fungi",
      TRUE ~ "unknown"
    ),
    degree = igraph::degree(g_bxf)
  )

  edge_bxf <- igraph::as_edgelist(g_bxf, names = FALSE) |>
    as.data.frame() |>
    dplyr::rename(from = V1, to = V2) |>
    dplyr::mutate(
      x = layout_bxf[from, 1],
      y = layout_bxf[from, 2],
      xend = layout_bxf[to, 1],
      yend = layout_bxf[to, 2]
    )

  p_bxf <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_bxf,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      alpha = 0.12,
      color = "grey60",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      data = node_bxf,
      ggplot2::aes(x = x, y = y, color = kingdom, size = degree),
      alpha = 0.8
    ) +
    ggplot2::scale_color_manual(
      values = kingdom_colors,
      name = "Kingdom",
      labels = c(bacteria = "Bacteria", fungi = "Fungi")
    ) +
    ggplot2::scale_size_continuous(range = c(1, 6), name = "Degree") +
    ggplot2::labs(
      title = "Cross-kingdom co-occurrence network -- LAMPS 2022",
      subtitle = paste0(
        igraph::vcount(g_bxf),
        " nodes  |  ",
        igraph::ecount(g_bxf),
        " edges"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40"),
      legend.position = "right"
    )

  save_plot(p_bxf, "06_bxf_network_lamps2022.pdf", width = 9, height = 7)
}

# SECTION 6a: Bacteria centrality shifts ----
load("data/output/networks/centrality_shifts.rda")
p_bact_shift <- all_bact_shifts |>
  dplyr::mutate(site_label = site_labels[site]) |>
  ggplot2::ggplot(ggplot2::aes(x = delta_degree, y = delta_eigenvector)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_point(
    ggplot2::aes(color = site_label, size = deg_bact_only),
    alpha = 0.7
  ) +
  ggplot2::scale_color_manual(values = site_colors, name = "Site") +
  ggplot2::scale_size_continuous(
    range = c(1, 5),
    name = "Within-kingdom degree"
  ) +
  ggplot2::facet_wrap(~site_label) +
  ggplot2::labs(
    title = "Bacteria: centrality shift -- within vs. cross-kingdom",
    subtitle = "Delta = cross - within  |  negative = reduced centrality in cross-kingdom network",
    x = "Delta Degree",
    y = "Delta Eigenvector centrality"
  ) +
  theme_network

save_plot(p_bact_shift, "07_bact_centrality_shift.pdf", width = 9, height = 5)

# SECTION 6b: Fungi centrality shifts ----

p_fungi_shift <- all_fungi_shifts |>
  dplyr::mutate(site_label = site_labels[site]) |>
  ggplot2::ggplot(ggplot2::aes(x = delta_degree, y = delta_eigenvector)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  ggplot2::geom_point(
    ggplot2::aes(color = site_label, size = deg_fungi_only),
    alpha = 0.7
  ) +
  ggplot2::scale_color_manual(values = site_colors, name = "Site") +
  ggplot2::scale_size_continuous(
    range = c(1, 5),
    name = "Within-kingdom degree"
  ) +
  ggplot2::facet_wrap(~site_label) +
  ggplot2::labs(
    title = "Fungi: centrality shift -- within vs. cross-kingdom",
    subtitle = "Delta = cross - within  |  negative = reduced centrality in cross-kingdom network",
    x = "Delta Degree",
    y = "Delta Eigenvector centrality"
  ) +
  theme_network

save_plot(p_fungi_shift, "08_fungi_centrality_shift.pdf", width = 9, height = 5)

# SECTION 6c: Cross-kingdom degree distribution ----

p_cross_deg <- all_cross_centrality |>
  dplyr::mutate(site_label = site_labels[site]) |>
  ggplot2::ggplot(ggplot2::aes(x = deg_cross)) +
  ggplot2::geom_histogram(
    ggplot2::aes(fill = site_label),
    bins = 20,
    color = "white",
    alpha = 0.85
  ) +
  ggplot2::scale_fill_manual(values = site_colors, guide = "none") +
  ggplot2::facet_wrap(~site_label) +
  ggplot2::labs(
    title = "Degree distribution in cross-kingdom network",
    x = "Degree",
    y = "Count"
  ) +
  theme_network

save_plot(p_cross_deg, "09_cross_degree_dist.pdf", width = 8, height = 4)

# SECTION 7: Bipartite properties (valid sites only) ----
load("data/output/networks/bipartite_properties.rda")
bpt_valid <- purrr::keep(bpt_properties, ~ !all(is.na(.x$bpt_result)))

if (length(bpt_valid) > 0) {
  bpt_df <- purrr::imap_dfr(bpt_valid, function(x, net_name) {
    site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
    tibble::tibble(
      site = site,
      metric = names(x$bpt_result),
      value = as.numeric(x$bpt_result)
    )
  }) |>
    dplyr::mutate(site_label = site_labels[site])

  p_bpt <- bpt_df |>
    dplyr::filter(
      metric %in%
        c(
          "connectance",
          "nestedness",
          "NODF",
          "weighted connectance",
          "cluster coefficient"
        )
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x = site_label,
      y = value,
      fill = site_label
    )) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = site_colors) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = "Bipartite network properties -- cross-kingdom",
      x = NULL,
      y = NULL
    ) +
    theme_network

  save_plot(p_bpt, "10_bipartite_properties.pdf", width = 9, height = 5)
}
