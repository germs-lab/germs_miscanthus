## 011_network_visual.R
## Network visualizations for GERMS Miscanthus project
## Reference objects from: 010_networks_RMT.R

source("R/utils/00_setup.R")

# 0. Shared aesthetics & helpers ----

site_labels <- c(
  ef = "Energy Farm",
  lamps_2018 = "LAMPS 2018",
  lamps_2022 = "LAMPS 2022"
)

net_type_labels <- c(
  bact_net = "Bacteria",
  fungi_net = "Fungi",
  bxf_net = "Cross-kingdom",
  asym_net = "Asymmetric"
)

site_colors <- c(
  ef = "#4E9A8C",
  lamps_2018 = "#D4845A",
  lamps_2022 = "#7B6BB5"
)

kingdom_colors <- c(
  Bacteria = "#3A7AB5",
  Fungi = "#C45E5E"
)

theme_network <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

outdir <- "data/output/figures/networks"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width = 7, height = 5) {
  ggplot2::ggsave(
    filename = file.path(outdir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}


## 1. Network size and edge counts ----

net_props <- purrr::imap_dfr(
  all_networks,
  function(net, net_name) {
    props <- network_properties(net$graph, power_law_engine = "igraph")
    nms <- names(props)
    nms[is.null(nms) | is.na(nms) | nms == ""] <- paste0(
      "V",
      which(is.null(nms) | is.na(nms) | nms == "")
    )
    names(props) <- make.unique(nms)
    site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
    type <- stringr::str_remove(net_name, paste0("^", site, "_"))
    tibble::as_tibble_row(as.list(props)) |>
      dplyr::mutate(site = site, net_type = type, .before = 1)
  }
) |>
  dplyr::mutate(
    site_label = site_labels[site],
    type_label = net_type_labels[net_type]
  )
## 1a. Bar chart: nodes and edges per network
p_size <- net_props |>
  dplyr::filter(net_type %in% c("bact_net", "fungi_net", "bxf_net")) |>
  tidyr::pivot_longer(
    cols = c(node_size, links),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = dplyr::recode(metric, node_size = "Nodes", links = "Edges")
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = site_label, y = value, fill = site_label)) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::facet_grid(metric ~ type_label, scales = "free_y") +
  ggplot2::scale_fill_manual(values = site_colors) +
  ggplot2::labs(
    title = "Network size by site and type",
    x = NULL,
    y = "Count"
  ) +
  theme_network

save_plot(p_size, "01_network_size.pdf", width = 8, height = 5)

## 1b. Heatmap: topological properties (z-scored) across networks
prop_cols <- c(
  "avgK",
  "transitivity_avgCC",
  "GD",
  "modules",
  "grdy_modularity",
  "largest_connected"
)

p_heatmap <- net_props |>
  dplyr::filter(!is.na(avgK)) |>
  dplyr::select(site_label, type_label, dplyr::all_of(prop_cols)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(prop_cols),
    names_to = "property",
    values_to = "value"
  ) |>
  dplyr::group_by(property) |>
  dplyr::mutate(value_scaled = as.numeric(scale(value))) |>
  dplyr::ungroup() |>
  ggplot2::ggplot(ggplot2::aes(
    x = type_label,
    y = site_label,
    fill = value_scaled
  )) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::facet_wrap(~property, ncol = 3) +
  ggplot2::scale_fill_gradient2(
    low = "#3A7AB5",
    mid = "white",
    high = "#C45E5E",
    midpoint = 0,
    na.value = "grey85",
    name = "Z-score"
  ) +
  ggplot2::labs(
    title = "Topological properties across networks (z-scored)",
    x = NULL,
    y = NULL
  ) +
  theme_network +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

save_plot(
  p_heatmap,
  "02_network_properties_heatmap.pdf",
  width = 10,
  height = 7
)


## 2. RMT cutoffs ----

p_rmt <- rmt_cutoff_summary |>
  tidyr::pivot_longer(
    cols = c(bacteria_cutoff, fungi_cutoff, bxf_cutoff),
    names_to = "type",
    values_to = "cutoff"
  ) |>
  dplyr::mutate(
    type = dplyr::recode(
      type,
      bacteria_cutoff = "Bacteria",
      fungi_cutoff = "Fungi",
      bxf_cutoff = "Cross-kingdom"
    ),
    site_label = site_labels[site]
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = type, y = cutoff, fill = site_label)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::scale_fill_manual(values = site_colors, name = "Site") +
  ggplot2::labs(
    title = "RMT correlation cutoffs by site and network type",
    x = NULL,
    y = "Correlation cutoff"
  ) +
  ggplot2::ylim(0, 1) +
  theme_network

save_plot(p_rmt, "03_rmt_cutoffs.pdf", width = 6, height = 4)


## 3. Degree distributions ----

degree_df <- purrr::imap_dfr(
  all_networks,
  function(net, net_name) {
    g <- net$graph
    site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
    type <- stringr::str_remove(net_name, paste0("^", site, "_"))
    if (igraph::vcount(g) == 0) {
      return(NULL)
    }
    tibble::tibble(
      degree = igraph::degree(g),
      site = site,
      net_type = type
    )
  }
)

p_degree <- degree_df |>
  dplyr::filter(net_type %in% c("bact_net", "fungi_net", "bxf_net")) |>
  dplyr::mutate(
    site_label = site_labels[site],
    type_label = net_type_labels[net_type]
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = degree, fill = site_label)) +
  ggplot2::geom_histogram(bins = 30, alpha = 0.8, color = "white") +
  ggplot2::facet_grid(type_label ~ site_label, scales = "free") +
  ggplot2::scale_fill_manual(values = site_colors, guide = "none") +
  ggplot2::labs(
    title = "Degree distributions by network type and site",
    x = "Degree",
    y = "Count"
  ) +
  theme_network

save_plot(p_degree, "04_degree_distributions.pdf", width = 9, height = 7)


## 4. Within-kingdom network layouts ----
## Fruchterman-Reingold layout; node color = module, size = degree

plot_kingdom_network <- function(graph, title, seed = 42) {
  if (igraph::vcount(graph) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No network") +
        ggplot2::theme_void() +
        ggplot2::ggtitle(title)
    )
  }

  set.seed(seed)
  layout_fr <- igraph::layout_with_fr(graph)
  membership <- igraph::cluster_fast_greedy(graph)$membership

  edge_df <- igraph::as_edgelist(graph, names = FALSE) |>
    as.data.frame() |>
    dplyr::rename(from = V1, to = V2) |>
    dplyr::mutate(
      x = layout_fr[from, 1],
      y = layout_fr[from, 2],
      xend = layout_fr[to, 1],
      yend = layout_fr[to, 2]
    )

  node_df <- tibble::tibble(
    x = layout_fr[, 1],
    y = layout_fr[, 2],
    module = factor(membership),
    degree = igraph::degree(graph)
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey70",
      alpha = 0.4,
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = x, y = y, color = module, size = degree),
      alpha = 0.8
    ) +
    ggplot2::scale_size_continuous(range = c(1, 5), guide = "none") +
    ggplot2::guides(color = "none") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, hjust = 0.5))
}

kingdom_plots <- purrr::imap(
  all_networks[stringr::str_detect(names(all_networks), "bact_net|fungi_net")],
  function(net, net_name) {
    site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
    type <- stringr::str_remove(net_name, paste0("^", site, "_"))
    title <- paste0(site_labels[site], "\n", net_type_labels[type])
    plot_kingdom_network(net$graph, title = title)
  }
)

p_kingdom_grid <- patchwork::wrap_plots(kingdom_plots, ncol = 3) +
  patchwork::plot_annotation(
    title = "Within-kingdom co-occurrence networks",
    subtitle = "Node size = degree  |  Color = module",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 13))
  )

save_plot(p_kingdom_grid, "05_kingdom_networks.pdf", width = 12, height = 8)


## 5. Cross-kingdom network layout — LAMPS 2022 (ef + lamps_2018 are empty) ----
## Node color = kingdom (bacteria vs. fungi)

g_bxf <- all_networks[["lamps_2022_bxf_net"]]$graph
bact_ids <- rownames(aligned_matrices[["lamps_2022"]][["matx_b"]])
fungi_ids <- rownames(aligned_matrices[["lamps_2022"]][["matx_f"]])

if (igraph::vcount(g_bxf) > 0) {
  set.seed(42)
  layout_bxf <- igraph::layout_with_fr(g_bxf)
  node_names <- igraph::V(g_bxf)$name

  edge_bxf <- igraph::as_edgelist(g_bxf, names = FALSE) |>
    as.data.frame() |>
    dplyr::rename(from = V1, to = V2) |>
    dplyr::mutate(
      x = layout_bxf[from, 1],
      y = layout_bxf[from, 2],
      xend = layout_bxf[to, 1],
      yend = layout_bxf[to, 2]
    )

  node_bxf <- tibble::tibble(
    x = layout_bxf[, 1],
    y = layout_bxf[, 2],
    name = node_names,
    kingdom = dplyr::case_when(
      name %in% bact_ids ~ "Bacteria",
      name %in% fungi_ids ~ "Fungi",
      TRUE ~ "Unknown"
    ),
    degree = igraph::degree(g_bxf)
  )

  p_bxf <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_bxf,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey65",
      alpha = 0.35,
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      data = node_bxf,
      ggplot2::aes(x = x, y = y, color = kingdom, size = degree),
      alpha = 0.85
    ) +
    ggplot2::scale_color_manual(values = kingdom_colors, name = "Kingdom") +
    ggplot2::scale_size_continuous(range = c(1.5, 6), name = "Degree") +
    ggplot2::labs(
      title = "Cross-kingdom co-occurrence network — LAMPS 2022",
      subtitle = paste0(
        igraph::vcount(g_bxf),
        " nodes  |  ",
        igraph::ecount(g_bxf),
        " edges  |  ",
        sum(node_names %in% bact_ids),
        " bacteria  |  ",
        sum(node_names %in% fungi_ids),
        " fungi"
      )
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey40")
    )

  save_plot(p_bxf, "06_bxf_network_lamps2022.pdf", width = 9, height = 7)
}


## 6. Centrality shifts — LAMPS 2022 only ----
## Δ = cross-kingdom centrality − within-kingdom centrality

## 6a. Bacteria: Δ degree vs Δ eigenvector
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
    title = "Bacteria: centrality shift — within vs. cross-kingdom",
    subtitle = "Δ = cross − within  |  negative = reduced centrality in cross-kingdom",
    x = "Δ Degree",
    y = "Δ Eigenvector centrality"
  ) +
  theme_network

save_plot(p_bact_shift, "07_bact_centrality_shift.pdf", width = 9, height = 5)

## 6b. Fungi: Δ degree vs Δ eigenvector
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
    title = "Fungi: centrality shift — within vs. cross-kingdom",
    subtitle = "Δ = cross − within  |  negative = reduced centrality in cross-kingdom",
    x = "Δ Degree",
    y = "Δ Eigenvector centrality"
  ) +
  theme_network

save_plot(p_fungi_shift, "08_fungi_centrality_shift.pdf", width = 9, height = 5)

## 6c. Cross-kingdom degree histogram
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


## 7. Bipartite network metrics — LAMPS 2022 only ----

bpt_valid <- purrr::keep(bpt_properties, ~ !all(is.na(.x$bpt_result)))

if (length(bpt_valid) > 0) {
  bpt_df <- purrr::imap_dfr(
    bpt_valid,
    function(x, net_name) {
      site <- stringr::str_extract(net_name, "^ef|^lamps_2018|^lamps_2022")
      tibble::tibble(
        site = site,
        metric = names(x$bpt_result),
        value = as.numeric(x$bpt_result)
      )
    }
  ) |>
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
      title = "Bipartite network properties — cross-kingdom",
      x = NULL,
      y = NULL
    ) +
    theme_network

  save_plot(p_bpt, "10_bipartite_properties.pdf", width = 9, height = 5)
}
