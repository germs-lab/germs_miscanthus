#source("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/fit_power_law.R") #fit power law degree distribution
# Previously called global
network_properties <- function(graph, power_law_engine) {
  # Guard against empty graphs (no vertices or no edges)
  if (vcount(graph) == 0 || ecount(graph) == 0) {
    cli::cli_alert_warning(
      "Graph is empty (no vertices or no edges). Returning NA values for all properties."
    )
    values <- rep(NA_real_, 9)
    names(values) <- c(
      "node_size",
      "links",
      "r2",
      "avgK",
      "transitivity_avgCC",
      "GD",
      "modules",
      "grdy_modularity",
      "largest_connected"
    )
    return(values)
  }

  node_size <- vcount(graph) # == gorder(graph)
  links <- ecount(graph)
  r2 <- fit_power_law_rmt(graph, engine = power_law_engine)
  avgK <- mean(centr_degree(graph)$res)
  transitivity_avgCC <- transitivity(graph, type = "average", isolates = "zero") # average clustering coefficient. Zero for bipartite network (no triangular subnetwork).
  GD <- mean_distance(graph, directed = FALSE, unconnected = TRUE)
  grdy <- cluster_fast_greedy(graph)
  modules <- length(grdy)
  grdy_modularity <- modularity(grdy)
  largest_connected <- round(
    max(component_distribution(graph) * node_size),
    0
  ) # node number in the largest connected component

  values <- c(
    node_size,
    links,
    r2,
    avgK,
    transitivity_avgCC,
    GD,
    modules,
    grdy_modularity,
    largest_connected
  )
  names(values) <- c(
    "node_size",
    "links",
    "r2",
    "avgK",
    "transitivity_avgCC",
    "GD",
    "modules",
    "grdy_modularity",
    "largest_connected"
  )
  return(values)
}
