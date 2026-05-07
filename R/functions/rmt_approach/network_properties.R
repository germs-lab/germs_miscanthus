library(igraph)
library(brainGraph)

#source("/Users/maggieyuan/Documents/!AMFnetwork/GitHub/fit_power_law.R") #fit power law degree distribution
# Previously called global
network_properties <- function(graph, power_law_engine) {
  size = gorder(graph)
  links = ecount(graph)
  r2 = fit_power_law_rmt(graph, engine = power_law_engine)
  avgK = mean(centr_degree(graph)$res)
  avgCC = transitivity(graph, type = "average", isolates = "zero")
  GD = mean_distance(graph, directed = F, unconnected = T)
  grdy = cluster_fast_greedy(graph)
  modules = length(grdy)
  M = modularity(grdy)
  largest_connected <- round(
    max(component_distribution(my_graph) * conn_nodes),
    0
  ) # node number in the largest connected component

  values = c(size, links, r2, avgK, avgCC, GD, modules, M)
  names(values) = c(
    "size",
    "links",
    "r2",
    "avgK",
    "avgCC",
    "GD",
    "modules",
    "M",
    "largest_connected"
  )
  return(values)
}
