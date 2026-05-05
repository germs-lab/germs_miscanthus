library(igraph)
library(brainGraph)
library(bipartite)

bpt <- function(adj, bact_ids, fungi_ids) {
  adj <- as.matrix(adj)
  id_bact <- which(rownames(adj) %in% bact_ids)
  id_fungi <- which(rownames(adj) %in% fungi_ids)
  bpt_matx <- adj[id_bact, id_fungi]

  if (length(dim(bpt_matx)) > 0 && nrow(bpt_matx) > 0 && ncol(bpt_matx) > 0) {
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
  } else {
    bipartite_result <- rep(NA_real_, 13)
  }

  return(bipartite_result)
}
