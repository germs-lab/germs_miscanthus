# this function generates adjacency matrix for random networks based on the adjacency matrix for the empirical network.
# bact_ids / fungi_ids: character vectors of row/col names belonging to each kingdom

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
	} # X is the longer edge of the rectangle

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
