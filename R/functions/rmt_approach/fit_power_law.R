# https://github.com/Mengting-Maggie-Yuan/Fungal-bacterial-network-construction-and-analysis/blob/main/fit_power_law.R
fit_power_law = function(graph) {
    # calculate degree
    d = centr_degree(graph)$res
    dd = degree_distribution(graph)
    degree = 1:max(d)
    probability = dd[-1]
    # delete blank values
    nonzero.position = which(probability != 0)
    # probability = probability[nonzero.position]
    # degree = degree[nonzero.position]
    degree <- igraph::degree(graph)
    degree <- degree[degree > 0] # remove zero-degree nodes

    if (length(unique(degree)) < 2) {
        return(NA) # can't fit a power law
    }

    probability <- degree / sum(degree)
    lm(log(probability) ~ log(degree))

    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    # print(paste("Alpha =", round(alpha, 3)))
    return(round(R.square, 3))
    # plot
    # plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", col = 1, main = "Degree Distribution")
    # curve(power.law.fit, col = "red", add = T, n = length(d))
}
