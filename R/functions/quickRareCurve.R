# you will need to install the parallel package before hand
# and the vegan package if you dont have it
# https://dave-clark.github.io/post/speeding-up-rarefaction-curves-for-microbial-community-ecology/
quickRareCurve <- function(
  x,
  step = 1,
  sample,
  xlab = "Sample Size",
  ylab = "Species",
  label = TRUE,
  col,
  lty,
  max.cores = T,
  nCores = 1,
  ...
) {
  require(future)
  require(future.apply)

  plan(multisession) # Initiating  background R sessions on CURRENT machine

  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) {
    stop("function accepts only integers (counts)")
  }
  if (missing(col)) {
    col <- par("col")
  }
  if (missing(lty)) {
    lty <- par("lty")
  }
  library_size <- rowSums(x) # calculates library sizes
  species_num <- specnumber(x) # calculates n species for each sample
  if (any(species_num <= 0)) {
    message("empty rows removed")
    x <- x[species_num > 0, , drop = FALSE]
    library_size <- library_size[species_num > 0]
    species_num <- species_num[species_num > 0]
  } # removes any empty rows

  nr <- nrow(x) # number of samples
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)

  # Future_lapply
  # set number of cores
  mc <- ifelse(max.cores, parallelly::availableCores() - 1L, nCores)

  message(paste(
    "Using",
    mc,
    "cores.",
    "Max cores available:",
    parallelly::availableCores()
  ))

  out <- future_lapply(seq_len(nr), function(i) {
    n <- seq(1, library_size[i], by = step)
    if (n[length(n)] != library_size[i]) {
      n <- c(n, library_size[i])
    }
    drop(rarefy(x[i, ], n))
  })

  Nmax <- future_sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- future_sapply(out, max)
  p <- plot(
    c(1, max(Nmax)),
    c(1, max(Smax)),
    xlab = xlab,
    ylab = ylab,
    type = "n",
    ...
  )

  if (!missing(sample)) {
    abline(v = sample)
    rare <- future_sapply(out, function(z) {
      approx(x = attr(z, "Subsample"), y = z, xout = sample, rule = 1)$y
    })
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) {
    ordilabel(cbind(library_size, species_num), labels = rownames(x), ...)
  }

  return(list(
    rarefy_calcs = out,
    plots = p
  ))

  plan(sequential) # Explicit closing of R sessions

  message("Concurrent R sessions closed")
}
