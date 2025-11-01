p_iNEXT <- function(
  x,
  q = c(0, 1, 2),
  datatype = "abundance",
  endpoint = NULL,
  knots = 40,
  conf = 0.95,
  nboot = 100,
  xlab = "Sample Size",
  ylab = "Diversity",
  label = TRUE,
  col,
  lty,
  max.cores = TRUE,
  nCores = 1,
  ...
) {
  require(future)
  require(future.apply)
  require(iNEXT)
  require(ggplot2)

  plan(multisession) # Initiating background R sessions on CURRENT machine
  x = spider
  x <- do.call(rbind, x)
  # x <- as.matrix(x)
  # if (!identical(all.equal(x, round(x)), TRUE)) {
  #   stop("function accepts only integers (counts)")
  # }

  library_size <- rowSums(x) # calculates library sizes
  species_num <- vegan::specnumber(x) # calculates n species for each sample
  if (any(species_num <= 0)) {
    message("empty rows removed")
    x <- x[species_num > 0, , drop = FALSE]
    library_size <- library_size[species_num > 0]
    species_num <- species_num[species_num > 0]
  } # removes any empty rows

  nr <- nrow(x) # number of samples
  if (missing(col)) {
    col <- par("col")
  }
  if (missing(lty)) {
    lty <- par("lty")
  }
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)

  # Set number of cores
  mc <- ifelse(max.cores, parallelly::availableCores() - 1L, nCores)

  message(paste(
    "Using",
    mc,
    "cores.",
    "Max cores available:",
    parallelly::availableCores()
  ))

  # Set endpoint for rarefaction/extrapolation
  if (is.null(endpoint)) {
    endpoint <- max(library_size) * 2 # Default: 2x max library size
  }

  # Parallel iNEXT computation
  out <- future_lapply(
    seq_len(nr),
    function(i) {
      iNEXT::iNEXT(
        x = x[i, ],
        q = q,
        datatype = datatype,
        endpoint = endpoint,
        knots = knots,
        conf = conf,
        nboot = nboot
      )
    },
    future.seed = TRUE
  )

  names(out) <- rownames(x)

  # Create combined data frame for ggplot
  inext_data <- purrr::map_dfr(
    out,
    function(result) {
      result$iNextEst %>%
        as.data.frame() %>%
        rownames_to_column(var = "index")
    },
    .id = "sample"
  )

  # # Generate ggplot visualization
  # p <- ggplot(
  #   inext_data,
  #   aes(
  #     x =coverage_based.m ,
  #     y = coverage_based.qD,
  #     color = sample,
  #     linetype = factor(coverage_based.Method)
  #   )
  # ) +
  #   geom_line(linewidth = 0.8) +
  #   # facet_wrap(
  #   #   ~q,
  #   #   labeller = labeller(
  #   #     q = c(
  #   #       "0" = "Richness (q=0)",
  #   #       "1" = "Shannon (q=1)",
  #   #       "2" = "Simpson (q=2)"
  #   #     )
  #   #   )
  #   # ) +
  #   labs(
  #     title = "iNEXT: Rarefaction/Extrapolation Curves",
  #     x = xlab,
  #     y = ylab,
  #     color = "Sample",
  #     linetype = "Method"
  #   ) +
  #   theme_minimal() +
  #   theme(
  #     legend.position = "bottom",
  #     strip.text = element_text(face = "bold")
  #   )

  # # Optional: geom_ribbon for confidence intervals
  # if (!is.na(conf)) {
  #   p <- p +
  #     geom_ribbon(
  #       aes(ymin = qD.LCL, ymax = qD.UCL, fill = sample),
  #       alpha = 0.1,
  #       color = NA
  #     )
  # }

  plan(sequential) # Explicit closing of R sessions
  message("Concurrent R sessions closed")

  return(list(
    inext_results = out,
    inext_data = combined_inext,
    plot = p
  ))
}


##########################
ggplotColors <- function(g) {
  d <- 360 / g # Calculate the distance between colors in HCL color space
  h <- cumsum(c(15, rep(d, g - 1))) # Create cumulative sums to define hue values
  hcl(h = h, c = 100, l = 65) # Convert HCL values to hexadecimal color codes
}

gg_inext_custom <- function(
  inext_data,
  type = 1,
  se = TRUE,
  facet.var = "None",
  color.var = "Assemblage",
  grey = FALSE
) {
  require(ggplot2)

  # Select columns based on type
  if (type == 1) {
    # Size-based rarefaction/extrapolation
    data <- inext_data %>%
      select(
        sample,
        size_based.m,
        size_based.qD,
        size_based.qD.LCL,
        size_based.qD.UCL,
        size_based.Method,
        size_based.Order.q
      ) %>%
      rename(
        m = size_based.m,
        qD = size_based.qD,
        qD.LCL = size_based.qD.LCL,
        qD.UCL = size_based.qD.UCL,
        Method = size_based.Method,
        Order.q = size_based.Order.q
      )
    x_label <- "Number of individuals"
    y_label <- "Species diversity"
  } else if (type == 2) {
    # Sample completeness
    data <- inext_data %>%
      select(
        sample,
        size_based.m,
        size_based.SC,
        size_based.SC.LCL,
        size_based.SC.UCL,
        size_based.Method,
        size_based.Order.q
      ) %>%
      rename(
        m = size_based.m,
        SC = size_based.SC,
        SC.LCL = size_based.SC.LCL,
        SC.UCL = size_based.SC.UCL,
        Method = size_based.Method,
        Order.q = size_based.Order.q
      )
    x_label <- "Number of individuals"
    y_label <- "Sample coverage"
  } else if (type == 3) {
    # Coverage-based rarefaction/extrapolation
    data <- inext_data %>%
      select(
        sample,
        coverage_based.SC,
        coverage_based.qD,
        coverage_based.qD.LCL,
        coverage_based.qD.UCL,
        coverage_based.Method,
        coverage_based.Order.q
      ) %>%
      rename(
        SC = coverage_based.SC,
        qD = coverage_based.qD,
        qD.LCL = coverage_based.qD.LCL,
        qD.UCL = coverage_based.qD.UCL,
        Method = coverage_based.Method,
        Order.q = coverage_based.Order.q
      )
    x_label <- "Sample coverage"
    y_label <- "Species diversity"
  }

  observed_endpoint <- data %>%
    filter(Method == "Observed") %>%
    mutate(
      Order.q = as.factor(Order.q),
      Method = factor(Method, levels = c("Rarefaction", "Extrapolation"))
    )

  data <- data %>%
    filter(Method %in% c("Rarefaction", "Extrapolation")) %>%
    mutate(
      Order.q = as.factor(Order.q),
      Method = factor(Method, levels = c("Rarefaction", "Extrapolation"))
    )

  # Colors
  # Check if the number of unique 'sample' is 8 or less
  if (length(unique(inext_data$sample)) <= 8) {
    cbPalette <- rev(c(
      "#999999",
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#330066",
      "#CC79A7",
      "#0072B2",
      "#D55E00"
    ))
  } else {
    # If there are more than 8 assemblages, start with the same predefined color palette
    # Then extend the palette by generating additional colors using the 'ggplotColors' function
    cbPalette <- rev(c(
      "#999999",
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#330066",
      "#CC79A7",
      "#0072B2",
      "#D55E00"
    ))
    cbPalette <- c(
      cbPalette,
      ggplotColors(length(unique(inext_data$sample)) - 8)
    )
  }

  # Determine color variable
  if (facet.var == "Order.q") {
    color_col <- "sample"
  }
  if (facet.var == "Assemblage") {
    color_col <- "Order.q"
  }
  if (facet.var == "Both") {
    data$col <- paste(data$sample, data$Order.q, sep = "-")
    color_col <- "col"
  }
  if (color.var == "Assemblage") {
    color_col <- "sample"
  }

  if (color.var == "Order.q") {
    color_col <- "Order.q"
  }

  if (color.var == "Both") {
    data$col <- paste(data$sample, data$Order.q, sep = "-")
    color_col <- "col"
  }

  # Base plot
  if (type == 2) {
    p <- ggplot(
      data,
      aes(
        x = m,
        y = SC,
        color = .data[[color_col]],
        linetype = Method
      )
    )
  } else {
    p <- ggplot(
      data,
      aes(
        x = if (type == 3) SC else m,
        y = if (type == 2) SC else qD,
        color = .data[[color_col]],
        linetype = Method
      )
    )
  }

  p <- p +
    geom_line(linewidth = 1.5) +
    geom_point(data = observed_endpoint, aes(shape = sample), size = 5) +
    scale_colour_manual(values = cbPalette) +
    labs(x = x_label, y = y_label) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      text = element_text(size = 12),
      legend.title = element_blank()
    )

  # Add confidence interval ribbons
  if (se) {
    if (type == 2) {
      p <- p +
        geom_ribbon(
          aes(ymin = SC.LCL, ymax = SC.UCL, fill = .data[[color_col]]),
          alpha = 0.2,
          color = NA
        )
    } else {
      p <- p +
        geom_ribbon(
          aes(ymin = qD.LCL, ymax = qD.UCL, fill = .data[[color_col]]),
          alpha = 0.2,
          color = NA
        )
    }
  }

  # Faceting
  if (facet.var == "Order.q") {
    p <- p +
      facet_wrap(
        ~Order.q,
        labeller = labeller(
          Order.q = c(`q = 0` = "q = 0", `q = 1` = "q = 1", `q = 2` = "q = 2")
        )
      )
  }

  if (facet.var == "Assemblage") {
    p <- p + facet_wrap(~sample, nrow = 1)
  }

  if (facet.var == "Both") {
    p <- p + facet_wrap(~ sample + Order.q)
  }

  if (grey) {
    p <- p + theme_bw()
  }

  return(p)
}
