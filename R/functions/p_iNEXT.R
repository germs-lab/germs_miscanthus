#' Parallel iNEXT Computation
#'
#' @description
#' Computes iNEXT diversity estimates in parallel across samples.
#' This is a parallelized wrapper around iNEXT::iNEXT() that processes
#' multiple samples simultaneously using the future framework.
#'
#' @param x A numeric matrix or data.frame where rows are samples and columns are species/OTUs.
#'   Row names should be sample identifiers.
#' @param q Numeric vector of diversity orders (default: c(0, 1, 2))
#' @param datatype Type of data: "abundance" or "incidence_freq" (default: "abundance")
#' @param endpoint Extrapolation endpoint (default: 2x maximum library size)
#' @param knots Number of knots for rarefaction/extrapolation curve (default: 40)
#' @param conf Confidence level for intervals (default: 0.95)
#' @param nboot Number of bootstrap replicates (default: 100)
#' @param max.cores Logical; if TRUE, uses all available cores minus 1 (default: FALSE)
#' @param nCores Number of cores to use if max.cores=FALSE (default: 1)
#' @param combine Logical; if TRUE, combines results into a single iNEXT object
#'   compatible with ggiNEXT(). If FALSE, returns a list of individual iNEXT objects (default: TRUE)
#' @param plan_strategy Future plan strategy: "multisession", "multicore", or "sequential" (default: "multisession")
#' @param verbose Logical; print progress messages (default: TRUE)
#' @param ... Additional arguments (currently unused)
#'
#' @return
#' If combine=TRUE: A single iNEXT object that can be plotted with ggiNEXT()
#' If combine=FALSE: A named list where each element is an iNEXT object for one sample
#'
#' @details
#' This function parallelizes iNEXT computation at the sample level. Each sample
#' is processed independently by a worker process. The bottleneck in iNEXT is
#' the bootstrap computation (controlled by nboot), which runs within each worker.
#'
#' The function automatically:
#' - Removes empty samples (rows with zero abundance)
#' - Sets up parallel workers using the future framework
#' - Manages memory with garbage collection
#' - Combines results into a format compatible with ggiNEXT (if combine=TRUE)
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' library(iNEXT)
#'
#' # Extract OTU table from phyloseq object
#' otu_mat <- as.matrix(otu_table(physeq_obj))
#'
#' # Run parallel iNEXT
#' results <- p_iNEXT(otu_mat, q = c(0, 1, 2), nCores = 4)
#'
#' # Plot with ggiNEXT
#' ggiNEXT(results, type = 1, facet.var = "Order.q")
#'
#' # Or get individual sample results
#' results_list <- p_iNEXT(otu_mat, combine = FALSE, nCores = 4)
#' }
#'
#' @export
p_iNEXT <- function(
  x,
  q = c(0, 1, 2),
  datatype = "abundance",
  endpoint = NULL,
  knots = 40,
  conf = 0.95,
  nboot = 100,
  max.cores = FALSE,
  nCores = 1,
  combine = TRUE,
  plan_strategy = "multisession",
  verbose = TRUE,
  ...
) {
  # Load required packages
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required. Please install it.")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required. Please install it.")
  }
  if (!requireNamespace("iNEXT", quietly = TRUE)) {
    stop("Package 'iNEXT' is required. Please install it.")
  }

  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  # Convert to matrix
  x <- as.matrix(x)

  # Ensure row names exist
  if (is.null(rownames(x))) {
    rownames(x) <- paste0("Sample", seq_len(nrow(x)))
    if (verbose) {
      message(
        "No row names found. Assigning default names: Sample1, Sample2, ..."
      )
    }
  }

  # Calculate library sizes and species counts
  library_size <- rowSums(x)
  species_num <- apply(x, 1, function(row) sum(row > 0))

  # Remove empty samples
  if (any(species_num <= 0)) {
    empty_samples <- rownames(x)[species_num <= 0]
    if (verbose) {
      message(sprintf(
        "Removing %d empty sample(s): %s",
        length(empty_samples),
        paste(empty_samples, collapse = ", ")
      ))
    }
    x <- x[species_num > 0, , drop = FALSE]
    library_size <- library_size[species_num > 0]
    species_num <- species_num[species_num > 0]
  }

  nr <- nrow(x)
  if (nr == 0) {
    stop("No samples remaining after removing empty samples")
  }

  # Set endpoint if not specified
  if (is.null(endpoint)) {
    endpoint <- max(library_size) * 2
    if (verbose) {
      message(sprintf("Endpoint set to 2x max library size: %d", endpoint))
    }
  }

  # Determine number of cores
  available_cores <- parallelly::availableCores()
  mc <- if (max.cores) {
    max(1, available_cores - 1)
  } else {
    min(nCores, available_cores)
  }

  if (verbose) {
    message(sprintf(
      "Using %d core(s). Available cores: %d",
      mc,
      available_cores
    ))
    message(sprintf(
      "Processing %d sample(s) with q = [%s]",
      nr,
      paste(q, collapse = ", ")
    ))
  }

  # Set up parallel plan
  if (mc > 1 && nr > 1) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    if (plan_strategy == "multisession") {
      future::plan(future::multisession, workers = mc)
    } else if (plan_strategy == "multicore") {
      future::plan(future::multicore, workers = mc)
    } else if (plan_strategy == "sequential") {
      future::plan(future::sequential)
    } else {
      stop(
        "Invalid plan_strategy. Choose 'multisession', 'multicore', or 'sequential'"
      )
    }
  } else {
    if (verbose && (mc == 1 || nr == 1)) {
      message("Running sequentially (only 1 core or 1 sample)")
    }
    future::plan(future::sequential)
  }

  # Parallel computation
  if (verbose) {
    message("Starting parallel iNEXT computation...")
  }

  out <- future.apply::future_lapply(
    seq_len(nr),
    function(i) {
      # Force garbage collection to free memory
      gc(verbose = FALSE, full = FALSE)

      # Run iNEXT for this sample
      tryCatch(
        {
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
        error = function(e) {
          list(
            error = TRUE,
            sample_name = rownames(x)[i],
            message = as.character(e)
          )
        }
      )
    },
    future.seed = TRUE,
    future.scheduling = 2.0 # Dynamic load balancing
  )

  # Name the results
  names(out) <- rownames(x)

  # Check for errors
  errors <- sapply(out, function(x) !is.null(x$error) && x$error)
  if (any(errors)) {
    error_samples <- names(out)[errors]
    error_messages <- sapply(out[errors], function(x) x$message)
    warning(sprintf(
      "iNEXT failed for %d sample(s): %s\nMessages: %s",
      sum(errors),
      paste(error_samples, collapse = ", "),
      paste(error_messages, collapse = "; ")
    ))
    # Remove failed samples
    out <- out[!errors]
  }

  if (length(out) == 0) {
    stop("iNEXT failed for all samples")
  }

  if (verbose) {
    message(sprintf(
      "Completed processing %d sample(s) successfully",
      length(out)
    ))
  }

  # Return results
  if (combine) {
    if (verbose) {
      message("Combining results into single iNEXT object...")
    }
    combined <- combine_iNEXT_list(out)
    return(combined)
  } else {
    return(out)
  }
}


#' Combine List of iNEXT Objects
#'
#' @description
#' Combines a list of individual iNEXT objects (one per sample) into a single
#' iNEXT object that is compatible with ggiNEXT() plotting.
#'
#' @param inext_list A named list of iNEXT objects
#'
#' @return A single iNEXT object combining all samples
#'
#' @details
#' This function merges the size_based and coverage_based results from multiple
#' iNEXT objects, preserving the structure expected by ggiNEXT().
#'
#' @keywords internal
combine_iNEXT_list <- function(inext_list) {
  if (!is.list(inext_list) || length(inext_list) == 0) {
    stop("inext_list must be a non-empty list")
  }

  # Extract components from each iNEXT object
  all_size_based <- lapply(names(inext_list), function(sample_name) {
    obj <- inext_list[[sample_name]]
    if (is.null(obj$iNextEst$size_based)) {
      return(NULL)
    }
    # Update Assemblage name to sample name
    df <- obj$iNextEst$size_based
    df$Assemblage <- sample_name
    return(df)
  })
  all_size_based <- do.call(
    rbind,
    all_size_based[!sapply(all_size_based, is.null)]
  )

  all_coverage_based <- lapply(names(inext_list), function(sample_name) {
    obj <- inext_list[[sample_name]]
    if (is.null(obj$iNextEst$coverage_based)) {
      return(NULL)
    }
    df <- obj$iNextEst$coverage_based
    df$Assemblage <- sample_name
    return(df)
  })
  all_coverage_based <- do.call(
    rbind,
    all_coverage_based[!sapply(all_coverage_based, is.null)]
  )

  # Combine AsyEst (asymptotic estimates)
  all_asy_est <- lapply(names(inext_list), function(sample_name) {
    obj <- inext_list[[sample_name]]
    if (is.null(obj$AsyEst)) {
      return(NULL)
    }
    df <- obj$AsyEst
    if (!"Assemblage" %in% colnames(df)) {
      df$Assemblage <- sample_name
    } else {
      df$Assemblage <- sample_name
    }
    return(df)
  })
  all_asy_est <- do.call(rbind, all_asy_est[!sapply(all_asy_est, is.null)])

  # Combine DataInfo
  all_data_info <- lapply(names(inext_list), function(sample_name) {
    obj <- inext_list[[sample_name]]
    if (is.null(obj$DataInfo)) {
      return(NULL)
    }
    df <- obj$DataInfo
    df$Assemblage <- sample_name
    return(df)
  })
  all_data_info <- do.call(
    rbind,
    all_data_info[!sapply(all_data_info, is.null)]
  )
  rownames(all_data_info) <- all_data_info$Assemblage

  # Create combined iNEXT object
  combined <- list(
    DataInfo = all_data_info,
    iNextEst = list(
      size_based = all_size_based,
      coverage_based = all_coverage_based
    ),
    AsyEst = all_asy_est
  )

  class(combined) <- "iNEXT"
  return(combined)
}


#' Extract and Flatten iNEXT Results to Data Frame
#'
#' @description
#' Extracts size-based and coverage-based results from iNEXT object(s)
#' and returns a long-format data frame suitable for custom plotting.
#'
#' @param inext_obj Either a single iNEXT object or a list of iNEXT objects
#'
#' @return A data frame with columns for sample, diversity metrics, and metadata
#'
#' @details
#' This function is useful when you want to create custom plots using ggplot2
#' rather than using ggiNEXT() or gg_inext_custom().
#'
#' @examples
#' \dontrun{
#' results <- p_iNEXT(otu_mat, combine = FALSE)
#' df <- flatten_iNEXT_results(results)
#' # Now use df for custom ggplot2 plotting
#' }
#'
#' @export
flatten_iNEXT_results <- function(inext_obj) {
  if ("iNEXT" %in% class(inext_obj)) {
    # Single iNEXT object
    size_based <- inext_obj$iNextEst$size_based
    coverage_based <- inext_obj$iNextEst$coverage_based

    # Add prefix to distinguish size vs coverage based
    colnames(size_based) <- paste0("size_based.", colnames(size_based))
    colnames(coverage_based) <- paste0(
      "coverage_based.",
      colnames(coverage_based)
    )

    # Merge by Assemblage and Order.q
    result <- merge(
      size_based,
      coverage_based,
      by.x = c("size_based.Assemblage", "size_based.Order.q"),
      by.y = c("coverage_based.Assemblage", "coverage_based.Order.q"),
      all = TRUE,
      suffixes = c("", "")
    )

    # Rename for clarity
    result$sample <- result$size_based.Assemblage
    result$size_based.Assemblage <- NULL

    return(result)
  } else if (is.list(inext_obj)) {
    # List of iNEXT objects - combine first
    combined <- combine_iNEXT_list(inext_obj)
    return(flatten_iNEXT_results(combined))
  } else {
    stop("inext_obj must be an iNEXT object or list of iNEXT objects")
  }
}


##########################
# Plotting functions
##########################

#' Generate ggplot2 Color Palette
#'
#' @param g Number of colors to generate
#' @return Vector of hexadecimal color codes
#' @keywords internal
ggplotColors <- function(g) {
  d <- 360 / g
  h <- cumsum(c(15, rep(d, g - 1)))
  hcl(h = h, c = 100, l = 65)
}

#' Custom ggplot for iNEXT Results
#'
#' @description
#' Custom plotting function for flattened iNEXT results.
#' Similar to ggiNEXT but works with the flattened data frame format.
#'
#' @param inext_data Flattened iNEXT data frame (from flatten_iNEXT_results)
#' @param type Plot type: 1 (size-based), 2 (sample completeness), 3 (coverage-based)
#' @param se Show confidence intervals
#' @param facet.var Faceting variable: "None", "Order.q", "Assemblage", "Both"
#' @param color.var Color variable: "None", "Order.q", "Assemblage", "Both"
#' @param grey Use grey theme
#'
#' @return A ggplot2 object
#'
#' @export
gg_inext_custom <- function(
  inext_data,
  type = 1,
  se = TRUE,
  facet.var = "None",
  color.var = "Assemblage",
  grey = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required")
  }

  require(ggplot2)
  require(dplyr)

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

  # Apply facet.var -> color.var logic
  if (facet.var == "Order.q") {
    color.var <- "Assemblage"
  }
  if (facet.var == "Assemblage") {
    color.var <- "Order.q"
  }

  # Handle color.var
  if (color.var == "None") {
    if (length(unique(data$sample)) > 1 & length(unique(data$Order.q)) > 1) {
      warning("invalid color.var setting, changing to 'Both'")
      color.var <- "Both"
      data$col <- paste(data$sample, data$Order.q, sep = "-")
      observed_endpoint$col <- paste(
        observed_endpoint$sample,
        observed_endpoint$Order.q,
        sep = "-"
      )
      color_col <- "col"
    } else if (length(unique(data$sample)) > 1) {
      warning("invalid color.var setting, changing to 'Assemblage'")
      color.var <- "Assemblage"
      color_col <- "sample"
    } else if (length(unique(data$Order.q)) > 1) {
      warning("invalid color.var setting, changing to 'Order.q'")
      color.var <- "Order.q"
      color_col <- "Order.q"
    } else {
      color_col <- "sample"
    }
  } else if (color.var == "Order.q") {
    color_col <- "Order.q"
  } else if (color.var == "Assemblage") {
    if (length(unique(data$sample)) == 1) {
      warning("invalid color.var setting, changing to 'Order.q'")
      color_col <- "Order.q"
    } else {
      color_col <- "sample"
    }
  } else if (color.var == "Both") {
    if (length(unique(data$sample)) == 1) {
      warning("invalid color.var setting, changing to 'Order.q'")
      color_col <- "Order.q"
    } else {
      data$col <- paste(data$sample, data$Order.q, sep = "-")
      observed_endpoint$col <- paste(
        observed_endpoint$sample,
        observed_endpoint$Order.q,
        sep = "-"
      )
      color_col <- "col"
    }
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
          Order.q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2")
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
