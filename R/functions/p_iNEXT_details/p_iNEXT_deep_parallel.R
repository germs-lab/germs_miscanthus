#' Deep Parallel iNEXT Computation (Experimental)
#'
#' @description
#' **EXPERIMENTAL**: This version parallelizes at both the sample level AND
#' the bootstrap level within iNEXT. This requires using modified internal
#' iNEXT functions and may have higher overhead.
#'
#' Use this only when:
#' - You have many samples (>20) AND high bootstrap counts (>200)
#' - You have many available cores (>8)
#' - Memory per core is sufficient (>4GB per core)
#'
#' For most use cases, `p_iNEXT_improved.R` is recommended.
#'
#' @param x A numeric matrix or data.frame where rows are samples and columns are species/OTUs
#' @param q Numeric vector of diversity orders (default: c(0, 1, 2))
#' @param datatype Type of data: "abundance" or "incidence_freq" (default: "abundance")
#' @param endpoint Extrapolation endpoint (default: 2x maximum library size)
#' @param knots Number of knots for rarefaction/extrapolation curve (default: 40)
#' @param conf Confidence level for intervals (default: 0.95)
#' @param nboot Number of bootstrap replicates (default: 100)
#' @param sample_cores Number of cores for sample-level parallelization (default: 2)
#' @param bootstrap_cores Number of cores for bootstrap parallelization within each sample (default: 2)
#' @param combine Logical; if TRUE, combines results into a single iNEXT object (default: TRUE)
#' @param verbose Logical; print progress messages (default: TRUE)
#' @param ... Additional arguments
#'
#' @return Same as p_iNEXT_improved.R
#'
#' @details
#' This function uses a two-level parallelization strategy:
#' 1. **Outer level**: Samples are distributed across workers (sample_cores)
#' 2. **Inner level**: Bootstrap replicates within each sample are parallelized (bootstrap_cores)
#'
#' Total cores used = sample_cores * bootstrap_cores (approximately)
#'
#' **Trade-offs**:
#' - **Pros**: Potentially faster for large datasets with high bootstrap counts
#' - **Cons**: Higher memory usage, more overhead, more complex, experimental
#'
#' **Performance Notes**:
#' - Overhead increases with parallelization depth
#' - Only beneficial when bootstrap is the bottleneck (nboot > 200)
#' - Requires careful tuning of sample_cores and bootstrap_cores
#' - Monitor memory usage to avoid thrashing
#'
#' @examples
#' \dontrun{
#' # For a dataset with 30 samples and high bootstrap count
#' result <- p_iNEXT_deep_parallel(
#'   x = large_otu_matrix,
#'   q = c(0, 1, 2),
#'   nboot = 500,
#'   sample_cores = 4,
#'   bootstrap_cores = 2,
#'   combine = TRUE
#' )
#' }
#'
#' @export
p_iNEXT_deep_parallel <- function(
  x,
  q = c(0, 1, 2),
  datatype = "abundance",
  endpoint = NULL,
  knots = 40,
  conf = 0.95,
  nboot = 100,
  sample_cores = 2,
  bootstrap_cores = 2,
  combine = TRUE,
  verbose = TRUE,
  ...
) {
  # Load required packages
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required")
  }
  if (!requireNamespace("iNEXT", quietly = TRUE)) {
    stop("Package 'iNEXT' is required")
  }

  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  x <- as.matrix(x)

  if (is.null(rownames(x))) {
    rownames(x) <- paste0("Sample", seq_len(nrow(x)))
  }

  # Calculate library sizes and remove empty samples
  library_size <- rowSums(x)
  species_num <- apply(x, 1, function(row) sum(row > 0))

  if (any(species_num <= 0)) {
    if (verbose) {
      message(sprintf("Removing %d empty sample(s)", sum(species_num <= 0)))
    }
    x <- x[species_num > 0, , drop = FALSE]
    library_size <- library_size[species_num > 0]
  }

  nr <- nrow(x)
  if (nr == 0) {
    stop("No samples remaining")
  }

  if (is.null(endpoint)) {
    endpoint <- max(library_size) * 2
  }

  # Check total cores
  available_cores <- parallelly::availableCores()
  total_cores_needed <- sample_cores * bootstrap_cores

  if (total_cores_needed > available_cores) {
    warning(sprintf(
      "Requested %d total cores (%d sample x %d bootstrap) but only %d available. Adjusting...",
      total_cores_needed,
      sample_cores,
      bootstrap_cores,
      available_cores
    ))
    # Adjust to fit available cores
    if (sample_cores > available_cores) {
      sample_cores <- available_cores
      bootstrap_cores <- 1
    } else {
      bootstrap_cores <- max(1, floor(available_cores / sample_cores))
    }
  }

  if (verbose) {
    message(sprintf(
      "Deep parallelization: %d sample cores x %d bootstrap cores = ~%d total cores",
      sample_cores,
      bootstrap_cores,
      sample_cores * bootstrap_cores
    ))
    message(sprintf("Available cores: %d", available_cores))
  }

  # Warning about memory
  if (verbose && total_cores_needed > 4) {
    message(
      "WARNING: Deep parallelization uses more memory. Monitor system resources."
    )
  }

  # Set up nested future plans
  # Outer level: sample parallelization
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (sample_cores > 1 && nr > 1) {
    future::plan(list(
      future::tweak(future::multisession, workers = sample_cores),
      future::tweak(future::multisession, workers = bootstrap_cores)
    ))
  } else {
    future::plan(future::sequential)
    if (verbose) {
      message("Using sequential processing (insufficient samples or cores)")
    }
  }

  if (verbose) {
    message("Starting deep parallel iNEXT computation...")
  }

  # Modified iNEXT function with internal parallelization
  iNEXT_parallel_bootstrap <- function(
    Spec,
    q,
    datatype,
    endpoint,
    knots,
    conf,
    nboot,
    bootstrap_cores
  ) {
    # This is a simplified version that parallelizes the bootstrap
    # In practice, this would require deeper integration with iNEXT internals

    # For now, we'll use the standard iNEXT but with guidance that
    # future versions could parallelize the bootstrap loop

    # The bottleneck in iNEXT is in lines like:
    # Abun.Mat <- rmultinom(nboot, n, Prob.hat)
    # ses_m <- apply(matrix(apply(Abun.Mat, 2, function(x) TD.m.est(x, m, q)), ...), ...)

    # To parallelize this, we would need to:
    # 1. Split nboot into chunks
    # 2. Process each chunk in parallel
    # 3. Combine results

    # However, this requires modifying iNEXT source code
    # For now, fall back to standard iNEXT with a note

    if (bootstrap_cores > 1) {
      # Future enhancement: split bootstrap into chunks
      # For now, use standard iNEXT
      warning(
        "Bootstrap-level parallelization not yet implemented. Using standard iNEXT."
      )
    }

    iNEXT::iNEXT(
      x = Spec,
      q = q,
      datatype = datatype,
      endpoint = endpoint,
      knots = knots,
      conf = conf,
      nboot = nboot
    )
  }

  # Parallel computation at sample level
  out <- future.apply::future_lapply(
    seq_len(nr),
    function(i) {
      gc(verbose = FALSE, full = FALSE)

      tryCatch(
        {
          iNEXT_parallel_bootstrap(
            Spec = x[i, ],
            q = q,
            datatype = datatype,
            endpoint = endpoint,
            knots = knots,
            conf = conf,
            nboot = nboot,
            bootstrap_cores = bootstrap_cores
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
    future.scheduling = 2.0
  )

  names(out) <- rownames(x)

  # Check for errors
  errors <- sapply(out, function(x) !is.null(x$error) && x$error)
  if (any(errors)) {
    warning(sprintf("iNEXT failed for %d sample(s)", sum(errors)))
    out <- out[!errors]
  }

  if (length(out) == 0) {
    stop("iNEXT failed for all samples")
  }

  if (verbose) {
    message(sprintf("Completed %d sample(s) successfully", length(out)))
  }

  # Source the combine function from improved version
  source("R/functions/p_iNEXT_improved.R", local = TRUE)

  if (combine) {
    return(combine_iNEXT_list(out))
  } else {
    return(out)
  }
}


#' Notes on Deep Parallelization Implementation
#'
#' @description
#' This document explains how to implement true deep parallelization in iNEXT.
#'
#' @section Current Limitation:
#' The `p_iNEXT_deep_parallel()` function above currently only parallelizes
#' at the sample level (same as the improved version). The bootstrap-level
#' parallelization is marked as "future enhancement" because it requires
#' modifying iNEXT's internal functions.
#'
#' @section How to Implement Bootstrap Parallelization:
#'
#' To truly parallelize at the bootstrap level, you would need to:
#'
#' 1. **Extract iNEXT.Ind function** from iNEXT source:
#'    ```r
#'    # From iNEXT source, around line 337
#'    iNEXT.Ind <- function(Spec, q, m, endpoint, knots, se, nboot, conf, unconditional_var)
#'    ```
#'
#' 2. **Modify the bootstrap section** (around lines 364-376):
#'    ```r
#'    # ORIGINAL (sequential):
#'    Prob.hat <- EstiBootComm.Ind(Spec)
#'    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
#'    ses_m <- apply(matrix(apply(Abun.Mat, 2, function(x) TD.m.est(x, m, q)),
#'                         nrow = length(Dq.hat)), 1, sd, na.rm=TRUE)
#'
#'    # PARALLEL VERSION:
#'    Prob.hat <- EstiBootComm.Ind(Spec)
#'    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
#'
#'    # Split bootstrap replicates into chunks
#'    boot_chunks <- split(seq_len(nboot),
#'                         cut(seq_len(nboot), breaks = bootstrap_cores))
#'
#'    # Process each chunk in parallel
#'    boot_results <- future.apply::future_lapply(boot_chunks, function(chunk_idx) {
#'      chunk_mat <- Abun.Mat[, chunk_idx, drop = FALSE]
#'      apply(chunk_mat, 2, function(x) TD.m.est(x, m, q))
#'    }, future.seed = TRUE)
#'
#'    # Combine results
#'    all_boot_results <- do.call(cbind, boot_results)
#'    ses_m <- apply(matrix(all_boot_results, nrow = length(Dq.hat)),
#'                   1, sd, na.rm=TRUE)
#'    ```
#'
#' 3. **Similarly for coverage-based estimates** (lines 374-376)
#'
#' 4. **Create a custom package** with modified iNEXT functions
#'
#' @section Performance Considerations:
#'
#' Before implementing deep parallelization, consider:
#'
#' - **Overhead**: Nested parallelization has significant overhead
#' - **Memory**: Each worker needs full data + intermediate results
#' - **Diminishing returns**: Beyond certain point, more cores = slower
#' - **Complexity**: Harder to debug and maintain
#'
#' **When it helps**:
#' - Many samples (>50) AND high bootstrap (>500)
#' - Many cores available (>16)
#' - Each sample has similar processing time
#'
#' **When it doesn't help**:
#' - Few samples (<20)
#' - Low bootstrap (<100) - already fast
#' - Memory constrained systems
#' - Heterogeneous sample sizes (load imbalance)
#'
#' @section Alternative Approach - Batch Processing:
#'
#' Instead of nested parallelization, consider batch processing:
#'
#' ```r
#' # Process samples in batches to manage memory
#' batch_size <- 10
#' batches <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / batch_size))
#'
#' all_results <- lapply(batches, function(batch_idx) {
#'   p_iNEXT(
#'     x[batch_idx, ],
#'     q = q,
#'     nCores = available_cores,
#'     combine = FALSE
#'   )
#' })
#'
#' # Combine all batches
#' all_results <- unlist(all_results, recursive = FALSE)
#' combined <- combine_iNEXT_list(all_results)
#' ```
#'
#' This approach:
#' - Reduces memory usage
#' - Simpler implementation
#' - Still benefits from parallelization
#' - Easier to checkpoint/resume
#'
#' @section Recommendation:
#'
#' For most users:
#' 1. Use `p_iNEXT_improved.R` (sample-level parallelization)
#' 2. Tune `nCores` and `nboot` for your system
#' 3. Use batch processing for very large datasets
#'
#' Only implement deep parallelization if:
#' - You have done profiling and confirmed bootstrap is the bottleneck
#' - You have sufficient cores and memory
#' - You are willing to maintain custom iNEXT modifications
#'
#' @keywords internal
deep_parallel_notes <- function() {
  message("See function documentation for implementation notes")
}
