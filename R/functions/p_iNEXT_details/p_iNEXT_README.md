# Parallel iNEXT Implementation Guide

This directory contains multiple implementations of parallel iNEXT computation with different trade-offs and use cases.

## Quick Start

**For most users**, use `p_iNEXT_improved.R`:

```r
source("R/functions/p_iNEXT_improved.R")

# Run parallel iNEXT
result <- p_iNEXT(
  x = otu_matrix,  # Samples as rows, species as columns
  q = c(0, 1, 2),
  nCores = 4,
  combine = TRUE   # Returns iNEXT object compatible with ggiNEXT()
)

# Plot with standard ggiNEXT
library(iNEXT)
ggiNEXT(result, type = 1, facet.var = "Order.q")
```

## File Overview

### 1. `p_iNEXT.R` (Original)
**Status**: Legacy, maintained for backward compatibility

**Features**:
- Sample-level parallelization
- Returns list format: `$inext_results`
- Requires custom plotting with `gg_inext_custom()`

**Use when**:
- Existing code depends on this format
- Already have custom plotting workflows

### 2. `p_iNEXT_improved.R` (Recommended)
**Status**: **RECOMMENDED for all new code**

**Key Features**:
- ✅ Direct `ggiNEXT()` compatibility (set `combine = TRUE`)
- ✅ Robust error handling
- ✅ Comprehensive documentation
- ✅ Helper functions for data manipulation
- ✅ Flexible output formats
- ✅ Better parallel plan management

**Use when**:
- Starting a new project
- Want standard iNEXT/ggiNEXT workflow
- Need better error handling
- Want well-documented code

**Example**:
```r
source("R/functions/p_iNEXT_improved.R")

# Standard workflow
result <- p_iNEXT(otu_mat, q = c(0, 1, 2), nCores = 4)
ggiNEXT(result, type = 1)  # Works directly!

# Or get individual samples
samples <- p_iNEXT(otu_mat, combine = FALSE, nCores = 4)

# Or flatten for custom ggplot2
df <- flatten_iNEXT_results(result)
ggplot(df, aes(x = size_based.m, y = size_based.qD, color = sample)) +
  geom_line()
```

### 3. `p_iNEXT_deep_parallel.R` (Experimental)
**Status**: Experimental, not fully implemented

**Concept**: Two-level parallelization (samples AND bootstrap)

**Current State**:
- Sample-level parallelization implemented
- Bootstrap-level parallelization documented but not implemented
- Requires modifying iNEXT source code for full implementation

**Use when**:
- You have >50 samples AND >500 bootstrap replicates
- You have >16 cores available
- You're willing to maintain custom iNEXT modifications
- Standard parallelization is insufficient after profiling

**Note**: For most users, the overhead of nested parallelization outweighs benefits. Use `p_iNEXT_improved.R` instead.

### 4. `p_iNEXT_comparison.md`
Detailed comparison of approaches, migration guide, and recommendations.

### 5. `p_iNEXT_examples.R`
Comprehensive examples showing all features and use cases.

## Parallelization Strategy Explained

### What Gets Parallelized?

All versions use **sample-level parallelization**:

```
Sample 1 → Worker 1 → iNEXT() [includes bootstrap] → Result 1
Sample 2 → Worker 2 → iNEXT() [includes bootstrap] → Result 2
Sample 3 → Worker 3 → iNEXT() [includes bootstrap] → Result 3
...
```

Each worker runs a complete `iNEXT()` computation for one sample, including all bootstrap replicates.

### Why Not Deeper Parallelization?

The bootstrap loop in iNEXT (the main bottleneck) is:
1. Already vectorized with `apply()`
2. Has low per-iteration computation time
3. Would have high parallelization overhead
4. Tightly coupled with other iNEXT internals

**Result**: Sample-level parallelization provides the best balance of:
- Speedup
- Simplicity
- Memory efficiency
- Maintainability

### Performance Expectations

With `nCores = 4` and `nboot = 100`:
- **Ideal speedup**: 4x
- **Realistic speedup**: 2.5-3.5x (due to overhead)
- **Memory usage**: ~4x (one copy per worker)

Performance scales best when:
- Number of samples ≥ number of cores
- Samples have similar processing time
- Sufficient RAM per core (~2-4GB minimum)
- Bootstrap count is reasonable (50-200)

## Choosing the Right Tool

### Decision Tree

```
Are you starting a new project?
├─ YES → Use p_iNEXT_improved.R
│
└─ NO → Do you need ggiNEXT() compatibility?
    ├─ YES → Migrate to p_iNEXT_improved.R (see comparison.md)
    │
    └─ NO → Do you have existing custom plotting?
        ├─ YES → Keep p_iNEXT.R (legacy)
        │
        └─ NO → Consider p_iNEXT_improved.R for better features

Do you have >50 samples AND >500 bootstrap?
└─ Consider profiling first, then maybe deep parallelization
```

## Common Workflows

### Workflow 1: Standard Analysis (Recommended)

```r
library(iNEXT)
source("R/functions/p_iNEXT_improved.R")

# Run analysis
result <- p_iNEXT(
  x = otu_matrix,
  q = c(0, 1, 2),
  nCores = 4,
  nboot = 100,
  combine = TRUE
)

# Plot with ggiNEXT
p1 <- ggiNEXT(result, type = 1, facet.var = "Order.q")
p2 <- ggiNEXT(result, type = 2)
p3 <- ggiNEXT(result, type = 3, facet.var = "Order.q")

print(p1)
print(p2)
print(p3)
```

### Workflow 2: Custom Plotting

```r
source("R/functions/p_iNEXT_improved.R")

# Run analysis
result <- p_iNEXT(otu_matrix, q = c(0, 1, 2), nCores = 4)

# Flatten to data frame
df <- flatten_iNEXT_results(result)

# Custom ggplot2
library(ggplot2)
ggplot(df, aes(x = size_based.m, y = size_based.qD, 
               color = sample, group = sample)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = size_based.qD.LCL, ymax = size_based.qD.UCL,
                  fill = sample), alpha = 0.2) +
  facet_wrap(~size_based.Order.q) +
  labs(x = "Sample Size", y = "Diversity", title = "Custom Plot") +
  theme_bw()
```

### Workflow 3: Phyloseq Integration

```r
library(phyloseq)
source("R/functions/p_iNEXT_improved.R")

# Extract from phyloseq
otu_mat <- as.matrix(t(otu_table(physeq)))  # Samples as rows

# Calculate endpoint
endpoint <- max(rowSums(otu_mat)) * 2

# Run iNEXT
result <- p_iNEXT(
  x = otu_mat,
  q = c(0, 1, 2),
  endpoint = endpoint,
  nCores = 4,
  verbose = TRUE
)

# Plot
ggiNEXT(result, type = 1, facet.var = "Order.q") +
  theme_bw() +
  labs(title = "Rarefaction Curves")
```

### Workflow 4: Processing Many Datasets

```r
source("R/functions/p_iNEXT_improved.R")

# Process multiple phyloseq objects
results <- lapply(phyloseq_list, function(ps) {
  otu_mat <- as.matrix(t(otu_table(ps)))
  p_iNEXT(otu_mat, q = c(0, 1, 2), nCores = 4, combine = TRUE)
})

# Plot each
plots <- lapply(seq_along(results), function(i) {
  ggiNEXT(results[[i]], type = 1) +
    ggtitle(names(phyloseq_list)[i])
})

# Display
library(patchwork)
wrap_plots(plots, ncol = 2)
```

### Workflow 5: Batch Processing Large Datasets

```r
source("R/functions/p_iNEXT_improved.R")

# For very large datasets, process in batches
batch_size <- 50
n_samples <- nrow(large_otu_matrix)
batches <- split(seq_len(n_samples), ceiling(seq_len(n_samples) / batch_size))

# Process each batch
all_results <- lapply(batches, function(idx) {
  cat("Processing batch:", min(idx), "-", max(idx), "\n")
  p_iNEXT(
    large_otu_matrix[idx, ],
    q = c(0, 1, 2),
    nCores = 8,
    combine = FALSE,  # Keep as list for now
    verbose = FALSE
  )
})

# Flatten and combine all
all_results <- unlist(all_results, recursive = FALSE)
final_result <- combine_iNEXT_list(all_results)

# Now can use with ggiNEXT
ggiNEXT(final_result, type = 1)
```

## Performance Tuning

### Choosing Number of Cores

```r
# Check available cores
parallelly::availableCores()

# General rule: use N-1 cores (leave one for system)
nCores <- parallelly::availableCores() - 1

# For small datasets (<10 samples), sequential may be faster
if (nrow(otu_mat) < 10) {
  nCores <- 1
}

# For memory-constrained systems
# Estimate: ~2-4GB RAM per core needed
available_ram_gb <- 16  # Your system RAM
cores_by_ram <- floor(available_ram_gb / 3)
nCores <- min(nCores, cores_by_ram)
```

### Choosing Bootstrap Count

```r
# Quick exploratory analysis
nboot <- 50

# Standard analysis
nboot <- 100

# Publication-quality
nboot <- 200

# Note: Higher bootstrap = longer computation but more stable CI
# Diminishing returns beyond ~200
```

### Monitoring Performance

```r
# Time your analysis
system.time({
  result <- p_iNEXT(otu_mat, q = c(0, 1, 2), nCores = 4)
})

# Compare sequential vs parallel
time_seq <- system.time({
  r1 <- p_iNEXT(otu_mat, nCores = 1, verbose = FALSE)
})

time_par <- system.time({
  r2 <- p_iNEXT(otu_mat, nCores = 4, verbose = FALSE)
})

cat("Speedup:", time_seq["elapsed"] / time_par["elapsed"], "x\n")
```

## Troubleshooting

### Issue: Out of Memory

**Solution 1**: Reduce cores
```r
result <- p_iNEXT(otu_mat, nCores = 2)  # Instead of 4
```

**Solution 2**: Reduce bootstrap
```r
result <- p_iNEXT(otu_mat, nboot = 50)  # Instead of 100
```

**Solution 3**: Process in batches (see Workflow 5)

### Issue: Slow Performance

**Check 1**: Are you using enough cores?
```r
p_iNEXT(otu_mat, nCores = parallelly::availableCores() - 1)
```

**Check 2**: Is parallelization overhead the issue? (few samples)
```r
# If < 5-10 samples, sequential might be faster
if (nrow(otu_mat) < 10) {
  result <- p_iNEXT(otu_mat, nCores = 1)
}
```

**Check 3**: Profile to find bottleneck
```r
library(profvis)
profvis({
  result <- p_iNEXT(otu_mat, nCores = 4)
})
```

### Issue: Failed Samples

The improved version handles this gracefully:
```r
result <- p_iNEXT(otu_mat, nCores = 4, verbose = TRUE)
# Will show which samples failed and continue with successful ones
```

### Issue: Not Compatible with ggiNEXT

Make sure `combine = TRUE`:
```r
result <- p_iNEXT(otu_mat, combine = TRUE)  # Default
class(result)  # Should be "iNEXT"
```

## Migration from Original

If you have code using original `p_iNEXT.R`:

```r
# OLD:
source("R/functions/p_iNEXT.R")
result <- p_iNEXT(otu_mat, nCores = 4)
samples <- result$inext_results
# ... custom plotting ...

# NEW (option 1 - equivalent):
source("R/functions/p_iNEXT_improved.R")
samples <- p_iNEXT(otu_mat, nCores = 4, combine = FALSE)
# ... keep custom plotting ...

# NEW (option 2 - recommended):
source("R/functions/p_iNEXT_improved.R")
result <- p_iNEXT(otu_mat, nCores = 4, combine = TRUE)
ggiNEXT(result, type = 1)  # Use standard plotting
```

See `p_iNEXT_comparison.md` for detailed migration guide.

## Testing

Run the examples:
```r
source("R/functions/p_iNEXT_examples.R")
```

This will test all major features with the spider dataset.

## References

- Original iNEXT package: https://github.com/AnneChao/iNEXT
- Chao, A., et al. (2014). Rarefaction and extrapolation with Hill numbers. Methods in Ecology and Evolution, 5(5), 451-456.
- future package: https://future.futureverse.org/

## Support

For issues or questions:
1. Check `p_iNEXT_comparison.md` for detailed explanations
2. Run `p_iNEXT_examples.R` to see working examples
3. Review troubleshooting section above
4. Check original iNEXT documentation: `?iNEXT::iNEXT`

## License

Same as the parent repository.
