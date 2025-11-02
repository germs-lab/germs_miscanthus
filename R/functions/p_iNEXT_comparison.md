# p_iNEXT Improvements - Comparison and Usage Guide

## Summary of Changes

### Original `p_iNEXT.R` Issues
1. **Limited output format**: Only returns list structure, not directly compatible with `ggiNEXT()`
2. **Poor error handling**: No graceful handling of failed samples
3. **Suboptimal parallel configuration**: Hardcoded plan setup, no cleanup
4. **Missing features**: No option to combine results, limited documentation
5. **Code duplication**: Has plotting code mixed with computation

### Improved `p_iNEXT_improved.R` Features

#### 1. **Better Output Compatibility**
- **New `combine` parameter**: Merge results into single iNEXT object compatible with `ggiNEXT()`
- **Flexible output**: Choose between list of iNEXT objects or combined object
- **Helper functions**: `combine_iNEXT_list()` and `flatten_iNEXT_results()` for data manipulation

#### 2. **Enhanced Parallel Processing**
- **Configurable strategy**: Choose between `multisession`, `multicore`, or `sequential`
- **Proper plan management**: Restores previous plan on exit
- **Dynamic load balancing**: Uses `future.scheduling = 2.0` for better distribution
- **Smart core allocation**: Only parallelizes when it makes sense (>1 core, >1 sample)

#### 3. **Robust Error Handling**
- **Per-sample error catching**: Failed samples don't crash entire computation
- **Informative warnings**: Reports which samples failed and why
- **Continues processing**: Returns successful samples even if some fail

#### 4. **Better Memory Management**
- **Explicit garbage collection**: Reduces memory footprint in workers
- **Efficient data structures**: Avoids unnecessary copies

#### 5. **Comprehensive Documentation**
- **Full roxygen2 docs**: Complete parameter descriptions
- **Usage examples**: Clear examples for common use cases
- **Return value documentation**: Explains output structure

#### 6. **Cleaner Architecture**
- **Separated concerns**: Plotting code remains separate from computation
- **Modular design**: Helper functions can be used independently
- **Input validation**: Checks inputs before expensive computation

## Usage Comparison

### Original Usage
```r
# Original p_iNEXT
result <- p_iNEXT(
  x = otu_matrix,
  q = c(0, 1, 2),
  max.cores = FALSE,
  nCores = 4
)

# Returns list with $inext_results
sample_list <- result$inext_results

# Cannot directly use ggiNEXT() - need to manually extract and combine
# Would need custom plotting with gg_inext_custom()
```

### Improved Usage - Combined Output (Default)
```r
# Improved p_iNEXT with combined output (compatible with ggiNEXT)
result <- p_iNEXT(
  x = otu_matrix,
  q = c(0, 1, 2),
  nCores = 4,
  combine = TRUE  # Default
)

# Result is an iNEXT object - directly compatible!
library(iNEXT)
ggiNEXT(result, type = 1, facet.var = "Order.q", color.var = "Assemblage")
ggiNEXT(result, type = 2)
ggiNEXT(result, type = 3)
```

### Improved Usage - List Output
```r
# Get list of individual iNEXT objects
result_list <- p_iNEXT(
  x = otu_matrix,
  q = c(0, 1, 2),
  nCores = 4,
  combine = FALSE
)

# Access individual samples
sample1 <- result_list$Sample1
ggiNEXT(sample1, type = 1)

# Or combine later if needed
combined <- combine_iNEXT_list(result_list)
ggiNEXT(combined, type = 1)
```

### Improved Usage - Custom Plotting
```r
# For custom ggplot2 plots
result <- p_iNEXT(otu_matrix, q = c(0, 1, 2), nCores = 4)

# Flatten to data frame
df <- flatten_iNEXT_results(result)

# Custom ggplot
library(ggplot2)
ggplot(df, aes(x = size_based.m, y = size_based.qD, color = sample)) +
  geom_line() +
  facet_wrap(~size_based.Order.q) +
  theme_minimal()

# Or use the custom plotting function
gg_inext_custom(df, type = 1, facet.var = "Order.q")
```

## Performance Comparison

### Parallelization Strategy
Both versions use the **same parallelization approach**:
- Parallelize at the **sample level** (row-wise)
- Each worker runs one complete `iNEXT::iNEXT()` call
- Bootstrap computations (the bottleneck) run within each worker

**Why not deeper parallelization?**
- iNEXT's internal functions are tightly coupled
- Bootstrap loops are vectorized with `apply()` - already efficient
- Overhead of deeper parallelization would outweigh benefits
- Sample-level parallelization is simpler and more maintainable

### When to Use Parallel Processing
- **Effective when**: Many samples (>10), reasonable bootstrap (100-200)
- **Less effective when**: Few samples (<5), very high bootstrap (>500)
- **Memory consideration**: Each worker needs ~RAM/nCores available

## Migration Guide

### For Existing Code Using Original `p_iNEXT()`

#### Quick Migration (Minimal Changes)
```r
# OLD:
result <- p_iNEXT(otu_mat, q = c(0), nCores = 4)
sample_list <- result$inext_results

# NEW (equivalent):
sample_list <- p_iNEXT(otu_mat, q = c(0), nCores = 4, combine = FALSE)
```

#### Recommended Migration (Use New Features)
```r
# OLD workflow:
result <- p_iNEXT(otu_mat, q = c(0, 1, 2), nCores = 4)
# ... manual data extraction and plotting ...

# NEW workflow:
result <- p_iNEXT(otu_mat, q = c(0, 1, 2), nCores = 4, combine = TRUE)
ggiNEXT(result, type = 1, facet.var = "Order.q")
```

## Complete Example Workflow

```r
library(phyloseq)
library(iNEXT)

# Load phyloseq object
physeq <- main_mxg_physeq_list$mxg_lamps_2022$lamps_2022_AMF_physeq

# Extract OTU table (samples as rows, OTUs as columns)
otu_mat <- as.matrix(t(otu_table(physeq)))

# Calculate endpoint
max_lib_size <- max(rowSums(otu_mat))
endpoint <- max_lib_size * 2

# Run parallel iNEXT with improved version
result <- p_iNEXT(
  x = otu_mat,
  q = c(0, 1, 2),
  endpoint = endpoint,
  nboot = 100,
  nCores = 4,
  combine = TRUE,
  verbose = TRUE
)

# Plot with ggiNEXT (original iNEXT plotting function)
library(ggplot2)

# Size-based rarefaction curves
p1 <- ggiNEXT(result, type = 1, facet.var = "Order.q") +
  theme_bw() +
  labs(title = "Rarefaction/Extrapolation Curves")

# Sample completeness
p2 <- ggiNEXT(result, type = 2) +
  theme_bw() +
  labs(title = "Sample Completeness")

# Coverage-based
p3 <- ggiNEXT(result, type = 3, facet.var = "Order.q") +
  theme_bw() +
  labs(title = "Coverage-based Curves")

print(p1)
print(p2)
print(p3)

# Or use custom plotting
df <- flatten_iNEXT_results(result)
p_custom <- gg_inext_custom(df, type = 1, facet.var = "Order.q", color.var = "Assemblage")
print(p_custom)
```

## Recommendations

### When to Use Original `p_iNEXT.R`
- Legacy code that depends on specific output format
- Already have custom plotting code built around it

### When to Use Improved `p_iNEXT_improved.R`
- **New projects** - Better design and features
- **Want ggiNEXT compatibility** - Set `combine = TRUE`
- **Need better error handling** - More robust with many samples
- **Complex workflows** - More flexible with helper functions
- **Production code** - Better documentation and maintainability

### Future Improvements (Not Implemented)
1. **Deep parallelization**: Parallelize bootstrap within iNEXT
   - Requires modifying iNEXT source code
   - Complex implementation
   - May not provide significant speedup due to overhead
   
2. **Streaming results**: Process samples in batches for memory efficiency
   - Useful for very large datasets (>1000 samples)
   - More complex implementation
   
3. **Progress reporting**: Real-time progress bar
   - Can use `progressr` package with future
   - Adds dependency

4. **Caching**: Save/load intermediate results
   - Useful for expensive computations
   - Requires careful cache key management

## Conclusion

The improved version maintains the same parallelization strategy (sample-level) while adding:
- Better output compatibility with standard iNEXT workflow
- More robust error handling
- Better documentation
- More flexible usage options
- Cleaner code architecture

The key innovation is the `combine` parameter that produces output directly compatible with `ggiNEXT()`, eliminating the need for custom plotting code in most cases.
