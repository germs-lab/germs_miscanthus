# Quick Start Guide - Parallel iNEXT

## TL;DR - Just Want to Run iNEXT in Parallel?

```r
# 1. Source the improved version
source("R/functions/p_iNEXT_improved.R")

# 2. Prepare your data (samples as rows, species as columns)
otu_mat <- as.matrix(t(otu_table(your_phyloseq)))

# 3. Run parallel iNEXT
result <- p_iNEXT(
  x = otu_mat,
  q = c(0, 1, 2),          # Hill numbers
  nCores = 4,              # Use 4 cores
  combine = TRUE           # Merge results (default)
)

# 4. Plot with standard ggiNEXT
library(iNEXT)
ggiNEXT(result, type = 1, facet.var = "Order.q")
```

**That's it!** ✅

## What You Get

```
Input: OTU Matrix (samples × species)
  ↓
[p_iNEXT with 4 cores]
  ↓
Output: iNEXT object (compatible with ggiNEXT)
  ↓
[ggiNEXT]
  ↓
Beautiful rarefaction plots! 📊
```

## Files You Need

| File | What It Does | Do You Need It? |
|------|--------------|-----------------|
| `p_iNEXT_improved.R` | **Main implementation** | ✅ YES - This is what you use |
| `p_iNEXT_README.md` | Complete user guide | 📖 Read for detailed workflows |
| `p_iNEXT_examples.R` | Working examples | 🧪 Run to test |
| `p_iNEXT_comparison.md` | Technical comparison | 📊 If migrating from original |
| `IMPLEMENTATION_SUMMARY.md` | Technical deep dive | 🔬 If interested in design |
| `p_iNEXT_deep_parallel.R` | Experimental version | ⚠️ Probably not |
| `p_iNEXT.R` (original) | Legacy version | 🗂️ For backward compatibility |

## Common Use Cases

### Use Case 1: Standard Analysis
```r
source("R/functions/p_iNEXT_improved.R")
result <- p_iNEXT(otu_matrix, q = c(0, 1, 2), nCores = 4)
ggiNEXT(result, type = 1, facet.var = "Order.q")
```

### Use Case 2: From Phyloseq
```r
library(phyloseq)
source("R/functions/p_iNEXT_improved.R")

otu_mat <- as.matrix(t(otu_table(physeq)))
result <- p_iNEXT(otu_mat, q = c(0, 1, 2), nCores = 4)
ggiNEXT(result, type = 1)
```

### Use Case 3: Custom Plotting
```r
source("R/functions/p_iNEXT_improved.R")
result <- p_iNEXT(otu_matrix, nCores = 4)
df <- flatten_iNEXT_results(result)

# Now use ggplot2 for custom plots
library(ggplot2)
ggplot(df, aes(x = size_based.m, y = size_based.qD, color = sample)) +
  geom_line()
```

### Use Case 4: Individual Sample Access
```r
source("R/functions/p_iNEXT_improved.R")
samples <- p_iNEXT(otu_matrix, nCores = 4, combine = FALSE)

# Access individual samples
sample1 <- samples$Sample1
ggiNEXT(sample1, type = 1)
```

## Key Parameters

| Parameter | Default | What It Does |
|-----------|---------|--------------|
| `x` | required | Your OTU matrix (samples × species) |
| `q` | `c(0,1,2)` | Hill numbers (0=richness, 1=Shannon, 2=Simpson) |
| `nCores` | 1 | Number of cores to use |
| `combine` | `TRUE` | Merge into single iNEXT object? |
| `nboot` | 100 | Bootstrap replicates |
| `endpoint` | 2x max | Extrapolation endpoint |

## How Many Cores Should I Use?

```r
# Check available cores
parallelly::availableCores()

# Good rule: use N-1 (leave one for system)
nCores <- parallelly::availableCores() - 1

# For small datasets (<10 samples), use 1
if (nrow(otu_mat) < 10) nCores <- 1
```

## Troubleshooting

### "Out of memory"
→ Reduce cores: `nCores = 2`
→ Or reduce bootstrap: `nboot = 50`

### "Taking too long"
→ Increase cores: `nCores = 8`
→ Or reduce bootstrap: `nboot = 50`

### "Not compatible with ggiNEXT"
→ Make sure `combine = TRUE`

### "Need individual samples"
→ Set `combine = FALSE`

## Need More Help?

1. **Examples**: Run `source("R/functions/p_iNEXT_examples.R")`
2. **User Guide**: Read `p_iNEXT_README.md`
3. **Technical Details**: Read `IMPLEMENTATION_SUMMARY.md`
4. **Original iNEXT**: See https://github.com/AnneChao/iNEXT

## Comparison with Original p_iNEXT.R

| Feature | Original | Improved |
|---------|----------|----------|
| Parallelization | ✅ Yes | ✅ Yes |
| ggiNEXT compatible | ❌ No | ✅ Yes |
| Error handling | ⚠️ Basic | ✅ Robust |
| Documentation | ⚠️ Basic | ✅ Comprehensive |
| Helper functions | ❌ No | ✅ Yes |
| Output flexibility | ⚠️ Limited | ✅ High |

**Recommendation**: Use improved version for all new projects.

## What's Next?

After you get results:
1. ✅ Plot with `ggiNEXT(result, type = 1)`
2. ✅ Try different plot types (type = 1, 2, 3)
3. ✅ Try different faceting (facet.var = "Order.q", "Assemblage", "Both")
4. ✅ Customize with ggplot2 themes
5. ✅ Save plots with `ggsave()`

## Full Example Workflow

```r
# Load packages
library(phyloseq)
library(iNEXT)
library(ggplot2)

# Source improved iNEXT
source("R/functions/p_iNEXT_improved.R")

# Load your phyloseq object
physeq <- readRDS("path/to/phyloseq.rds")

# Extract OTU table (transpose so samples are rows)
otu_mat <- as.matrix(t(otu_table(physeq)))

# Run parallel iNEXT
result <- p_iNEXT(
  x = otu_mat,
  q = c(0, 1, 2),
  nCores = 4,
  nboot = 100,
  verbose = TRUE
)

# Plot 1: Size-based rarefaction
p1 <- ggiNEXT(result, type = 1, facet.var = "Order.q") +
  theme_bw() +
  labs(title = "Rarefaction Curves")

# Plot 2: Sample completeness
p2 <- ggiNEXT(result, type = 2) +
  theme_bw() +
  labs(title = "Sample Completeness")

# Plot 3: Coverage-based
p3 <- ggiNEXT(result, type = 3, facet.var = "Order.q") +
  theme_bw() +
  labs(title = "Coverage-based Rarefaction")

# Save plots
ggsave("rarefaction.pdf", p1, width = 10, height = 6)
ggsave("completeness.pdf", p2, width = 10, height = 6)
ggsave("coverage.pdf", p3, width = 10, height = 6)

# Get asymptotic estimates
result$AsyEst

# Get data info
result$DataInfo
```

## Questions?

- **"How does this compare to original iNEXT?"** → Same computation, parallel execution, better output format
- **"Is it faster?"** → Yes, 2.5-3.5x speedup with 4 cores (typical)
- **"Can I still use ggiNEXT?"** → Yes! That's the main improvement
- **"Do I need to change my code?"** → Only if you want ggiNEXT compatibility (recommended)
- **"What about deep parallelization?"** → Usually not worth the overhead, but see `p_iNEXT_deep_parallel.R` if interested

---

**Bottom Line**: Use `p_iNEXT_improved.R` with `combine = TRUE` for hassle-free parallel iNEXT that works with standard ggiNEXT plotting. 🚀
