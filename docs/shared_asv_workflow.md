# Shared ASV Analysis Workflow

This documentation describes the workflow for identifying and analyzing shared ASVs (Amplicon Sequence Variants) between projects at different occupancy thresholds.

## Overview

The analysis consists of three main scripts that work together:
1. `07_abund_occupancy.R` - Identifies core ASVs at different thresholds
2. `08_upset_plots.R` - Creates UpSet plots to visualize ASV intersections
3. `09_asv_tables.R` - Generates taxonomic tables for shared ASVs

## Workflow

### Step 1: Identify Core ASVs (07_abund_occupancy.R)

This script:
- Analyzes ASV occupancy across projects at multiple thresholds (60%, 70%, 80%, 90%)
- Identifies ASVs shared between all projects
- Identifies ASVs shared by 2 or more projects
- Saves results to `data/output/rdata/shared_asvs_by_threshold_07.rda`

**Output:**
- `shared_asvs_by_threshold.rda` - R data object containing:
  - `core_asvs_by_threshold` - ASVs for each project at each threshold
  - `shared_asvs_by_threshold` - Intersection data at each threshold

**Data Structure:**
```r
shared_asvs_by_threshold$threshold_0.6$all_projects  # ASVs shared by all projects
shared_asvs_by_threshold$threshold_0.6$two_or_more  # ASVs shared by 2+ projects
shared_asvs_by_threshold$threshold_0.6$core_by_project   # ASVs for each project
shared_asvs_by_threshold$threshold_0.6$asv_project_counts  # Count vector
```

### Step 2: Create UpSet Plots (08_upset_plots.R)

This script:
- Loads the shared ASV data from Step 1
- Creates UpSet plots for each threshold showing ASV intersections
- Annotates plots with taxonomic information
- Saves plots to `data/output/figures/`

**Output:**
- `shared_asvs_threshold_0.6_upset.svg` (and similar for other thresholds)

**Features:**
- Shows intersection sizes between projects
- Highlights ASVs shared by all projects in green
- Includes taxonomic annotation
- High-resolution SVG format for publication

### Step 3: Generate Taxonomic Tables (09_asv_tables.R)

This script:
- Loads the shared ASV data from Step 1
- Extracts taxonomic information from phyloseq objects
- Creates formatted gt() tables with:
  - ASVs shared across all projects
  - ASVs shared by 2+ projects
  - Summary statistics across thresholds
- Saves tables as HTML files

**Output:**
- `shared_asvs_threshold_0.6_taxonomy_table.html` - All projects table
- `shared_asvs_2plus_threshold_0.6_taxonomy_table.html` - 2+ projects table
- `shared_asvs_summary_table.html` - Summary across all thresholds

**Table Features:**
- Full taxonomic classification (Kingdom → Genus)
- Number of projects sharing each ASV
- Color-coded project counts
- Publication-ready formatting

## Usage

Run the scripts in order:

```r
# Step 1: Identify core ASVs at different thresholds
source("R/07_abund_occupancy.R")

# Step 2: Create UpSet plots
source("R/08_upset_plots.R")

# Step 3: Generate taxonomic tables
source("R/09_asv_tables.R")
```

## Dependencies

Required R packages:
- phyloseq
- tidyverse
- vegan
- ComplexUpset
- gt
- gtExtras

These are loaded via `R/utils/00_setup.R`.

## Customization

### Modifying Thresholds

To analyze different occupancy thresholds, edit line 122 in `07_abund_occupancy.R`:

```r
thresholds <- c(0.6, 0.7, 0.8, 0.9)  # Modify as needed
```

### Filtering Projects

To analyze only specific projects, edit line 125 in `07_abund_occupancy.R`:

```r
# Current: analyzes all _DNA projects
project_names <- grep("_DNA$", names(main_16S_physeq_list), value = TRUE)

# Example: analyze only specific projects
project_names <- c("ef_16S_DNA", "lamps_2018_16S_DNA")
```

### Plot Appearance

Modify plot parameters in `08_upset_plots.R` (lines 359-386):
- `n_intersections` - Number of intersections to show
- `width`, `height` - Plot dimensions
- Color schemes - Modify `fill` and `color` parameters

### Table Formatting

Modify table appearance in `09_asv_tables.R`:
- Line 84-112: Modify gt() styling for "all projects" tables
- Line 156-184: Modify gt() styling for "2+ projects" tables
- Line 204-227: Modify summary table formatting

## Output Files

All outputs are saved to `data/output/`:

**R Data Objects:**
- `rdata/shared_asvs_by_threshold.rda`

**Figures:**
- `figures/shared_asvs_threshold_0.6_upset.svg`
- `figures/shared_asvs_threshold_0.7_upset.svg`
- `figures/shared_asvs_threshold_0.8_upset.svg`
- `figures/shared_asvs_threshold_0.9_upset.svg`

**Tables:**
- `figures/shared_asvs_threshold_0.6_taxonomy_table.html`
- `figures/shared_asvs_2plus_threshold_0.6_taxonomy_table.html`
- (and similar for other thresholds)
- `figures/shared_asvs_summary_table.html`

## Troubleshooting

### Error: "object 'main_16S_physeq_list' not found"

Ensure that the phyloseq objects are loaded. Run:
```r
source("R/utils/00_setup.R")
```

### Error: "package 'gt' not found"

Install required packages:
```r
install.packages("gt")
install.packages("gtExtras")
```

### No ASVs shared across all projects

This is expected if the occupancy threshold is too high. Try:
1. Lower the threshold values
2. Check that projects have sufficient overlap
3. Review the "2+ projects" results instead

## Notes

- The workflow focuses on DNA-based 16S projects (identified by "_DNA" suffix)
- Occupancy thresholds represent the proportion of samples where an ASV must be present
- Results are reproducible as long as input phyloseq objects remain unchanged
- All plots and tables are publication-ready with high resolution (600 DPI)

## Author

Bolívar Aponte Rolón  
Date: 2026-01-26
