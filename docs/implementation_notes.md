# Implementation Notes: Shared ASV Analysis

## Design Decisions

### 1. Data Structure (07_abund_occupancy.R)

The implementation creates a two-level nested list structure:

```r
core_asvs_by_threshold
├── threshold_0.6
│   ├── ef_16S
│   ├── lamps_2018_16S
│   └── lamps_2022_16S
├── threshold_0.7
│   └── ...
└── ...

shared_asvs_by_threshold
├── threshold_0.6
│   ├── all_projects       # ASVs in ALL projects
│   ├── two_or_more        # ASVs in 2+ projects
│   ├── by_project         # ASVs per project
│   └── asv_project_counts # Count vector
└── ...
```

**Rationale:**
- Easy to filter and subset for different analyses
- Compatible with downstream UpSet plot functions
- Preserves full information for multiple visualization types
- Efficient for programmatic access

### 2. Threshold Selection

Default thresholds: 60%, 70%, 80%, 90%

**Rationale:**
- 60%: Captures moderately prevalent ASVs (standard "core microbiome" definition)
- 70-80%: Progressively more stringent core definitions
- 90%: Very strict, near-ubiquitous ASVs
- Excludes lower thresholds (<60%) to focus on truly prevalent ASVs
- Excludes 100% as it's too stringent for real-world data

### 3. Project Filtering

Uses `grep("_DNA$", ...)` to select only DNA-based 16S projects.

**Rationale:**
- Excludes RNA-based projects which have different ASV profiles
- Ensures comparability between projects
- Follows repository's existing naming convention
- Easy to modify for different project subsets

### 4. Visualization Strategy (08_upset_plots.R)

Uses ComplexUpset instead of Venn diagrams.

**Rationale:**
- Venn diagrams become unreadable with 3+ sets
- UpSet plots scale better for multiple comparisons
- Shows exact intersection sizes clearly
- Consistent with existing repository visualizations
- Publication-quality output

### 5. Table Generation (09_asv_tables.R)

Creates three types of tables:
1. All-projects shared ASVs (per threshold)
2. 2+ projects shared ASVs (per threshold)
3. Summary across all thresholds

**Rationale:**
- "All projects" focuses on truly core microbiome
- "2+ projects" captures broader sharing patterns
- Summary table provides quick overview
- HTML format allows interactive viewing
- gt() tables are publication-ready

### 6. Taxonomic Annotation

Uses the first phyloseq object as taxonomy reference.

**Rationale:**
- All phyloseq objects have consistent taxonomy (reassigned together in 05_rebuild_and_transform.R)
- Avoids redundant taxonomy merging
- Fast lookup for large ASV sets
- Maintains consistency with existing workflow

### 7. File Naming Convention

Pattern: `shared_asvs_{threshold}_{type}.{ext}`

**Examples:**
- `shared_asvs_threshold_0.6_upset.svg`
- `shared_asvs_2plus_threshold_0.7_taxonomy_table.html`

**Rationale:**
- Clear and descriptive
- Sortable by threshold
- Easy to find specific outputs
- Consistent with repository conventions

## Key Features

### 1. Robustness

- Handles cases with no shared ASVs gracefully
- Validates input data existence
- Provides informative console output
- Clear error messages

### 2. Reproducibility

- All outputs are deterministic
- Results depend only on input phyloseq objects
- No random components
- Clear data provenance

### 3. Extensibility

- Easy to add new thresholds
- Simple to modify project selection
- Modular design for new analysis types
- Well-commented for future modifications

### 4. Performance

- Efficient use of purrr for iteration
- Minimal data copying
- Single-pass calculations where possible
- Reasonable for datasets with 1000s of ASVs

## Integration with Existing Workflow

### Dependencies

Uses existing repository functions:
- `summarize_abundance_occupancy()` - From R/functions/
- `create_binary_df_from_flat()` - From R/functions/upset_plot_helpers.R
- `create_asv_upset_data()` - From R/functions/

### Consistency

Follows repository patterns:
- purrr functional programming style
- tidyverse data manipulation
- Consistent naming conventions
- Standardized section headers

### Data Flow

```
main_16S_physeq_list (from 00_setup.R)
    ↓
07_abund_occupancy.R
    ↓
shared_asvs_by_threshold.rda
    ↓
├── 08_upset_plots.R → SVG files
└── 09_asv_tables.R → HTML tables
```

## Testing Recommendations

When running the scripts for the first time:

1. **Check phyloseq objects:**
   ```r
   names(main_16S_physeq_list)
   sapply(main_16S_physeq_list, ntaxa)
   ```

2. **Verify outputs:**
   ```r
   load("data/output/rdata/shared_asvs_by_threshold.rda")
   names(shared_asvs_by_threshold)
   sapply(shared_asvs_by_threshold, function(x) length(x$all_projects))
   ```

3. **Check plot generation:**
   - Look for SVG files in `data/output/figures/`
   - Verify they open correctly

4. **Review tables:**
   - Open HTML files in browser
   - Check taxonomic information is complete
   - Verify counts make sense

## Potential Extensions

Future enhancements could include:

1. **Interactive tables** using DT package
2. **Phylum-level analysis** in addition to ASV-level
3. **Statistical testing** for enrichment in shared ASVs
4. **Network visualization** of ASV co-occurrence
5. **Export to Excel** for non-R users
6. **Rarefaction curves** for shared ASVs

## Author

Bolívar Aponte Rolón  
Implementation Date: 2026-01-26
