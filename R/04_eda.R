###########################################################################
# Exploratory Data Analysis for Nested Phyloseq Objects
#
# This script performs EDA on nested phyloseq lists with the structure:
# main_physeq_list$project_list$sequencing_type_physeq
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-29
##########################################################################

source("R/utils/00_setup.R")

#--------------------------------------------------------
# SECTION 1: Basic Exploration of Nested Phyloseq
#--------------------------------------------------------
# Total abundance only

# Energy Farm Collab
# Examine the structure of your phyloseq list
str(main_mxg_physeq_list, max.level = 2)
names(main_mxg_physeq_list)
length(main_mxg_physeq_list)


# Explore phyloseq
purrr::iwalk(main_mxg_physeq_list, function(project_list, project_name) {
  cat("### PROJECT:", project_name, "###\n")
  purrr::iwalk(project_list, function(physeq_obj, region_name) {
    full_name <- paste(project_name, region_name, sep = "_")
    explore_phyloseq_list(physeq_obj, full_name)
  })
})

# Summary tables
physeq_summary <- purrr::imap_dfr(
  main_mxg_physeq_list,
  function(project_list, project_name) {
    purrr::imap_dfr(
      project_list,
      function(physeq_obj, region_name) {
        data.frame(
          #physeq_list = project_name,
          region = gsub(".*_([^_]+)_physeq$", "\\1", region_name),
          n_taxa = ntaxa(physeq_obj),
          n_samples = nsamples(physeq_obj),
          total_reads = sum(sample_sums(physeq_obj)),
          min_reads_per_sample = min(sample_sums(physeq_obj)),
          max_reads_per_sample = max(sample_sums(physeq_obj)),
          mean_reads_per_sample = mean(sample_sums(physeq_obj)),
          median_reads_per_sample = median(sample_sums(physeq_obj))
        )
      },
      .id = "project_name"
    )
  },
  .id = "project_list"
)


physeq_summary

# Details for each phyloseq object in list
nested_summary <- explore_nested_phyloseq(main_mxg_physeq_list)


#--------------------------------------------------------
# Read Count Analysis for Nested Phyloseq
#--------------------------------------------------------

# Get read count data
nested_read_counts <- analyze_read_counts(main_mxg_physeq_list)

# Visualize read counts by project and sequencing type
read_count_plots <- list(
  # Density plot
  density = ggplot(nested_read_counts, aes(x = n_seqs, fill = project)) +
    geom_density(alpha = 0.6, position = "identity") +
    facet_wrap(~project, scales = "free") +
    labs(
      title = "Read Count Distribution by Project",
      x = "Number of Sequences",
      y = "Density"
    ) +
    theme_minimal() +
    theme(legend.position = "none"),

  # Box plot comparison
  boxplot = ggplot(
    nested_read_counts,
    aes(x = project, y = n_seqs, fill = project)
  ) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_y_log10() +
    labs(
      title = "Read Count Comparison Across Projects",
      x = "Project",
      y = "Number of Sequences (log10)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ),

  # Summary line plot
  summary_line = nested_read_counts %>%
    group_by(project) %>%
    arrange(n_seqs) %>%
    mutate(sample_rank = row_number()) %>%
    ggplot(aes(x = sample_rank, y = n_seqs, color = project)) +
    geom_line(linewidth = 1) +
    facet_wrap(~project, scales = "free") +
    labs(
      title = "Read Count Distribution (Ranked)",
      x = "Sample Rank",
      y = "Number of Sequences"
    ) +
    theme_minimal() +
    theme(legend.position = "none"),

  # Goods coverage
  goods_cov = nested_read_counts |>
    ggplot(aes(x = n_seqs, y = goods, color = crop)) +
    geom_point() +
    geom_smooth(se = FALSE, color = "black") +
    facet_wrap(~project)
)

# Display plots
print(read_count_plots)

# Rarafection curves ---------------------------------

# Mapping for rarefaction curves - all datasets
# This errors out with:
# The total size of the 8 globals exported for future expression (‘FUN()’) is 600.41 MiB. This exceeds the maximum allowed size 100.00 MiB per by R option "future.globals.maxSize".

# options(future.globals.maxSize = 1000 * 1024^2)
# rarecurve_results <- main_mxg_physeq_list %>%
#   purrr::imap(function(project_list, project_name) {
#     purrr::imap(project_list, function(physeq_obj, physeq_name) {
#       # Get phyloseq object name for plot title (clean up formatting)
#       plot_title <- gsub("_physeq$", "", physeq_name) %>%
#         gsub("_", " ", .) %>%
#         tools::toTitleCase(.)

#       otu_mat <- physeq_obj %>%
#         otu_table() %>%
#         data.frame() %>%
#         as.matrix() %>%
#         t()

#       # Calculate endpoint
#       max_lib_size <- max(rowSums(otu_mat))
#       endpoint <- max_lib_size * 2

#       # Run parallel iNEXT
#       inext_result <- p_iNEXT(
#         x = otu_mat,
#         q = c(0, 1, 2),
#         endpoint = endpoint,
#         nboot = 100,
#         nCores = 6,
#         combine = TRUE,
#         verbose = TRUE
#       )

#       # Create ggiNEXT plot with phyloseq object name as title
#       inext_plot <- ggiNEXT(
#         inext_result,
#         type = 1,
#         facet.var = "Order.q",
#         color.var = "Assemblage"
#       ) +
#         theme_bw() +
#         labs(title = plot_title)

#       # Return both iNEXT result and plot
#       list(
#         inext_result = inext_result,
#         inext_plot = inext_plot
#       )
#     })
#   })

# Individual plotting of rare curves to avoid error on future.globals.maxSize
p_iNEXt_list <- function(physeq_obj, q = c(0, 1, 2), nCores = 1, type = 1) {
  # Get phyloseq object name for plot title

  plot_title <- as_label(enexpr(physeq_obj)) %>% #
    str_replace_all(".*\\$", "") %>%
    str_remove("_physeq$") %>%
    str_to_upper()

  otu_mat <- physeq_obj %>%
    otu_table() %>%
    data.frame() %>%
    as.matrix() %>%
    t()

  # Calculate endpoint
  max_lib_size <- max(rowSums(otu_mat))
  endpoint <- max_lib_size * 1.25

  # Helper function to make norrow lines
  set_layer_param <- function(plot, i, param, value) {
    if (!is.null(plot$layers[[i]]$aes_params)) {
      plot$layers[[i]]$aes_params[[param]] <- value
    }
    plot
  }

  # Run parallel iNEXT
  inext_result <- p_iNEXT(
    x = otu_mat,
    q = q,
    endpoint = endpoint,
    nboot = 100,
    nCores = nCores,
    combine = TRUE,
    verbose = TRUE
  )

  # Create ggiNEXT plot with phyloseq object name as title
  inext_plot <- ggiNEXT(
    inext_result,
    se = FALSE,
    type = type,
    facet.var = "Order.q",
    color.var = "Assemblage"
  ) +
    theme_bw() +
    labs(title = plot_title, x = "Number of sequences") +
    guides(color = "none", shape = "none", fill = "none")

  # We want narrow lines and no shape
  narrow <- set_layer_param(inext_plot, 1, "size", 0)
  narrow <- set_layer_param(inext_plot, 2, "linewidth", 0.5)

  # assign back into your results list
  inext_plot_narrow <- narrow

  # Return both iNEXT result and plot
  list(
    inext_result = inext_result,
    inext_plot = inext_plot
  )
}

# Individual rarefaction curves for each dataset
ef_16S_iNEXT <- p_iNEXt_list(
  main_mxg_physeq_list$mxg_ef$ef_16S_physeq,
  q = 0,
  nCores = 4,
  type = 1
)

ef_AMF_iNEXT <- p_iNEXt_list(
  main_mxg_physeq_list$mxg_ef$ef_AMF_physeq,
  q = 0,
  nCores = 4,
  type = 1
)

options(future.globals.maxSize = 1000 * 1024^2)
lamps_16S_2018_iNEXT <- p_iNEXt_list(
  main_mxg_physeq_list$mxg_lamps_2018$lamps_2018_16S_physeq,
  q = 0,
  nCores = 4,
  type = 1
)

lamps_ITS_2018_iNEXT <- p_iNEXt_list(
  main_mxg_physeq_list$mxg_lamps_2018$lamps_2018_ITS_physeq,
  q = 0,
  nCores = 4,
  type = 1
)
#TODO 2025-11-18
lamps_16S_2022_iNEXT <- p_iNEXt_list(
  main_mxg_physeq_list$mxg_lamps_2022$lamps_2022_16S_physeq,
  q = 0,
  nCores = 4,
  type = 1
)

lamps_AMF_2022_iNEXT <- p_iNEXt_list(
  main_mxg_physeq_list$mxg_lamps_2022$lamps_2022_AMF_physeq,
  q = 0,
  nCores = 4,
  type = 1
)

rarecurve_results <- list(
  mxg_ef = list(
    ef_16S_iNEXT = ef_16S_iNEXT,
    ef_AMF_iNEXT = ef_AMF_iNEXT
  ),
  mxg_lamps_2018 = list(
    lamps_16S_2018 = lamps_16S_2018_iNEXT,
    lamps_ITS_2018 = lamps_ITS_2018_iNEXT
  ),
  mxg_lamps_2022 = list(
    lamps_16S_2022 = lamps_16S_2022_iNEXT,
    lamps_AMF_2022 = lamps_AMF_2022_iNEXT
  )
)

# ggiNEXT(
#   test$inext_result,
#   type = 1,
#   facet.var = "Order.q",
#   color.var = "Assemblage"
# ) +
#   theme_bw() +
#   labs(title = "VVV", x = "Number of sequences")
#   guides(color = "none", shape = "none", fill = "none")

save(rarecurve_results, file = "data/output/rdata/rarecurves.rda")
load("data/output/rdata/rarecurves.rda")
