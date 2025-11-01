###########################################################################
# Exploratory Data Analysis for Nested Phyloseq Objects
#
# This script performs EDA on nested phyloseq lists with the structure:
# main_relab_physeq_list$project_list$sequencing_type_physeq
#
# Author: Bolívar Aponte Rolón
# Date: 2025-10-29
##########################################################################

source("R/utils/00_setup.R")

#--------------------------------------------------------
# SECTION 1: Basic Exploration of Nested Dataset
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

# Summary table
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


# Summary of nested phyloseq list
explore_nested_phyloseq <- function(nested_list) {
  results <- purrr::imap(nested_list, function(project_list, project_name) {
    cat("\n", rep("=", 60), "\n")
    cat("PROJECT:", toupper(project_name), "\n")
    cat(rep("=", 60), "\n")

    purrr::imap(project_list, function(physeq_obj, seq_type) {
      cat("\n", rep("-", 40), "\n")
      cat("Sequencing Type:", toupper(seq_type), "\n")
      cat(rep("-", 40), "\n")

      # Basic phyloseq summaries
      cat("Basic Summary:\n")
      print(metagMisc::phyloseq_summary(
        physeq_obj,
        more_stats = FALSE,
        long = FALSE
      ))

      cat("\nRead/Sequencing Summary:\n")
      print(microbiome::summarize_phyloseq(physeq_obj))

      # # Taxonomic distribution
      # if ("phylum" %in% colnames(tax_table(physeq_obj))) {
      #   cat("\nPhylum Distribution:\n")
      #   phyla_dist <- phyloseq_ntaxa_by_tax(
      #     physeq_obj,
      #     TaxRank = "phylum",
      #     relative = FALSE,
      #     add_meta_data = FALSE
      #   ) %>%
      #     as.data.frame() %>%
      #     mutate(sum = sum(N.OTU)) %>%
      #     group_by(phylum) %>%
      #     summarise(occurance_in_samples = n()) %>%
      #     arrange(desc(occurance_in_samples))

      #   print(phyla_dist)
      # }

      # Sample and taxa counts
      list(
        project = project_name,
        region = gsub(".*_([^_]+)_physeq$", "\\1", seq_type),
        n_samples = nsamples(physeq_obj),
        n_taxa = ntaxa(physeq_obj),
        sample_vars = sample_variables(physeq_obj),
        rank_names = rank_names(physeq_obj),
        total_reads = sum(sample_sums(physeq_obj)),
        min_reads_per_sample = min(sample_sums(physeq_obj)),
        max_reads_per_sample = max(sample_sums(physeq_obj)),
        mean_reads_per_sample = mean(sample_sums(physeq_obj)),
        median_reads_per_sample = median(sample_sums(physeq_obj))
      )
    })
  })

  return(results)
}


nested_summary <- explore_nested_phyloseq(main_mxg_physeq_list)


#--------------------------------------------------------
# Read Count Analysis for Nested Structure
#--------------------------------------------------------

# Function to analyze read counts across nested phyloseq objects

analyze_read_counts <- function(nested_list) {
  read_summaries <- purrr::imap(
    nested_list,
    function(project_list, project_name) {
      purrr::imap(project_list, function(physeq_obj, seq_type) {
        # Get read counts
        reads <- readcount(physeq_obj) %>%
          as.data.frame() %>%
          rownames_to_column(var = "sample_id") %>%
          rename(n_seqs = ".") %>%
          group_by(sample_id) %>%
          mutate(
            n_singletons = sum(n_seqs == 1),
            goods = 1 - (n_singletons / n_seqs),
            project = gsub("_physeq", "", seq_type),
            region = gsub(".*_([^_]+)_physeq$", "\\1", seq_type)
          )

        metadata <- physeq_obj %>%
          physeq2df() %>%
          select(sample_id, crop)

        new_df <- dplyr::left_join(reads, metadata, by = "sample_id")

        return(new_df)
      })
    }
  ) %>%
    purrr::list_flatten(name_spec = "{outer}_{inner}") %>%
    purrr::map_dfr(~.x)

  return(read_summaries)
}

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

# TODO
# Rarefacttion curves - all datasets
test <- quickRareCurve(
  t(as.matrix(data.frame(otu_table(
    main_mxg_physeq_list$mxg_lamps_2022$lamps_2022_AMF_physeq
  )))),
  max.cores = F,
  nCores = 4,
  label = FALSE
)

test2 <- iNEXT::iNEXT(
  ant,
  datatype = "incidence_freq"
)


iNEXT::ggiNEXT(test2, type = 1) +
  theme_bw(base_size = 18) +
  theme(legend.position = "None")


test3 <- p_iNEXT(
  t(as.matrix(data.frame(otu_table(
    main_mxg_physeq_list$mxg_lamps_2022$lamps_2022_AMF_physeq
  )))),
  datatype = "abundance",
  max.cores = FALSE,
  nCores = 4
)

out <- iNEXT(spider, q = c(0, 1, 2), datatype = "abundance", endpoint = 500)
ggiNEXT(out, type = 1, facet.var = "Assemblage")

# max_lib_size <- main_mxg_physeq_list$mxg_lamps_2022$lamps_2022_AMF_physeq %>%
#   otu_table() %>%
#   data.frame() %>%
#   as.matrix() %>%
#   rowSums() %>%
#   max()

# result <- p_iNEXT(
#   x = t(as.matrix(data.frame(otu_table(
#     main_mxg_physeq_list$mxg_lamps_2022$lamps_2022_AMF_physeq
#   )))),
#   q = c(0),
#   endpoint = max_lib_size * 2,
#   max.cores = F,
#   nCores = 4
# )

# rarecurve_results <- main_mxg_physeq_list %>%
#   purrr::map(function(project_list) {
#     purrr::map(project_list, function(physeq_obj) {
#       physeq_obj %>%
#         otu_table() %>%
#         data.frame() %>%
#         as.matrix() %>%
#         t() %>%
#         quickRareCurve(max.cores = FALSE, nCores = 4, label = FALSE)
#     })
#   })
