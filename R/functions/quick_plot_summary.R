# Function for basic plots
quick_plot_summary <- function(physeq_obj, obj_name) {
  # Sample sums barplot
  p1 <- plot_bar(physeq_obj, fill = "Sample") +
    ggtitle(paste("Sample Sums -", obj_name)) +
    theme(legend.position = "none")

  # Taxa abundance heatmap (top 20)
  top_taxa <- names(sort(taxa_sums(physeq_obj), decreasing = TRUE))[1:20]
  physeq_top <- prune_taxa(top_taxa, physeq_obj)

  p2 <- plot_heatmap(physeq_top) +
    ggtitle(paste("Top 20 Taxa -", obj_name))

  return(list(sample_plot = p1, taxa_plot = p2))
}
