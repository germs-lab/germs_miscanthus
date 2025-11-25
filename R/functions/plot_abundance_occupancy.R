plot_abundance_occupancy <- function(
  df,
  threshold,
  title = NULL,
  xtitle = NULL,
  ytitle = NULL
) {
  # Set defaults for any NULL parameters
  if (is.null(title)) {
    title <- "Abundance-Occupancy curve"
  }
  if (is.null(xtitle)) {
    xtitle <- "Log10(Mean Relative Abundance)"
  }
  if (is.null(ytitle)) {
    ytitle <- "Occupancy"
  }

  title_labs <- list(labs(
    title = title,
    x = xtitle,
    y = ytitle
  ))

  plot <- ggplot(
    df,
    aes(
      x = base::log10(mean_rel_abundance),
      y = occupancy_prop,
      fill = base::ifelse(occupancy_prop > threshold, "Core", "Not core")
    )
  ) +
    geom_point(
      shape = 21,
      size = 1.5,
      alpha = 0.9,
      aes(color = after_scale(fill))
    ) +
    geom_hline(yintercept = threshold, linetype = 2, linewidth = 1) +
    scale_fill_manual(
      values = c("Core" = "#CC2D35", "Not core" = "grey")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "bold"
      ),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "inside",
      legend.position.inside = c(0.92, 0.1),
      legend.background = element_rect(
        fill = ggplot2::alpha("white", 0.7),
        color = NA
      ),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.key.width = grid::unit(0.4, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    ) +
    title_labs

  return(plot)
}
