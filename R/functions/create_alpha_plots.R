create_alpha_plot <- function(
  data,
  metric,
  group_var = NULL,
  title = NULL,
  subtitle = NULL,
  y_title_add = NULL
) {
  if (is.null(group_var) || !group_var %in% colnames(data)) {
    # No grouping - single boxplot
    p <- ggplot(data, aes(x = 1, y = .data[[metric]])) +
      geom_boxplot(alpha = 0.7, fill = "steelblue") +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = ifelse(
          is.null(title),
          paste(str_to_title(metric), "Diversity"),
          title
        ),
        subtitle = subtitle,
        x = "",
        y = paste(str_to_title(metric), str_to_title(y_title_add))
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    # With grouping
    p <- ggplot(
      data,
      aes(
        x = .data[[group_var]],
        y = .data[[metric]],
        fill = .data[[group_var]]
      )
    ) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(
        title = ifelse(
          is.null(title),
          paste(str_to_title(metric), "Diversity"),
          title
        ),
        subtitle = subtitle,
        x = str_to_title(gsub("_", " ", group_var)),
        y = paste(str_to_title(metric), str_to_title(y_title_add))
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  }

  return(p)
}
