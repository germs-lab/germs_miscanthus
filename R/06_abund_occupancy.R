###########################################################################
# Abundance Occupancy Analysis on Miscanthus
#
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-24
##########################################################################

source("R/utils/00_setup.R")

# Getting our cores
main_core_list <- purrr::imap(
  main_physeq_list,
  function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      get_prevalent_rare(
        physeq_obj,
        thresholds = c(90, 80, 70, 60, 30),
        detection = 0 / 100,
        include.lowest = FALSE
      )
    })
  }
)

main_mxg_core_list <- purrr::imap(
  main_mxg_physeq_list,
  function(project_list, project_name) {
    purrr::imap(project_list, function(physeq_obj, physeq_name) {
      get_prevalent_rare(
        physeq_obj,
        thresholds = c(90, 80, 70, 60, 30),
        detection = 0 / 100,
        include.lowest = FALSE
      )
    })
  }
)

# TODO
# Polish in to a cohesive workflow for all 16S phyloseqs
# 2025-11-24
test <- main_physeq_list$ef_physeq_list$ef_16S_physeq %>%
  otu_table() %>%
  as.data.frame() %>%
  summarise(
    abundance = rowSums(.),
    occupancy = 1 * rowSums(. > 0),
    occupancy_prop = occupancy / ncol(.), #aka Prevalence
    mean_abundance = abundance / occupancy,
    mean_rel_abundance = apply(
      vegan::decostand(., method = "total", MARGIN = 2),
      1,
      mean
    )
    # membership = base::ifelse(occupancy_prop > 0.90, "Core", "Not core")
  )
rownames(test) <- rownames(otu_table(
  main_physeq_list$ef_physeq_list$ef_16S_physeq
))

# Test plot
ggplot(
  test,
  aes(
    x = base::log10(mean_rel_abundance),
    y = occupancy_prop,
    fill = base::ifelse(occupancy_prop > 0.60, "Core", "Not core")
  )
) +
  geom_point(
    shape = 21,
    size = 1.5,
    alpha = 0.9,
    aes(color = after_scale(fill))
  ) +
  scale_fill_manual(
    values = c("Core" = "#CC2D35", "Not core" = "grey")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 12,
      face = "bold"
    ),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(
      fill = ggplot2::alpha("white", 0.7),
      color = NA
    ),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.4, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  ) +
  labs(
    title = "Abundance-Occupancy curve",
    x = "Log10(mean abundance)",
    y = "Occupancy"
  )
