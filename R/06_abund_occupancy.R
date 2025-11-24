###########################################################################
# Abundance Occupancy Analysis on Miscanthus
#
#
# Author: Bolívar Aponte Rolón
# Date: 2025-11-24
##########################################################################

source("R/utils/00_setup.R")


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
