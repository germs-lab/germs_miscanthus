##########################################################################

# Author: Bolívar Aponte Rolón
# Date: 2025-12-11
##########################################################################
library(UpSetR)
movies <- read.csv(
  system.file("extdata", "movies.csv", package = "UpSetR"),
  header = T,
  sep = ";"
)

listInput <- list(
  EF = rownames(otu_table(main_16S_physeq_list$ef_16S_DNA)),
  LAMPS_2018 = rownames(otu_table(main_16S_physeq_list$lamps_2018_16S_DNA))
)


upset(fromList(listInput), order.by = "freq")

# TODO
# Focus on main_16S_physeq_list then subset to only the MXG within the 16S_DNA
