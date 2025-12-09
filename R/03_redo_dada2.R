###################################################################
# Reassignment of taxonomy via dada2::assignTaxonomy()
#
# After the union of all otu tables in 02_main_subset_fasta.R we reassing taxonomy to the combined otu table for consitency in database version and ASV names

# Author: Bolívar Aponte Rolón
# Date: 2025-12-09
###################################################################

library("Rcpp")
library("dada2") #v1.38.0

load("data/output/rdata/new_16S_DNS_seqtab.rda")

set.seed(755) # random number generator for reproducibility
silva.ref <- "data/input/databases/SILVA/silva_nr99_v138.2_toSpecies_trainset.fa" # CHANGE ME to location on your machine

cat("\n", rep("=", 40), "\n")
cat("Taxonomy assignment\n")
cat(rep("=", 40), "\n")

taxa <- assignTaxonomy(
    new_16S_DNA_seqtab,
    silva.ref,
    multithread = TRUE,
    minBoot = 50,
    tryRC = TRUE
) # Multithread = FALSE in Windows. TRUE in Mac/Linux.


cat("\n", rep("=", 40), "\n")
cat("Inspection\n")
cat(rep("=", 40), "\n")

# Inspecting the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
str(taxa.print)
dim(taxa.print)
colnames(taxa.print)


cat("\n", rep("=", 40), "\n")
cat("Saving\n")
cat(rep("=", 40), "\n")

write.csv(
    taxa.print,
    "data/output/taxonomy/new_16S_DNA_tax.csv"
)

cat("\n", rep("=", 40), "\n")
cat("Complete\n")
cat(rep("=", 40), "\n")
