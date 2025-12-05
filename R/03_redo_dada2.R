BiocManager::install("dada2", force = TRUE)


set.seed(755) # random number generator for reproducibility
silva.ref <- file.path(
    out_dir,
    "processed/taxonomy/sh_general_release_dynamic_29.11.2022.fasta"
) # CHANGE ME to location on your machine
taxa <- assignTaxonomy(
    seqtab.nochim,
    silva.ref,
    multithread = TRUE,
    minBoot = 50,
    tryRC = TRUE
) # Multithread = FALSE in Windows. TRUE in Mac/Linux.

taxa <- assignTaxonomy(
    seqtab.nochim,
    "/Users/loriensmacbook/Box Sync/Bioinformatics/energy_farm_16S/silva_nr99_v138.1_train_set.fa.gz",
    multithread = TRUE
)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

# Inspecting the taxonomic assignments:

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
str(taxa.print)
dim(taxa.print)
colnames(taxa.print)

write.csv(
    taxa.print,
    file.path(out_dir, "processed/taxonomy/assign_tax_mim2.csv")
)

# Done
