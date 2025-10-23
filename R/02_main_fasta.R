## Cleaning ASV names for FASTA file
rownames(taxa_16S) <- paste0("ASV_", 1:nrow(taxa_16S))

asv_fasta <- seqtab2fasta(seqtab_nochim_16S)

seqtab_nochim_16S <- t(seqtab_nochim_16S) # Retaining sequences and asigning shorthand ASV names

row.names(seqtab_nochim_16S) <- sub(">", "", asv_fasta$asv_headers)

asv_16S <- otu_table(seqtab_nochim_16S, taxa_are_rows = TRUE)


## Save FASTA file
dir.create("data/output/processed/sequences/", recursive = TRUE)
write(
  asv_fasta$asv_fasta,
  file.path("data/output/processed/sequences/energy_farm_collab_16S.fa")
)

## Cleaning ASV names for FASTA file

asv_fasta <- seqtab2fasta(seqtab_nochim_AMF)

seqtab_nochim_AMF <- t(seqtab_nochim_AMF) # Retaining sequences and asigning shorthand ASV names

row.names(seqtab_nochim_AMF) <- sub(">", "", asv_fasta$asv_headers)

asv_AMF <- otu_table(seqtab_nochim_AMF, taxa_are_rows = TRUE)


## Save FASTA file
dir.create("data/output/processed/sequences/", recursive = TRUE)
write(
  asv_fasta$asv_fasta,
  file.path("data/output/processed/sequences/energy_farm_collab_AMF.fa")
)
