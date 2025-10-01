##############################################
# DADA2 pipeline for 16S
#
#
#
# Project: Energy Farm Collab
# Author: Lorien Radmer
# Origin: ~/dada2.qmd
##############################################
##############################################
# Follow tutorial in https://benjjneb.github.io/dada2/tutorial.html

## libraries
## ---------------------------------------------------------------------------------------------------
library(dada2)


## load and define data
## ---------------------------------------------------------------------------------------------------
path <- "/Users/loriensmacbook/Box Sync/Bioinformatics/energy_farm_16S/demultiplexed"
#list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)


## read quality, filter, and trim
## ---------------------------------------------------------------------------------------------------
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = c(150, 150),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
head(out)


## Error rates
## ---------------------------------------------------------------------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)


## Samples
## ---------------------------------------------------------------------------------------------------
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
dadaFs[[1]]


## Merge
## ---------------------------------------------------------------------------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])


## Sequence table
## ---------------------------------------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


## Remove chimeras
## ---------------------------------------------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)


## Final checks
## ---------------------------------------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim"
)
rownames(track) <- sample.names
head(track)


## Taxonomy
## ---------------------------------------------------------------------------------------------------
taxa <- assignTaxonomy(
  seqtab.nochim,
  "/Users/loriensmacbook/Box Sync/Bioinformatics/energy_farm_16S/silva_nr99_v138.1_train_set.fa.gz",
  multithread = TRUE
)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


## Export progress
## ---------------------------------------------------------------------------------------------------
write.csv(seqtab.nochim, file = "seqtab_nochim.csv", row.names = FALSE)
write.csv(taxa, file = "taxatable.csv", row.names = FALSE)

save.image("dada2_done.rdata")

saveRDS(seqtab.nochim, "seqtab.nochim.rds")
saveRDS(taxa, "taxa.rds")
