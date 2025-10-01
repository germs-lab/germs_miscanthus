##############################################
# DADA2 pipeline for AMF
#
#
#
# Project: Energy Farm Collab
# Author: Lorien Radmer
# Origin: ~/dada2.qmd
##############################################
##############################################

# Follow tutorial in https://benjjneb.github.io/dada2/tutorial.html and https://benjjneb.github.io/dada2/ITS_workflow.html

## libraries
## -------------------------------------------------------------------------------
library(dada2)
library(ShortRead)
library(Biostrings)


## load and define data
## -------------------------------------------------------------------------------
path <- "/Users/loriensmacbook/Desktop/energy_farm_AMF/demultiplexed"
#list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)


## Identify primers
## -------------------------------------------------------------------------------
FWD <- "AATGATACGGCGACCACCGAGATCTACACTATGGTAATTCTTTGGAGGGCAAGTCTGGTGCC" # 5'- 3' Complete forward primer construct (NS31f_il) (AMF)
nchar(FWD) #Number of primer nucleotides.
REV <- "CAAGCAGAAGACGGCATACGAGATNNNNNNNNNNNNAGTCAGTCAGACGAACCCAAACACTTTGGTTTCC" # 5'- 3' Complete reverse primer construct (AML2r-il) (AMF)
nchar(REV)


## Verify orientation of the primers
## -------------------------------------------------------------------------------
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer) # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(
    Forward = dna,
    Complement = Biostrings::complement(dna),
    Reverse = Biostrings::reverse(dna),
    RevComp = Biostrings::reverseComplement(dna)
  )
  return(sapply(orients, toString)) # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


## Pre-filter to remove Ns
## -------------------------------------------------------------------------------
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


## Checking that primers were removed
## -------------------------------------------------------------------------------
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])
)


## Load cutadapt
## -------------------------------------------------------------------------------
cutadapt <- "/Users/loriensmacbook/miniconda3/envs/cutadapt/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R


## Run cutadapt
## -------------------------------------------------------------------------------
path.cut <- file.path(path, "cutadapt")
if (!dir.exists(path.cut)) {
  dir.create(path.cut)
}
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
for (i in seq_along(fnFs)) {
  system2(
    cutadapt,
    args = c(
      R1.flags,
      R2.flags,
      "-n",
      2, # -n 2 required to remove FWD and REV from reads
      "-o",
      fnFs.cut[i],
      "-p",
      fnRs.cut[i], # output files
      fnFs.filtN[i],
      fnRs.filtN[i]
    )
  ) # input files
}


## Check for primers
## -------------------------------------------------------------------------------
rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]])
)


## Formatting
## -------------------------------------------------------------------------------
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(
  path.cut,
  pattern = "_R1_001.fastq",
  full.names = TRUE
))
cutRs <- sort(list.files(
  path.cut,
  pattern = "_R2_001.fastq",
  full.names = TRUE
))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_S")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


## Inspect quality
## -------------------------------------------------------------------------------
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])


## Filter and trim
## -------------------------------------------------------------------------------
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(
  cutFs,
  filtFs,
  cutRs,
  filtRs,
  trimLeft = 10,
  truncLen = c(190, 190),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
) # on windows, set multithread = FALSE
head(out)


## Error rates
## -------------------------------------------------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)


## Sample Inference
## -------------------------------------------------------------------------------
# if some samples are missing data, add next two lines:
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
dadaFs[[1]]


## Merge
## -------------------------------------------------------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
head(mergers[[1]])


## Sequence table - for merged
## -------------------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


## Sequence table - F only
## -------------------------------------------------------------------------------
seqtabf <- makeSequenceTable(dadaFs)
dim(seqtabf)
table(nchar(getSequences(seqtabf)))


## Sequence table - R only
## -------------------------------------------------------------------------------
seqtabr <- makeSequenceTable(dadaRs)
dim(seqtabr)
table(nchar(getSequences(seqtabr)))


## Remove chimeras - merged
## -------------------------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)


## Remove chimeras - F only
## -------------------------------------------------------------------------------
seqtabf.nochim <- removeBimeraDenovo(
  seqtabf,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
dim(seqtabf.nochim)
sum(seqtabf.nochim) / sum(seqtabf)


## Remove chimeras - R only
## -------------------------------------------------------------------------------
seqtabr.nochim <- removeBimeraDenovo(
  seqtabr,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
dim(seqtabr.nochim)
sum(seqtabr.nochim) / sum(seqtabr)


## Final checks - for merged
## -------------------------------------------------------------------------------
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


## Final checks - F only
## -------------------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
trackf <- cbind(out, sapply(dadaFs, getN), rowSums(seqtabf.nochim))
colnames(trackf) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(trackf) <- sample.names
head(trackf)


## Final checks - R only
## -------------------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
trackr <- cbind(out, sapply(dadaRs, getN), rowSums(seqtabr.nochim))
colnames(trackr) <- c("input", "filtered", "denoisedR", "nonchim")
rownames(trackr) <- sample.names
head(trackr)


# seqtab.nochim <- readRDS("/Users/loriensmacbook/Box Sync/Bioinformatics/energy_farm_AMF/seqtab.nochim.rds")

## Taxonomy
## -------------------------------------------------------------------------------
taxa <- assignTaxonomy(
  seqtab.nochim,
  "/Users/loriensmacbook/Desktop/energy_farm_AMF/maarjam_database_SSU_TYPE.PD.2021.fasta",
  multithread = TRUE
)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


## Taxonomy - F only
## -------------------------------------------------------------------------------
taxaf <- assignTaxonomy(
  seqtabf.nochim,
  "/Users/loriensmacbook/Desktop/energy_farm_AMF/maarjam_database_SSU_TYPE.PD.2021.fasta",
  multithread = TRUE
)
taxaf.print <- taxaf
rownames(taxaf.print) <- NULL
head(taxaf.print)


## Taxonomy - R only
## -------------------------------------------------------------------------------
taxar <- assignTaxonomy(
  seqtabr.nochim,
  "/Users/loriensmacbook/Desktop/energy_farm_AMF/maarjam_database_SSU_TYPE.PD.2021.fasta",
  multithread = TRUE
)
taxar.print <- taxar
rownames(taxar.print) <- NULL
head(taxar.print)


## Export progress
## -------------------------------------------------------------------------------
write.csv(seqtab.nochim, file = "seqtab_nochim.csv", row.names = FALSE)
write.csv(taxa, file = "taxatable.csv", row.names = FALSE)

save.image("dada2_done.rdata")

saveRDS(seqtab.nochim, "seqtab.nochim.rds")
saveRDS(taxa, "taxa.rds")


## Export progress - F only
## -------------------------------------------------------------------------------
write.csv(seqtabf.nochim, file = "seqtabf_nochim.csv", row.names = FALSE)
write.csv(taxaf, file = "taxaftable.csv", row.names = FALSE)

save.image("dada2f_done.rdata")

saveRDS(seqtabf.nochim, "seqtabf.nochim.rds")
saveRDS(taxaf, "taxaf.rds")


## Export progress - R only
## -------------------------------------------------------------------------------
write.csv(seqtabr.nochim, file = "seqtabr_nochim.csv", row.names = FALSE)
write.csv(taxar, file = "taxartable.csv", row.names = FALSE)

save.image("dada2r_done.rdata")

saveRDS(seqtabr.nochim, "seqtabr.nochim.rds")
saveRDS(taxar, "taxar.rds")
