# From https://raw.githubusercontent.com/microsud/microbiomeutilities/046a9f928f35fc2a02dd6cfb9175a0c7ddf65c4b/R/extensions.R

add_refseq <- function(x, tag = "ASV") {
  if (!is(x, "phyloseq")) {
    stop("Input is not an object of phyloseq class")
  }

  nucl <- Biostrings::DNAStringSet(taxa_names(x))
  names(nucl) <- taxa_names(x)
  x <- merge_phyloseq(x, nucl)

  rm(nucl)

  if (is.na(tag) || is.null(tag)) {
    taxa_names(x) <- paste0("taxa", seq(ntaxa(x)))
    return(x)
  } else {
    taxa_names(x) <- paste0(tag, seq(ntaxa(x)))
    return(x)
  }
}
