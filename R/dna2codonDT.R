#' DNA to a condon data.table
#'
#' Convert a DNA coding sequence into a data.table of codons and amino acids.
#'
#' @param cdsSeq Character = A DNA coding sequence.
#'
#' @param compressTab Logical = Should the table be compressed? Default is TRUE,
#' in which case, each row is a unique codon. If FALSE, then each codon is represented
#' by 3 rows, one for each nucleotide comprising the codon.
#'
#' @param geneticCode Integer = A value relating to the \code{numcode} argument in
#' \code{seqinr::translate}.
#'
#' @details Assumes that the sequence is in its correct reading frame, that is,
#' the first DNA nucleotide is the first base of the first codon. Default is 1.
#'
#' @return
#' Returns a data.table with the following columns:
#'
#' \enumerate{
#'    \item \code{$CODON} = The codon number, 1:N.
#'    \item \code{$NUC} = The nucleotides comprising the codon.
#'    \item \code{$AMINO} = The amino acid residue.
#'    \item \code{$DNA} = The DNA bases.
#' }
#'
#' @examples
#' X <- 'ATGCGTACTTCA'
#'
#' dna2codonDT(X, compressTab=TRUE)
#'
#' dna2codonDT(X, compressTab=FALSE)
#'
#' @export

dna2codonDT <- function(cdsSeq, compressTab=FALSE, geneticCode=1){
  require(data.table)
  require(seqinr)
  require(tidyverse)

  # Split the sequence into codons
  codonList <- seqinr::splitseq(
    seq = cdsSeq  %>% as.character %>% s2c
  )

  # Iterate over ith the codons
  codonTab <- lapply(1:length(codonList), function(i){
    # Subset codon
    cod <- codonList[i]

    # Get the relative nucleotide positions within the codon
    n <- length(codonList[1:i])*3

    # Get the amino acid
    amino <- seqinr::translate(cod %>% s2c, numcode=geneticCode)

    # Make the table for this codon
    if(compressTab==TRUE){
      # Compressed table format
      tab <- data.table(
        CODON=i,
        NUC=paste(n-2, n-1, n, sep='|'),
        AMINO=amino,
        DNA=cod
      )
    } else if(compressTab==FALSE){
      # Uncompressed table format
      tab <- data.table(
        CODON=i,
        NUC=c(n-2, n-1, n),
        AMINO=amino,
        DNA=cod
      )
    }
  }) %>%
    do.call('rbind',.)

  # Output
  codonTab %>% return()
}
