#' DNA to a condon data.table
#'
#' Convert a DNA coding sequence into a data.table of codons and amino acids.
#'
#' @param cds.seq Character = A DNA coding sequence.
#'
#' @details Assumes that the sequence is in its correct reading frame, that is,
#' the first DNA nucleotide is the first base of the first codon.
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
#' dna2codonDT(X)
#'
#' @export

dna2codonDT <- function(cds.seq){
  require(data.table)
  require(seqinr)
  require(tidyverse)

  codon.list <- seqinr::splitseq(
    seq = cds.seq  %>% as.character %>% s2c
  )

  codon.tab <- lapply(1:length(codon.list), function(i){
    cod <- codon.list[i]
    n <- length(codon.list[1:i])*3
    amino <- seqinr::translate(cod %>% s2c)
    data.table(
      CODON=i,
      NUC=paste(n-2, n-1, n, sep=','),
      AMINO=amino,
      DNA=cod
    )
  }) %>%
    do.call('rbind',.)

  codon.tab %>% return()
}
