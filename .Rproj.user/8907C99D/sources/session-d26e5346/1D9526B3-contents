#' Convert a data.table of gene annotations to a list of exon sequences
#'
#' A data.table with information on the start and end of exons is
#' used to generate a list containing each exon as an indexed item.
#' The indexes exon sequences are in the order and orientation that
#' they would be transcribed in. That is, they can be pasted together
#' and translated to form the protein sequence. Note, this function has
#' been designed to work for a single gene on a single chromosome.
#' If you want to run multiple genes, you will need to loop the function.
#'
#' @param annotDT Data.table: Contains information on the exon positions.
#' Requires the following columns:
#' \enumerate{
#'    \item The chromosome ID (see param \code{chromCol}).
#'    \item The start position of the exon on the positive strand, that is,
#'    left to right on the genomic sequence (see param \code{startCol}).
#'    \item The end position of the exon on the positive strand, that is,
#'    left to right on the genomic sequence (see param \code{endCol}).
#'    \item The strand position of the exon, '+' for positive, and '-' for
#'    negative (see param \code{strandCol}).
#' }
#'
#' @param genomeSeq DNAStringSet: The loaded genome sequence as a
#' \code{DNAStringSet} object, as per R's \code{Biostrings} package
#' (Pagès et al.). The sequence names must match the chromosome name in \code{annotDT}.
#'
#' @param chromCol Character: The chromosome column name in \code{annotDT}.
#' Default is \code{'CHROM'}.
#'
#' @param startCol Character: The exon start position in \code{annotDT}.
#' Default is \code{'START'}.
#'
#' @param endCol Character: The exon end position in \code{annotDT}.
#' Default is \code{'END'}.
#'
#' @param strandCol Character: The exon strand position in \code{annotDT}.
#' Default is \code{'STRAND'}.
#'
#' @returns Returns a list, with each indexed item a character vector,
#' the exon extracted from the genome. Exons are ordered based on their order
#' of transcriptions, and are in the correct orientation.
#'
#' @export
annotDT2exonList <- function(annotDT, genomeSeq, chromCol='CHROM', startCol='START', endCol='END', strandCol='STRAND'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  # Required packages
  require(Biostrings); require(data.table); require(tidyverse)

  # Check columns
  annotDT <- as.data.table(annotDT)
  col_check <- sum(c(chromCol,startCol,endCol,strandCol) %in% colnames(annotDT))
  if(col_check!=4){
    stop('Columns specified in arguments `chromCol`, `startCol`, `endCol`, and `strandCol` must all be in `annotDT`. See ?annotDT2exonList.')
  }
  setnames(annotDT, c(chromCol,startCol,endCol,strandCol), c('CHROM','START','END','STRAND'))

  # Check that the genome is a DNAStringSet class.
  if(class(genomeSeq)!='DNAStringSet'){
    stop('Argument `genomeSeq` must be a DNAStringSet class object. See ?annotDT2exonList.')
  }

  # --------------------------------------------+
  # Execute
  # --------------------------------------------+
  strand <- annotDT$STRAND[1]

  chrom <- annotDT$CHROM[1]

  genome_sub <- genomeSeq[chrom]

  # Iterate over exons
  if(strand=='+'){
    result <- setorder(annotDT, START) %>%
      .[, EXON:=1:.N] %>%
      split(., by='EXON') %>%
      lapply(., function(X){
        Biostrings::subseq(genome_sub,X$START,X$END) %>%
          as.character()
      })
  } else if(strand=='-'){
    result <- setorder(annotDT, -START) %>%
      .[, EXON:=1:.N] %>%
      split(., by='EXON') %>%
      lapply(., function(X){
        Biostrings::reverseComplement(Biostrings::subseq(genome_sub,X$START,X$END)) %>%
          as.character()
      })
  }

  # Return a list of coding sequences
  return(result)
}
