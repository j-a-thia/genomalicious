#' DNA to a translated codon data.table
#'
#' Convert a DNA coding sequence into a data.table of codons and amino acids.
#'
#' @param dnaSeq Character: For a coding sequence (e.g., cDNA transcript), a single
#' character vector. For a group of exons, a list of character vectors.
#' See argument \code{type} and Details.
#'
#' @param type Character: one of \code{'cds'} (coding sequence) or \code{'exons'}
#' (exon regions). See Details.
#'
#' @param compressTab Logical: Should the table be compressed? Default is TRUE,
#' in which case, each row is a unique codon. If FALSE, then each codon is represented
#' by 3 rows, one for each nucleotide comprising the codon.
#'
#' @param geneticCode Integer: A value relating to the \code{numcode} argument in
#' \code{seqinr::translate}.
#'
#' @details The argument \code{type} dictates what to pass to the argument
#' \code{dnaSeq}. If you want to translate a coding sequence (cDNA transcript),
#' then \code{type=='cds'} and \code{dnaSeq} must recieve a single character,
#' the DNA sequence.
#'
#' If you want to translate a series of exons, then \code{type='exon'} and
#' \code{dnaSeq} must receive a list, where each indexed item in the list is
#' a character vector, the DNA exon sequence. Note, it is assumed that the
#' exon sequences are ordered correctly, from first to last.
#'
#' For both cases, the function assumes that the sequence is in the correct
#' reading frame.
#'
#' @return
#' Returns a data.table with the following columns when \code{compressTab==TRUE}:
#'
#' \enumerate{
#'    \item \code{$CODON} = The codon number, 1:N.
#'    \item \code{$NUC} = The nucleotides positions comprising the codon, from
#'    1 to the gene's length.
#'    \item \code{$DNA} = The DNA bases.
#'    \item \code{$AMINO} = The amino acid residue.
#'    \item \code{$EXON} = The exon number, but only when \code{type=='exon'}.
#' }
#'
#' Otherwise, if \code{compressTab==FALSE}:
#'
#' \enumerate{
#'    \item \code{$CODON} = The codon number, 1:N.
#'    \item \code{$NUC.GENE} = The nucleotides position within the gene, from
#'    from 1 to the gene's length.
#'    \item \code{$NUC.CODON} = The nucleotides positions within the codon, from
#'    1 to 3.
#'    \item \code{$DNA} = The DNA bases.
#'    \item \code{$AMINO} = The amino acid residue.
#'    \item \code{$EXON} = The exon number, but only when \code{type=='exon'}.
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

dna2codonDT <- function(dnaSeq, type, compressTab=FALSE, geneticCode=1){
  require(data.table)
  require(seqinr)
  require(tidyverse)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ###   ASSERTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if(!type %in% c('cds','exon')){
    stop('Argument `type` must be one of "cds" or "exon". See dna2codonDT.')
  }

  if(type=='cds' & class(dnaSeq)!='character'){
    stop('Argument `dnaSeq` must be a single vector of character class when
         argument `type=="cds"`. See ?dna2codonDT.')
  }

  if(type=='exon' & class(dnaSeq)!='list'){
    stop('Argument `dnaSeq` must be a list of character class vectors when
         argument `type=="exon"`. See ?dna2codonDT.')
  }

  if(class(compressTab)!='logical'){
    stop('Argument `type` must be a logical class. See ?dna2codonDT.')
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ###   INTERNAL FUNCTION   ###
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  FUN_dna_2_cod <- function(seq_as_char, compressTab){
    # Split the sequence into codons
    codonList <- seqinr::splitseq(
      seq = seq_as_char  %>% as.character %>% s2c
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
          DNA=cod,
          AMINO=amino
        )
      } else if(compressTab==FALSE){
        # Uncompressed table format
        tab <- data.table(
          CODON=i,
          NUC.GENE=c(n-2, n-1, n),
          NUC.CODON=1:3,
          DNA=strsplit(cod, '')[[1]],
          AMINO=amino
        )
      }
    }) %>%
      do.call('rbind',.)

    return(codonTab)
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ###   EXECUTE   ###
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # If only providing coding sequence as single character string
  if(type=='cds'){
    result <- FUN_dna_2_cod(seq_as_char=dnaSeq, compressTab=compressTab)
  }

  # If providing exons as list of character strings
  if(type=='exon'){
    # Number of exons
    n <- length(dnaSeq)

    # Iterate over exons to create an index table for codons + nucleotides
    nucGeneIndex <- lapply(1:n, function(i){
      data.table(EXON=i, DNA=str_split(dnaSeq[[i]], '')[[1]])
    }) %>%
      do.call('rbind', .) %>%
      .[, NUC.GENE:=1:.N]

    # Get the codon table for the combined exons
    result <- paste(nucGeneIndex$DNA, collapse='') %>%
      FUN_dna_2_cod(seq_as_char=., compressTab=compressTab)

    # Add in the codons + nucleotides indexes for uncompressed and compressed
    # trasnlated codon data tables
    if(compressTab==FALSE){
      result <- left_join(result, nucGeneIndex)
    }

    if(compressTab==TRUE){
      # Align in the compressed nucleotide indexes against exons and codons
      nucGeneIndex <- result %>%
        .[, .(NUC.GENE=as.integer(unlist(str_split(NUC, '\\|')))), by=c('CODON','NUC')] %>%
        left_join(., nucGeneIndex, by='NUC.GENE')

      exonByCodon <- nucGeneIndex %>%
        .[, .(EXON=paste(sort(unique(EXON)),collapse='|')), by=CODON]

      # Combine
      result <- left_join(result, exonByCodon)
    }
  }

  # Output
  return(result)
}
