#' Translate a list of coding sequences into a codon data.table
#'
#' Convert a DNA coding sequence into a data.table of codons, nucleotide indexes,
#' DNA triplet sequences, and amino acids. If provided with coding sequences split
#' into exons, these will also be incorporated into the table.
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
#'    \item \code{$NUC.GENE} = The nucleotides positions comprising the codon, from
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

  if(!type %in% c('cds','exons')){
    stop('Argument `type` must be one of "cds" or "exon". See dna2codonDT.')
  }

  if(type=='cds' & class(dnaSeq)!='character'){
    stop('Argument `dnaSeq` must be a single vector of character class when
         argument `type=="cds"`. See ?dna2codonDT.')
  }

  if(type=='exons' & class(dnaSeq)!='list'){
    stop('Argument `dnaSeq` must be a list of character class vectors when
         argument `type=="exons"`. See ?dna2codonDT.')
  }

  if(class(compressTab)!='logical'){
    stop('Argument `type` must be a logical class. See ?dna2codonDT.')
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ###   EXECUTE   ###
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ### IF SINGLE CODING SEQUENCE PROVIDED
  if(type=='cds'){
    codonIndex <- seqinr::s2c(dnaSeq) %>%
      seqinr::splitseq(., frame=0, word=3) %>%
      data.table(TRIPLET=.) %>%
      .[, CODON:=1:.N] %>%
      .[, .(DNA=seqinr::s2c(TRIPLET)), by=CODON] %>%
      .[, NUC.GENE:=1:.N] %>%
      .[, NUC.CODON:=1:.N, by=CODON] %>%
      .[, AMINO:=seqinr::translate(DNA, numcode=geneticCode), by=CODON]
  }

  ### IF EXON LIST OF SEQUENCES PROVIDED
  if(type=='exons'){
    # Number of exons
    n <- length(dnaSeq)

    # Iterate over exons to create an exon index table
    exonIndex <- lapply(1:n, function(i){
      data.table(EXON=i, DNA=seqinr::s2c(dnaSeq[[i]]))
    }) %>%
      do.call('rbind', .) %>%
      .[, NUC.GENE:=1:.N]

    # Iterate over DNA bases to create a codon index table (with translated)
    # amino acids
    codonIndex <- seqinr::splitseq(exonIndex$DNA, frame=0, word=3) %>%
      data.table(TRIPLET=.) %>%
      .[, CODON:=1:.N] %>%
      .[, .(DNA=seqinr::s2c(TRIPLET)), by=CODON] %>%
      .[, NUC.GENE:=1:.N] %>%
      .[, NUC.CODON:=1:.N, by=CODON] %>%
      .[, AMINO:=seqinr::translate(DNA, numcode=geneticCode), by=CODON]

    # Combine the exon and codon index table together
    result <- left_join(
      exonIndex[, c('EXON','NUC.GENE')],
      codonIndex[, c('CODON','NUC.GENE','NUC.CODON','DNA','AMINO')]
    )
  }

  # Compress the result if that is the desired format
  if(compressTab==TRUE){
    # Collapse nucleotide index, DNA triplets, and amino acids
    nucCollapse <- result %>%
      .[, .(NUC.GENE=paste(sort(NUC.GENE),collapse='|')),by=CODON]
    dnaCollapse <- result %>%
      copy %>%
      setorder(., NUC.GENE) %>%
      .[, .(DNA=paste(DNA,collapse='')),by=CODON]
    aminoUniq <- result %>%
      .[, c('CODON','AMINO')] %>%
      unique()

    # An extra step if exons provided
    if(type=='exon'){
      exonCollapse <- result %>%
        .[, .(EXON=paste(sort(unique(EXON)),collapse='|')),by=CODON]

      result <- left_join(exonCollapse, nucCollapse) %>%
        left_join(., dnaCollapse) %>%
        left_join(., aminoUniq) %>%
        as.data.table()
    } else if(type=='cds'){
      # Otherwise...
      result <- left_join(nucCollapse, dnaCollapse) %>%
        left_join(., aminoUniq) %>%
        as.data.table()
    }
  }

  ### OUTPUT
  return(result)
}
