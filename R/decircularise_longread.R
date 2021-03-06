#' Decircularise a circular chromosome generated from long reads
#'
#' Long-read technologies can provide complete read throughs of circular
#' chromosomes. However, when the number of bases sequenced is greater than the
#' chromosome size, this can lead to sequence replication within the long reads.
#' This function makes the circule readthrough linear. It requires a
#' user-specified motif to identify replicated regions and remove them.
#'
#' @param fastaIn Character, the path to a FASTA file. Must only contain
#' one sequence.
#'
#' @param fastaOut Character, the path to write out the cut FASTA.
#'
#' @param motif Character, the sequence used to identify cut points in the
#' circular chromosome and to remove replicated regions. See Details.
#'
#' @param newHead Character, an optional new header to use for the cut FASTA
#' sequence. Default is NULL, in which case, the same header is used as per
#' the first line in \code{fastaIn}. You do not need to specify ">".
#'
#' @details The argument \code{motif} is used to identify a starting and end
#' positions to cut the circular chromosome and remove replciated regions.
#' If for example \code{motif=='AATTGGCC'} and the sequence in question was: \cr
#' \cr
#' GCGCAA AATTGGCC ACTATCTGCTAGCTAGCATAGCATCGATCAGCATGAC GCGCAA AATTGGCC \cr
#' \cr
#' The function will cut the sequence like so (marked with '|'): \cr
#' \cr
#' GCGCAA | AATTGGCC ACTATCTGCTAGCTAGCATAGCATCGATCAGCATGAC GCGCAA | AATTGGCC
#'
#' @export
#'
decircularise_longread <- function(fastaIn, fastaOut, motif, newHead=NULL){

  # Packages
  require(stringr); require(tidyverse)

  # Read FASTA lines
  fa_read <- readLines(fastaIn)

  # The number of sequence headers
  fa_heads <- grep('>', fa_read)

  # Throw error if there are >1 sequqneces
  if(length(fa_heads)>1){
    stop('There are more than 1 sequences in file specified by fastaFile.
         See ?decircularise_longread.')
  }

  # Get the sequence sans header
  fa_seq <- fa_read[2:length(fa_read)]

  if(length(fa_seq)>1){ fa_seq <- paste(fa_seq, collapse='')}

  # Identify all start points of the motif.
  motif_locs <- str_locate_all(fa_seq, motif)[[1]]

  motif_st1 <- motif_locs[1, 'start']
  motif_st2 <- motif_locs[2, 'start']

  # Subset
  seq_sub <- str_sub(fa_seq, motif_st1, motif_st2-1)

  # Write out
  if(is.null(newHead)){ out_head <- fa_read[1]
  } else { out_head <- paste('>', newHead, sep='') }

  writeLines(c(out_head, seq_sub), fastaOut)
}
