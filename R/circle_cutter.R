#' Cut a replicated circular sequence to generate a linear sequence
#'
#' This function cuts a circular sequence using a character string motif to
#' generate a single linear sequence. This is primarily intended for "read-through"
#' type sequences generated, for example, from assembly or long-read sequencing
#' of circular genomes. Such sequences typically have replicated regions
#' that need to be removed to generate a singular un-replicated linear sequence.
#'
#' @param query_seq Character: The query sequence.
#'
#' @param motif Character, the sequence used to identify cut points in the
#' circular chromosome and to remove replicated regions. See Details.
#'
#' @details The argument \code{motif} is used to identify a starting and end
#' positions to cut the circular chromosome and remove replciated regions.
#' If for example \code{motif=='AATTGGCC'} and the sequence in question was: \cr
#' \cr
#' AATTGGCC ACTATCTGCTAGCTAGCATAGCATCGATCAGCATGACGCGCAA AATTGGCC \cr
#' \cr
#' The function will cut the sequence like so (marked with '|'): \cr
#' \cr
#' | AATTGGCC ACTATCTGCTAGCTAGCATAGCATCGATCAGCATGACGCGCAA | AATTGGCC
#'
#' @return Returns the subset sequence as a character string.
#'
#' @export
#'
#' @examples
#' x <- 'AATTGGCCACTATCTGCTAGCTAGCATAGCATCGATCAGCATGACGCGCAAAATTGGCC'
#'
#' # Find character motif that is repeated
#' motif_hits <- circularity_test(x, word_size = 8)
#' motif_seq <- substr(x, motif_hits[1,1], motif_hits[1,2])
#'
#' circle_cutter(query_seq = x, motif=motif_seq)
#'
circle_cutter <- function(query_seq, motif){

  # Packages
  require(stringr); require(tidyverse)

  # Get the sequence sans header
  query_seq <- as.character(query_seq)

  if(length(query_seq)>1){ query_seq <- paste(query_seq, collapse='')}

  # Identify all start points of the motif.
  motif_locs <- str_locate_all(query_seq, motif)[[1]]

  motif_st1 <- motif_locs[1, 'start']
  motif_st2 <- motif_locs[2, 'start']

  # Subset
  seq_sub <- str_sub(query_seq, motif_st1, motif_st2-1)

  # Write out
  return(seq_sub)
}
