#' Reposition a mitogenomic sequence
#' 
#' This function takes a linear sequence of a mitochondrial genome and 
#' repositions it to start at a new base position.
#' 
#' @param dnaSeq Character: The mitogenomic sequence.
#' 
#' @param newStartPos Integer: The new start position for the sequence.
#' 
#' @details The sequence is repositioned with respect to the value of 
#' \code{newStartPos}. For example, consider a value of \code{dnaSeq}: \cr
#' \cr
#'      GCGCGCGCGCATATGACTAA \cr
#' \cr
#' For a value of \code{newStartPos==10}, the 10th base would become the new
#' start position and the 9th position in the original sequence is a cutoff.
#' Below this cutoff, the sequence is translocated to the right. \cr
#' \cr
#'      GCGCGCGCG | CATATGACTAA (original sequence, cut at 9th base) \cr\cr
#'      CATATGACTAAGCGCGCGCG (new sequence, repositioned)
#'
#' @export
#' 
mitogenome_seq_repos <- function(dnaSeq, newStartPos){
  require(stringr)
  
  # Get the new left and right sides of the sequence.
  new.seq.L <- stringr::str_sub(dnaSeq, newStartPos, nchar(dnaSeq))
  new.seq.R <- stringr::str_sub(dnaSeq, 1, newStartPos-1)
  
  # Join the left and right sides
  new.seq.join <- paste(new.seq.L, new.seq.R, sep='')
  
  # The length of each side.
  new.len.L <- nchar(new.seq.L)
  new.len.R <- nchar(new.seq.R)
  
  # Output
  return(new.seq.join)
}
