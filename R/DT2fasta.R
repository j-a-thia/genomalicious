#' Data table to FASTA file
#'
#' Converts a data table of sequences (e.g. like that produced by 
#' \code{fasta2DT()}) and writes a FASTA file
#' 
#' @param dat Data table: Contains columns \code{$LOCUS} of locus names
#' and \code{$SEQ} of nucleotide sequences.
#'
#' @param fastaFile Character: The path to the output FASTA file.
#' 
#' @export
DT2fasta <- function(dat, fastaFile){
  # Add the '>' to locus column
  dat$LOCUS <- paste0('> ', dat$LOCUS)
  
  # Convert to lines. Transpose dat to read rows in.
  fastaLines <- unlist(t(dat))
  
  # Open connection and write.
  fileCon <- file(fastaFile)
  writeLines(fastaLines, fileCon, sep='\n')
  close(fileCon)
}