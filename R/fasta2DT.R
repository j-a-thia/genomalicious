#' FASTA file to data table
#'
#' Reads a FASTA file and converts to a data table.
#'
#' @param fastaFile Character: The path to the input FASTA file.
#'
#' @return A data table with the columns \code{$LOCUS} (the locus name) and \code{$SEQ} (the sequence).
#'
#' @export
fasta2DT <- function(fastaFile){

  # BEGIN ............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table)

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Read in the FASTA as lines
  fastaLines <- readLines(fastaFile)

  # Determine position of sequence names
  seqIDpos <- grep(pattern='>', x=fastaLines)

  # Turn into a data table
  fastaDT <- data.table(LOCUS=fastaLines[seqIDpos], SEQ=fastaLines[seqIDpos+1])

  # Adjust LOCUS names and return
  fastaDT$LOCUS <- gsub(pattern='>', x=fastaDT$LOCUS, replacement='')
  return(fastaDT)
}
