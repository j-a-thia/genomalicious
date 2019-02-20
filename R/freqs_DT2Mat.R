#' Convert a data table of allele frequencies into a matrix (or vice versa)
#'
#' Takes a data table of allele frequencies in long format and converts it into
#' a matrix in wide format (loci in columns and populations in rows). The reverse
#' can also be done.
#'
#' @param dat Data table or Matrix: The object to transform. If this is a long data table
#' # of allele frequencies, then three columns are required:
#' \enumerate{
#'    \item (1) The population ID (see param \code{popCol}).
#'    \item (2) The locus ID (see param \code{lociCol}).
#'    \item (3) The Ref allele frequency (see param \code{freqCol}).
#' }
#' The population pool column serves as the pivot point to convert the long data table into a wide matrix.
#' If convertin from a frequency matrix to a data table, see argument \code{flip}.
#'
#' @param popCol Character: The column name with the population information.
#' Default = \code{'POOL'} (e.g. for pool-seq data).
#'
#' @param lociCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param freqCol Character: The column name with the Ref allele frequency.
#' Default = \code{'FREQ'}.
#'
#' @param flip Logical: Instead of converting a (long) data table to a (wide) matrix,
#' should a (wide) matrix be converted into a (long) data table? Default = \code{FALSE}.
#' If \code{TRUE}, then param \code{dat} must be a matrix, with loci names as column headers,
#' population IDs in the row names, and frequencies in the cells. When \code{TRUE}, params
#' \code{popCol}, \code{lociCol}, and \code{freqCol} become void.
#'
#' @return When \code{flip=FALSE}, converts a data table into a frequency matrix. When
#' \code{flip=TRUE}, converts a matrix into a data table with three columns: (1) \code{$POOL},
#' the population ID (as for pool-seq data); (2) \code{$LOCUS}, the locus ID; and (3) \code{FREQ},
#' the Ref allele frequency.
#'
#' @examples
#' data(pgposerPi)
#'
#' # Convert a long data table to a wide matrix
#' freqMat <- freq_DT2Mat(pgposerPi, popCol='POOL', lociCol='LOCUS', freqCol='PI', flip=FALSE)
#'
#' # Convert a wide matrix back to a data table
#' freqDT <- freq_DT2Mat(freqMat, flip=TRUE)
#'
#' export
freqs_DT2Mat <- function(dat, popCol='POOL', lociCol='LOCUS', freqCol='FREQ', flip=FALSE){
  # BEGIN ............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('tidyr', 'data.table')){require(lib, character.only=TRUE)}

  # If going from long data table to wide matrix, make sure dat is a data table.
  if(flip==FALSE & !'data.table' %in% class(dat)){
    stop("Argument flip==FALSE, but class(dat) isn't 'data.table'.")
  }

  # If going from a wide matrix to long data table, make sure dat is a matrix.
  if(flip==TRUE & !'matrix' %in% class(dat)){
    stop("Argument flip==TRUE, but class(dat) isn't 'matrix'.")
  }

  # If providing a data table, check that popCol, lociCol, and freqCol are in dat.
  if(class(dat)[1]=='data.table'){
    if(length(which((c(popCol, lociCol, freqCol) %in% colnames(dat))==FALSE)) > 0){
      stop("Argument dat does not have columns specified in arguments popCol, lociCol, or freqCol.")
    }
  }

  # If providing a matrix, check that their are row names.
  if(class(dat)[1]=='matrix' & is.null(rownames(dat))==TRUE){
    stop("Argument dat is frequency matrix, but has no population IDs in the row names.")
  }


  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  if(flip==FALSE){
    freqDT <- spread(dat[, c(popCol, freqCol, lociCol), with=FALSE], key=lociCol, value=freqCol)
    freqMat <- as.matrix(freqDT[, !popCol, with=FALSE])
    rownames(freqMat) <- freqDT[[popCol]]
    return(freqMat)
  } else if(flip==TRUE){
      freqDT <- data.table(dat, keep.rownames=TRUE)
      colnames(freqDT)[which(colnames(freqDT)=='rn')] <- 'POOL'
      locusNames <- colnames(freqDT)[which(colnames(freqDT)!='POOL')]
      freqDT <- melt(freqDT, id.vars='POOL', variable='LOCUS', value='FREQ')
      return(freqDT)
  }

}
