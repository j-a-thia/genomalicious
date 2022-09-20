#' Convert a data table of allele frequencies into a matrix (or vice versa)
#'
#' Takes a data table of allele frequencies in long format and converts it into
#' a matrix in wide format (loci in columns and populations in rows). The reverse
#' can also be done. See also \code{DT2Mat_genos} for converting matrix of genotypes.
#'
#' @param dat Data table or Matrix: The object to transform. If this is a long data table
#' of allele frequencies, then three columns are required:
#' \enumerate{
#'    \item The population ID (see param \code{popCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The Ref allele frequency (see param \code{freqCol}).
#' }
#' The population pool column serves as the pivot point to convert the long data table into a wide matrix.
#' If converting from a frequency matrix to a data table, see argument \code{flip}.
#'
#' @param popCol Character: The column name with the population information. Default is 'POP'.
#'
#' @param locusCol Character: The column name with the locus information. Default is 'LOCUS'.
#'
#' @param freqCol Character: The column name with the Ref allele frequency. Default is 'FREQ'.
#'
#' @param flip Logical: Instead of converting a (long) data table to a (wide) matrix,
#' should a (wide) matrix be converted into a (long) data table? Default = \code{FALSE}.
#' If \code{TRUE}, then param \code{dat} must be a matrix, with loci names as column headers,
#' population IDs in the row names, and frequencies in the cells. When \code{TRUE}, params
#' \code{popCol}, \code{locusCol}, and \code{freqCol} become void.
#'
#' @return When \code{flip=FALSE}, converts a data table into a frequency matrix. When
#' \code{flip=TRUE}, converts a matrix into a data table with three columns: (1) \code{$POP},
#' the population ID (as for pool-seq data); (2) \code{$LOCUS}, the locus ID; and (3) \code{FREQ},
#' the Ref allele frequency.
#'
#' @examples
#' data(data_PoolFreqs)
#'
#' # Convert a long data table to a wide matrix
#' freqMat <- DT2Mat_freqs(
#'    dat=data_PoolFreqs,
#'    popCol='POOL',
#'    locusCol='LOCUS',
#'    freqCol='FREQ',
#'    flip=FALSE)
#'
#' freqMat
#'
#' # Convert a wide matrix back to a data table
#' freqDT <- DT2Mat_freqs(
#'    freqMat,
#'    popCol='POOL',
#'    locusCol='LOCUS',
#'    freqCol='FREQ',
#'    flip=TRUE)
#'
#' freqDT
#'
#' @export
DT2Mat_freqs <- function(dat, popCol='POP', locusCol='LOCUS', freqCol='FREQ', flip=FALSE){
  # BEGIN ............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('tidyr', 'data.table')){require(lib, character.only=TRUE)}

  # If going from long data table to wide matrix, make sure dat is a data table.
  if(flip==FALSE & !'data.table'%in%class(dat)){
    stop("Argument flip==FALSE, but class(dat) isn't 'data.table': see ?DT2Mat_freqs")
  }

  # If going from a wide matrix to long data table, make sure dat is a matrix.
  if(flip==TRUE & !'matrix'%in%class(dat)){
    stop("Argument flip==TRUE, but class(dat) isn't 'matrix': see ?DT2Mat_freqs")
  }

  # If providing a matrix, check that their are row names.
  if('matrix'%in%class(dat) & is.null(rownames(dat))==TRUE){
    stop("Argument dat is frequency matrix, but has no population IDs in the row names: see ?DT2Mat_freqs")
  }

  # If providing a data table, check that popCol, locusCol, and freqCol are in dat.
  if('data.table'%in%class(dat)){
    if(length(which((c(popCol, locusCol, freqCol) %in% colnames(dat))==FALSE)) > 0){
      stop("Argument dat does not have columns specified in arguments popCol, locusCol, or freqCol: see ?DT2Mat_freqs")
    }
  }

  # Check the column arguments are specified
  if(flip==FALSE){
    if(is.na(popCol)){
      stop("Argument popCol unspecified: see ?DT2Mat_freqs")
    }

    if(is.na(locusCol)){
      stop("Argument locusCol unspecified: see ?DT2Mat_freqs")
    }

    if(is.na(freqCol)){
      stop("Argument freqCol unspecified: see ?DT2Mat_freqs")
    }
  }


  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  if(flip==FALSE){
    # Spread out the data table
    freqDT <- spread(dat[, c(popCol, freqCol, locusCol), with=FALSE], key=locusCol, value=freqCol)
    # Get the population column values
    popVals <- freqDT[[popCol]]
    # Turn the data table into a matrix
    freqMat <- as.matrix(freqDT[, !popCol, with=FALSE])
    # Add population values as rows
    rownames(freqMat) <- popVals
    return(freqMat)

  } else if(flip==TRUE){
      # Turn matrix into data table, keep row names
      freqDT <- data.table(dat, keep.rownames=TRUE)
      # The row names are turned into a column 'rn', replace.
      colnames(freqDT)[which(colnames(freqDT)=='rn')] <- popCol
      # Rejig the data table
      freqDT <- melt(freqDT, id.vars=popCol, variable.name=locusCol, value=freqCol)
      return(freqDT)
  }

}
