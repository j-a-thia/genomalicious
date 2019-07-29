#' Convert a data table of genotypes into a matrix (or vice versa)
#'
#' Takes a data table of genotypes in long format and converts it into
#' a matrix in wide format (loci in columns and individuals in rows). The reverse
#' can also be done. See also \code{DT2Mat_freqs} for converting matrix of frequencies.
#'
#' @param dat Data table or Matrix: The object to transform. If this is a long data table
#' # of genotypes coded as per VCF specifications, e.g. ('0/0', '0/1', '1/1'), or counts
#' of the Ref alleles e.g. (0, 1, 2).
#' Three columns are required:
#' \enumerate{
#'    \item The sampled individual ID (see param \code{sampCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The genotype (see param \code{genoCol}).
#' }
#' The sampled individual ID column serves as the pivot point to convert the long data table into a wide matrix.
#' If converting from a genotypes matrix to a data table, see argument \code{flip}.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#'
#' @param locusCol Character: The column name with the locus information.
#'
#' @param genoCol Character: The column name with the genotype information.
#'
#' @param genoScore Character: Should the returned genotypes by represented as
#' \code{'counts'} of the Ref allele (integer from 0 -> 2), or separated alleles,
#' \code{'sep'} (character, '0/0', '0/1', or '1/1')? Default = \code{'sep'}.
#'
#' @param flip Logical: Instead of converting a (long) data table to a (wide) matrix,
#' should a (wide) matrix be converted into a (long) data table? Default = \code{FALSE}.
#' If \code{TRUE}, then param \code{dat} must be a matrix, with loci names as column headers,
#' sample IDs in the row names, and genotypes in the cells. When \code{TRUE}, params
#' \code{sampCol}, \code{locusCol}, and \code{genoCol} are used to structure the new matrix.
#'
#' @return When \code{flip=FALSE}, converts a data table into a genotype matrix. When
#' \code{flip=TRUE}, converts a matrix into a data table with three columns: (1) \code{$SAMPLE},
#' the sample ID; (2) \code{$LOCUS}, the locus ID; and (3) \code{GT},
#' the Ref allele frequency.
#'
#' @examples
#' data(genomalicious_4pops)
#' datGt <- genomalicious_4pops[LOCUS %in% unique(genomalicious_4pops$LOCUS)[1:8]]
#'
#' # Convert a long data table to a wide matrix
#' genoMatSep <- DT2Mat_genos(datGt
#'               , sampCol='SAMPLE'
#'               , locusCol='LOCUS'
#'               , genoCol='GT'
#'               , flip=FALSE)
#'
#' genoMatCounts <- DT2Mat_genos(datGt
#'               , sampCol='SAMPLE'
#'               , locusCol='LOCUS'
#'               , genoCol='GT'
#'               , genoScore='counts'
#'               , flip=FALSE)
#'
#' # Convert a wide matrix back to a data table
#' genoDTSep <- DT2Mat_genos(genoMatSep
#'               , sampCol='SAMPLE'
#'               , locusCol='LOCUS'
#'               , genoCol='GT'
#'               , genoScore='sep'
#'               , flip=TRUE)
#'
#' @export
DT2Mat_genos <- function(dat, sampCol=NA, locusCol=NA, genoCol=NA, genoScore='sep', flip=FALSE){

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

  # If providing a matrix, check that their are row names.
  if(class(dat)[1]=='matrix' & is.null(rownames(dat))==TRUE){
    stop("Argument dat is a genotype matrix, but has no individual IDs in the row names.")
  }

  # If providing a data table, check that sampCol, locusCol, and genoCol are in dat.
  if(class(dat)[1]=='data.table'){
    if(length(which((c(sampCol, locusCol, genoCol) %in% colnames(dat))==FALSE)) > 0){
      stop("Argument dat does not have columns specified in arguments sampCol, locusCol, or genoCol.")
    }
  }

  # Check that genoScore option specified properly
  if(!genoScore %in% c('counts', 'sep')){
    stop("Argument genoScore must take the values of either 'counts' or 'sep'.")
  }

  # Check the column arguments are specified
  if(flip==FALSE){
    if(is.na(sampCol)){
      stop("Argument sampCol unspecified.")
    }

    if(is.na(locusCol)){
      stop("Argument locusCol unspecified.")
    }

    if(is.na(genoCol)){
      stop("Argument genoCol unspecified.")
    }
  }

  # Evaluate the genotype data type.
  # genoInitClass stores the initial genotype score information.
  if(flip==FALSE){
    genoInitClass <- class(dat[[genoCol]])
  } else if(flip==TRUE){
    genoInitClass <- class(dat[, 1])
  }

  if(!genoInitClass %in% c('character', 'integer', 'numeric')){
    stop('Genotypes must be either a character or integer class.')
  }

  # If genotypes counts are numerics, convert to integers.
  if(genoInitClass=='numeric'){
    if(class(dat)[1]=='data.table'){
      dat[[genoCol]] <- as.integer(dat[[genoCol]])
    } else if(class(dat=='matirx')){
      dat <- apply(dat, 2, as.integer)
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  if(flip==FALSE){
    genoMat <- spread(dat[, c(sampCol, locusCol, genoCol), with=FALSE], key=locusCol, value=genoCol)
    sampVals <- genoMat[[sampCol]]
    genoMat <- as.matrix(genoMat[, !(sampCol), with=FALSE])
    row.names(genoMat) <- sampVals

    # Return genos as separated alleles, per VCF format?
    # Or return as counts of the Ref alleles?
    if(genoScore=='sep' & genoInitClass=='character'){
      return(genoMat)
    } else if(genoScore=='counts' & genoInitClass=='integer'){
      return(genoMat)
    } else if(genoScore=='counts' & genoInitClass=='character'){
      genoMat <- apply(genoMat, 2, genoscore_converter)
      return(genoMat)
    } else if(genoScore=='sep' & genoInitClass=='integer'){
      genoMat <- apply(genoMat, 2, genoscore_converter)
      return(genoMat)
    }

  }

  if(flip==TRUE){
    # Convert the matrix into a data table, keeping row names
    genoDT <- data.table(dat, keep.rownames=TRUE)
    # The row names are turned into a column 'rn', replace.
    colnames(genoDT)[which(colnames(genoDT)=='rn')] <- sampCol
    # Rejig data table
    genoDT <- melt(genoDT, id.vars=sampCol, variable.name=locusCol, value.name=genoCol)
    # Adjust genotype score?
    if(genoScore=='sep' & genoInitClass=='character'){
        return(genoDT)
      } else if(genoScore=='counts' & genoInitClass=='integer'){
        return(genoDT)
      } else if(genoScore=='counts' & genoInitClass=='character'){
        genoDT$GT <- genoscore_converter(genoDT$GT)
        return(genoDT)
      } else if(genoScore=='sep' & genoInitClass=='integer'){
        genoDT$GT <- genoscore_converter(genoDT$GT)
        return(genoDT)
      }
    }
}


