#' Filter missing data by loci
#'
#' Parses a data table of genotypes/allele frequencies and returns a list of
#' loci that conform to a desired missing data threshold.
#'
#' Note, it is assumed that missing data values have already been put is as
#' an \code{NA}. If this is not done in advance, this function will not produce
#' the expected results.
#'
#' @param dat Data table: The must contain the columns:
#' \enumerate{
#'    \item The sample ID (see param \code{sampCol}), for genotype datasets only.
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The population ID (see param \code{popCol}).
#'    \item The genotypes (see param \code{genoCol}), or the allele frequencies
#'    (see param \code{freqCol})
#' }
#'
#' @param missing Numeric: The proportion of missing data a locus, a value between
#' 0 and 1.
#'
#' @param type Character: Is \code{dat} a data table of genotypes (\code{'genos'})
#' or a data table of allele frequencies (\code{'freqs'})? Default = \code{'genos'}.
#'
#' @param method Character: The method by which missingness filtering is performed.
#' Only valid when filtering is performed on genotypes (\code{type=='genos'}).
#' One of \code{'samples'}, or \code{'pops'}. Default = \code{'samples'}.
#' For \code{'samples'}, missingness is calculated across all sampled individuals
#' (irrespective of their populations) for genotypes at each locus. If the
#' missingness summed across samples is greater than the threshold, the locus
#' will be discarded. For \code{'pops'}, missingness at a locus is calculated
#' per population. If any populations has missingness above the threshold at
#' a locus, then that locus will be removed.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default = \code{'SAMPLE'}. Only needed when \code{type=='genos'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param popCol Character: The column name with population information.
#' Default = \code{'POP'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Missing genotypes are encoded with an \code{NA}. Default = \code{'GT'}.
#' Only needed when \code{type=='genos'}.
#'
#' @param freqCol Character: The column name with the allele frequency information.
#' Missing frequencies are encoded with an \code{NA}. Default = \code{'freqCol'}.
#' Only needed when \code{type=='freqs'}.
#'
#' @details If \code{type=='genos'}, then your output will depend on how you
#' specify the \code{method} argument. If \code{type=='freqs'}, then there is
#' just one output, those loci with missing data less than the \code{missing} threshold.
#'
#' @return Returns a character vector of locus names in \code{dat$LOCUS} that conform
#' to the missingness threshold (<= to the value of \code{missing}).
#'
#' @examples
#' library(genomalicious)
#'
#' simMiss <- data_4pops %>% copy()
#' simMiss$GT[sample(1:nrow(simMiss), 0.1*nrow(simMiss), replace=FALSE)] <- NA
#'
#' filter_missing_loci(simMiss, 0.10)
#'
#' @export

filter_missing_loci <- function(
    dat, missing, type='genos', method='samples', sampCol='SAMPLE',
    locusCol='LOCUS', popCol='POP', genoCol='GT', freqCol='FREQ'
){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- library(data.table)
  for(L in libs){ require(L, character.only=TRUE)}

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){
    stop("Argument `dat` isn't a data table. See ?filter_maf.")
  }

  # Check columns for data.table
  dat.cols <- colnames(dat)
  if(type=='genos'){

    if(sum(c(sampCol, locusCol, genoCol, popCol) %in% colnames(dat))!=4){
      stop("Not all specified columns (`sampCol`, `locusCol`, `genoCol`, `popCol`) are in data.table dat. See ?filter_missing_loci.")
    }

    colnames(dat)[match(c(sampCol, locusCol, genoCol, popCol),dat.cols)] <-
      c('SAMPLE','LOCUS','GT','POP')
  } else if (type=='freqs'){

    if(sum(c(locusCol, freqCol, popCol) %in% colnames(dat))!=3){
      stop("Not all specified columns (`locusCol`, `freqCol`, `popCol`) are in data.table dat. See ?filter_missing_loci.")
    }

    colnames(dat)[match(c(locusCol, freqCol, popCol),dat.cols)] <-
      c('LOCUS','FREQ','POP')
  }

  # Check that missing is between 0 and 1
  if(missing < 0 | missing > 1){
    stop("Argument `missing` must be between 0 and 1. See ?filter_missing_loci")
  }

  # Check that type is specified correctly.
  if(!type %in% c('genos','freqs')){
    stop("Argument `type' must be one of 'genos' or 'freqs'. See ?filter_missing_loci.")
  }

  # Check that method is specified correctly.
  if(type=='genos' & !method %in% c('samples','pops')){
    stop("Argument `method` must be one of 'samples' or 'pops'. See ?filter_missing_loci.")
  }

  # If genotypes, convert characters to integers
  if(class(dat$GT)=='character'){
    dat[, GT:=genoscore_converter(GT)]
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  if(type=='genos'){
    if(method=='samples'){
      dMiss <- dat[, .(MISS=sum(is.na(GT))/length(GT)), by='LOCUS']
      good.loci <- dMiss[MISS <= missing]$LOCUS
    } else if (method=='pops'){
      num.pops <- length(unique(dat$POP))
      dMiss <- dat[, .(MISS=sum(is.na(GT))/length(GT)), by=c('LOCUS','POP')]
      good.loci <- dMiss %>%
        .[MISS <= missing] %>%
        .[, .(NUM.POPS=sum(length(unique(POP)))), by='LOCUS'] %>%
        .[NUM.POPS == num.pops,] %>%
        .[['LOCUS']]
    }
  } else if(type=='freqs'){
    dMiss <- dat[, .(MISS=sum(is.na(FREQ))/length(FREQ)), by='LOCUS']
    good.loci <- dMiss[MISS <= missing]$LOCUS
  }

  # Output
  return(good.loci)

  # ........ END
}
