#' Calculate expected heterozygosity from allele frequencies
#'
#' Takes a matrix or data table/frame of allele frequencies for different populations,
#' at different loci, and returns various statistics for expected heterozygosity (He).
#' Assumes biallelic data.
#'
#' @param dat Matrix or Data table/frame: Can take a matrix of allele frequencies (row names
#' are populations, column names are loci, cells are frequencies of Ref alleles) or
#' a data table/frame in long format (requires)
#'
#' @param popCol Character: The column name in the long format data table/frame \code{dat}
#' corresponding to the population ID.
#'
#' @param locusCol Character: The column name in the long format data table/frame \code{dat}
#' corresponding to the locus ID.
#'
#' @param freqCol Character: The column name in the long format data table/frame \code{dat}
#' corresponding to the frequency of the Ref allele for a locus in a population..
#'
#' @details Calculates He as 2 * p * q; where p = the Ref allele frequency, and
#' q = the Alt allele frequency.
#'
#' @return Returns a list with the following indexed items (data tables): \cr
#' \enumerate{
#'     \item \code{popXlocus} = He for each population x locus combination (columns:
#'     \code{$POP}, \code{$LOCUS}, \code{$HE}).
#'     \item \code{popAv} = Average He per populations (columns: \code{$POP}, \code{$HE.AV}).
#'     \item \code{locusAv} = Average He per locus (columns: \code{$LOCUS}, \code{$HE.AV})
#' }
#'
#' @examples
#' # Data as matrix
#' data(data_FreqsMat)
#' He1 <- He_freqs(data_FreqsMat)
#'
#' # Data as long data table
#' data(data_PoolFreqs)
#' He2 <- He_freqs(data_PoolFreqs, popCol='POOL', locusCol='LOCUS', freqCol='PI')
#'
#' @export
He_freqs <- function(dat, popCol='POP', locusCol='LOCUS', freqCol='FREQ'){
  # BEGIN .............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table)

  if(!class(dat)[1] %in% c('matrix', 'data.table', 'data.frame')){
    stop("Argument dat isn't a matrix or a data table/frame")
  }

  if('data.frame' %in% class(dat)){ dat <- as.data.table(dat) }

  if('data.table'%in%class(dat) & is.na(popCol)){
    stop("Argument popCol unspecified.")
  }

  if('data.table'%in%class(dat) & is.na(locusCol)){
    stop("Argument locusCol unspecified.")
  }

  if('data.table'%in%class(dat) & is.na(freqCol)){
    stop("Argument freqCol unspecified.")
  }

  # --------------------------------------------+
  # Internal function
  # --------------------------------------------+
  # Use this to calculate the expected heterozygosity from a single Ref value, p.
  FUN_2pq <- function(p){
    return(2 * p * (1 - p))
  }

  # --------------------------------------------+
  # Analysis
  # --------------------------------------------+
  # Get the He for each population x locus combination
  if(class(dat)[1]=='matrix'){
    popXlocus <- as.data.table(apply(dat, 2, FUN_2pq))
    popXlocus$POP <- rownames(dat)
    popXlocus <- melt(popXlocus, id='POP', variable.name='LOCUS', value.name='HE')
  }

  if('data.table'%in%class(dat)){
    colReass <- match(c(popCol, locusCol, freqCol), colnames(dat))
    colnames(dat)[colReass] <- c('POP', 'LOCUS', 'FREQ')
    popXlocus <- dat[, FUN_2pq(FREQ), by=c('POP', 'LOCUS')]
    setnames(popXlocus, 'V1', 'HE')
  }

  # Get the population and locus averages
  popAv <- popXlocus[, mean(HE), by=POP]; setnames(popAv, 'V1', 'HE.AV')
  locusAv <- popXlocus[, mean(HE), by=LOCUS]; setnames(locusAv, 'V1', 'HE.AV')

  # Return a list
  return(list(popXlocus=popXlocus, popAv=popAv, locusAv=locusAv))

  # ............ END
}
