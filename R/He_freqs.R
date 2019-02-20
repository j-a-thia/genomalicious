#' Calculate expected heterozygosity from allele frequencies
#'
#' Takes a matrix of allele frequencies for different population pools,
#' at different loci, and returns various statistics for expected heterozygosity (He).
#' Assumes biallelic data.
#'
#' @param dat Matrix: Rows are population pools, columns are loci, cells are the Ref allele
#' frequencies.
#'
#' @details Calculates He as 2 * p * q; where p = the Ref allele frequency, and
#' q = the Alt allele frequency.
#'
#' @return Returns a list with the following indexed data tables: \cr
#' \enumerate{
#'     \item \code{poolAv} = The overall pool mean across loci.
#'     \item \code{locusAv} = The He for each locus and pool combination.
#' }
#' For each itemed index, there is a data table with columns \code{$POOL} (= population
#' pool ID) and \code{$HE} (= the expected heterozygosity). Additionally, in the
#' index \code{$locusAv}, there is a column \code{$LOCUS} (= the locus ID).
#'
#' @examples
#' data(pgposerFreqs)
#'
#' He <- He_freqs(pgposerFreqs)
#'
#' @export
He_freqs <- function(dat){
  # BEGIN .............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table)

  if(!'matrix' %in% class(dat)){ stop("Argument dat isn't a matrix.")}

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
  locusAv <- as.data.frame(t(apply(dat, 1, FUN_2pq)))
  locusAv$POOL <- rownames(locusAv)
  locusAv <- data.table(melt(locusAv, id.vars='POOL', variable.name='LOCUS', value.name='HE'))
  poolAv <- locusAv[, mean(HE), by='POOL']
  colnames(poolAv)[colnames(poolAv)=='V1'] <- 'HE'

  return(list(poolAv=poolAv, locusAv=locusAv))

  # ............ END
}
