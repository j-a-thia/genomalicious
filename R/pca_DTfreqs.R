#' Conduct a PCA on allele frequences stored in a (long) data table
#'
#' Takes a long data table of biallelic population allele frequencies
#' and conducts a PCA using R's \code{prcomp()} function.
#' Allele frequencies are scaled to a mean = 0 and variance = 1.
#'
#' @param dat Data table: A long data table with the following columns:
#' \enumerate{
#'    \item (1) The population ID (see param \code{popCol}).
#'    \item (2) The locus ID (see param \code{locusCol}).
#'    \item (3) The Ref allele frequency (see param \code{freqCol}).
#' }
#'
#' @param popCol Character: An optional argument. The column name with the
#' population information. Default is \code{'POP'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default is \code{'LOCUS'}.
#'
#' @param freqCol Character: The column name with the genotype information.
#' Default is \code{'GT'}.
#'
#' @return Returns a \code{prcomp} object with an additional index of \code{$pops},
#' the population ID.
#'
#' @references
#' Patterson et al. (2006) Population structure and eigenanalysis. PLOS Genetics.
#'
#' @examples
#' data(genomaliciousFreqsLong)
#'
#' # Conduct the PCA
#' pca <- pca_DTfreqs(genomaliciousFreqsLong)
#'
#' # The PCA scores
#' pca$x
#'
#' # The populations, i.e. rows in pca$x
#' pca$pops
#'
#' @export
pca_DTfreqs <- function(dat, popCol='POP', locusCol='LOCUS', freqCol='FREQ'){
  freqMat <- DT2Mat_freqs(dat, popCol=popCol, locusCol=locusCol, freqCol=freqCol)

  pca <- prcomp(freqMat, center=TRUE, scale=TRUE)

  pca$pops <- rownames(freqMat)

  return(pca)
}
