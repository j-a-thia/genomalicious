#' Convert an allele frequency data table (long format) to \code{OutFLANK} input
#'
#' Converts a long format data table of allele frequency information into
#' a data frame fit for analysis by OutFLANK (Whitlock & Lotterhos 2015).
#'
#' @param dat Data table: The biallelic allele frequency data. Requires the
#' following columns:
#' \enumerate{
#'    \item The population ID.
#'    \item The locus ID.
#'    \item The allele frequency of the Reference allele.
#'    \item The number of individuals (sample size) used to estimate allele frequency.
#' }
#'
#' @param popCol Character: The column name with the sampled individual information. Default = \code{'POP'}.
#'
#' @param locusCol Character: The column name with the locus information. Default = \code{'LOCUS'}.
#'
#' @param freqCol Character: The column name with the allele frequency information. Default = \code{'FREQ'}.
#'
#' @param indsCol Character: The column name with the number of individuals (sample size)
#' used to estimate the allele frequency. Default = \code{'INDS'}.
#'
#' @param HoCol Character: An optional argument. The column name with the observed number of
#' heterozygotes. Affects estimation of FST, see details. Default = \code{NULL}.
#'
#' @details FST is estimated using modified versions of the \code{OutFLANK} functions
#' code{WC_FST_FiniteSample_Diploids_2Alleles_NoCorr()} and \code{WC_FST_FiniteSample_Diploids_2Alleles().
#' The observed heterozygosity, Ho, is required to calculate the variance components. However,
#' allele frequencies alone do not carry this information. Users can specify a column
#' containing the Ho (argument \code{HoCol}) if this is known, otherwise, a naive assumption
#' that Ho = He (the expected heterozygosity) is made.
#'
#' @return A data frame in the required format for analysis by \code{OutFLANK}.
#'
#' @references Whitlock & Lotterhos (2015) Reliable detection of loci responsible
#' for local adaptation: Inference of a null model through trimming the distribution
#' of FST. The American Naturalist 186.
#'
#' @examples
#' data('genomalicious_PoolPi')
#'
#' freqs <- genomalicious_PoolPi
#'
#' ofDF <- outflank_input_DTfreqs(freqs, 'POOL', 'LOCUS', 'PI', 'INDS')
#'
#' @export
#'
outflank_input_DTfreqs <- function(dat, popCol='POP', locusCol='LOCUS', freqCol='FREQ', indsCol='INDS', HoCol=NULL){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  library(data.table)

  if(sum(c(popCol, locusCol, freqCol, indsCol) %in% colnames(dat)) < 4){
    stop('Columns specified in arguments popCol, locusCol, freqCol
         , and indsCol must all be present in colnames(dat)')
  }

  if(is.null(HoCol)==FALSE){
    if(!HoCol %in% colnames(dat)){
      stop('Argument HoCol is specified, but is not in colnames(dat)')
      }
  }

  # Rename the population and locus column
  colnames(dat)[match(c(popCol, locusCol, freqCol, indsCol), colnames(dat))] <- c(
    'POP', 'LOCUS', 'P', 'INDS')

  if(is.null(HoCol)==FALSE){
    colnames(dat)[match(HoCol, colnames(dat))] <- c('HO')}

  # Number of populations
  numPops <- length(unique(dat$POP))

  # Split data by loci
  dat <- split(dat, dat$LOCUS)

  # Iterate through loci
  fstDF <- lapply(dat, function(locus){
    if('HO' %in% colnames(locus)){
      Ho <- locus$HO
    } else{ Ho <- NULL }

    # Corrected and uncorrected (for sample size) FST
    fst.cor <- outflank_mod_fst_correct(locus$P, locus$INDS, Ho=Ho)
    fst.nocor <- outflank_mod_fst_nocorrect(locus$P, locus$INDS, Ho=Ho)

    # An FST data frame
    fstDF <- data.frame(LocusName=locus$LOCUS[1]
                        , cbind(as.data.frame(fst.nocor[c('FSTNoCorr', 'T1NoCorr', 'T2NoCorr')])
                                , as.data.frame(fst.cor)
                        )
    )
  })

  # Combine FST data frames across loci
  fstDF <- do.call('rbind', fstDF)

  return(fstDF)
}
