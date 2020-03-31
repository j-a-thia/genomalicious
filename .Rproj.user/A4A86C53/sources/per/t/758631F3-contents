#' Calculate Weir and Cockerham's FST from genotypes
#'
#' Function takes a long format genotype data table and calculates
#' Weir & Cockerhams FST, i.e. theta (Weir & Cockerham, 1984). Can
#' also conduct permutation testing and calculate bootstrap
#' confidence intervals.
#'
#' @param dat Data table: A long format data table of biallelic genotypes,
#' coded as '/' separated alleles ('0/0', '0/1', '1/1'), or counts
#' of the Alt alleles (0, 1, 2, repsectively).
#' Four columns are required:
#' \enumerate{
#'    \item The population ID (see param \code{popCol}).
#'    \item The sampled individual ID (see param \code{sampCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The genotype (see param \code{genoCol}).
#' }
#'
#' @param popCol Character: The column name with the population information.
#' Default = \code{'POP'}.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default = \code{'GT'}.
#'
#' @param iters Integer: The number of permutations or bootstrap replicates.
#' Default = 100.
#'
#' @param permTest Logical: Should a permutation test be conducted?
#' Default = \code{TRUE}.
#'
#' @param bootCI Logical: Should bootstrap confidence intervals be estimated?
#' Default = \code{TRUE}.
#'
#' @param doPairs Logical: Should pairwise FST be calculated? Default = \code{FALSE},
#' which calculates the among population FST.
#'
#' @param doDist Logical: Should a a distance matrix of FST be returned?
#' Default = \code{FALSE}. Only applied when \code{doPairs==TRUE}.
#'
#' @param perLocus Logical: Should the per locus FST be returned?
#' Default = \code{FALSE}.
#'
#' @details Permutation tests involve random shuffly of individuals
#' among populations, recalculating FST, and testing the hypothesis
#' that the permuted FST > observed FST. The p-value represents the
#' proportion of iterations that were "TRUE" to this expression.
#' NOTE: in personal trials, randomisation of individuals among populations
#' consistently produces negative FST, hence, typically any FST > 0
#' will be "significant" in the permutation procedure,
#' based on these criteria. \cr \cr
#' Bootstrapping randomly samples the loci with replacement
#' to estimate the distribution of FST values around the across-locus
#' mean FST. Confidence intervals are calculated as the 2.5% and
#' 97.5% percentiles of the bootstrapped FST values obtained.
#' This might be a more conservative test of significance, to
#' test the hypothesis that the observed FST is > 0.
#'
#' @return Returns a list, the contents vary depending on argument choices. \cr
#' When \code{doPairs==FALSE}, the analysis is among populations.
#' \enumerate{
#'    \item \code{$fst.mean}: the mean FST across loci.
#'    \item \code{$fst.locus}: the per locus FST.
#'    \item \code{$fst.perm.pval}: the p-value from permutation tests.
#'    \item \code{$fst.boot.ci}: the bootstrap confidence intervals.
#' }
#' When \code{doPairs==TRUE}, the analysis is between population pairs.
#' \enumerate{
#'    \item \code{$fst.mean}: data table, population pairs and their statistics for mean FST.
#'    \item \code{$fst.locus}: data table, population pairs and per locus FST.
#'    \item \code{$fst.dist}: dist, a distance matrix of mean FST values.
#' }
#'
#' @examples
#' data(genomalicious_4pops)
#'
#' # Among popuation FST, with permutation tests
#' # and bootstrap confidence intervals.
#' fst_among <- fstWC_genos(genomalicious_4pops
#'                 , permTest=TRUE
#'                 , bootCI=TRUE, iters=50)
#'
#' # Pairwise FST, with per locus estimates
#' # and a distance matrix
#' fst_pairs <- fstWC_genos(genomalicious_4pops
#'                 , doPairs=TRUE
#'                 , doDist=TRUE
#'                 , perLocus=TRUE)
#'
#' @export
fstWC_genos <- function(dat, popCol='POP', sampCol='SAMPLE'
                        , locusCol='LOCUS', genoCol='GT'
                        , permTest=FALSE, bootCI=FALSE, iters=100
                        , doPairs=FALSE, doDist=FALSE, perLocus=FALSE){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  # Check the class of dat
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Check that all columns are specified correctly
  if(sum(c(popCol, sampCol, locusCol, genoCol) %in% colnames(dat)) != 4){
    stop('Arguments popCol, sampCol, locusCol, and genoCol must be
         columns in dat. See ?genos2freqs.')
  }

  # Rename columns
  colnames(dat)[
    match(c(sampCol, popCol, locusCol, genoCol), colnames(dat))
    ] <- c('SAMPLE', 'POP', 'LOCUS', 'GT')

  # Get genotype class
  genoClass <- class(dat$GT)

  # Convert to integer
  if(genoClass=='character'){
    dat[, GT:=genoscore_converter(GT)]
  } else if(genoClass=='numeric'){
    dat[, GT:=as.integer(GT)]
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Frequency matrix
  freqMat <- genos2freqs(dat, popCol='POP', sampCol='SAMPLE'
                         , locusCol='LOCUS', genoCol='GT', returnMat=TRUE)

  # Sample size matrix
  sampMat <- as.matrix(spread(dat[, length(SAMPLE), by=c('POP', 'LOCUS')]
                              , key=LOCUS, value=V1)
                       , rownames='POP')

  # If pairwise FST is not desired:
  if(doPairs==FALSE){
    # START NORMAL

    # The observed FST
    fst_obs <- fstWC_freqs(freqMat=freqMat, sampMat=sampMat[, colnames(freqMat)]
                           , perLocus=perLocus)

    if(permTest==TRUE){
      fst_perm <- fstWC_permgenos(dat, popCol='POP', sampCol='SAMPLE'
                                  , locusCol='LOCUS', genoCol='GT', perms=iters)

      # p-value for permutations
      fst_obs$fst.perm.pval <- sum(fst_perm > fst_obs$fst.mean) / iters
    }

    if(bootCI==TRUE){
      fst_boot <- fstWC_boot(freqMat, sampMat, boots=iters)

      # Percentiles
      fst_obs$fst.boot.ci <- quantile(fst_boot, c(0.025, 0.975))
    }

    return(fst_obs)
    # END NORMAL
  }

  # If pairwise FST is desired
  if(doPairs==TRUE){
    # START PAIRWISE
    fst_obs <- fstWC_freqs(freqMat=freqMat, sampMat=sampMat
                           , doPairs=TRUE, doDist=doDist, perLocus=perLocus)


    # If either permutaton testing is desired:
    if(permTest==TRUE){
      fst_obs$fst.mean$PVAL <- 0

      # Iterate through each ith population pair
      for(i in 1:nrow(fst_obs$fst.mean)){
        pop_pair <- unlist(fst_obs$fst.mean[i, c('POP1', 'POP2')])
        perm_pair <- fstWC_permgenos(dat=dat[POP %in% pop_pair]
                                     , popCol='POP', sampCol='SAMPLE'
                                     , locusCol='LOCUS', genoCol='GT', perms=iters)
        fst_obs$fst.mean$PVAL[i] <- sum(perm_pair > fst_obs$fst.mean$FST[i]) / iters
        rm(pop_pair, perm_pair, i)
      }
    }

    # If bootstrap confidence intervals are desired:
    if(bootCI==TRUE){
      fst_obs$fst.mean$CI.025 <- 0
      fst_obs$fst.mean$CI.975 <- 0

      # Iterate through each ith population pair
      for(i in 1:nrow(fst_obs$fst.mean)){
        pop_pair <- unlist(fst_obs$fst.mean[i, c('POP1', 'POP2')])
        boot_pair <- fstWC_boot(freqMat[pop_pair,], sampMat[pop_pair,], boots=iters)

        # Percentiles
        boot_perc <- quantile(boot_pair, c(0.025, 0.975))
        fst_obs$fst.mean$CI.025[i] <- boot_perc['2.5%']
        fst_obs$fst.mean$CI.975[i] <- boot_perc['97.5%']

        rm(pop_pair, boot_pair, i)
      }
    }

    return(fst_obs)
    # END PAIRWISE
  }

  # ............ END
  }
