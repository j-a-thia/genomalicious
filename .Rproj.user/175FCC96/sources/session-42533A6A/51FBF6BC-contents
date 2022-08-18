#' Parametric bootstrap of \code{poolne_estim} allele frequencies
#'
#' Takes the results of \code{poolne_estim} (Gautier et al. 2013) and
#' performs a parametric bootstrap of Ref allele frequencies based on their
#' pi (estimated Ref allele frequency) and sd (the standard deviation).
#'
#' @param dat Data table: the \code{poolne_estim} data. For example, the
#' output from \code{genomalicious::poolne_estim_output()}. Requires 4 columns: \cr
#' \enumerate{
#'    \item \code{$POOL}, the population pool ID. \cr
#'    \item \code{$LOCUS}, the locus ID. \cr
#'    \item \code{$PI}, the estimated population frequency for the Ref allele. \cr
#'    \item \code{$SD}, the standard deviation for PI.
#' }
#'
#' @param num.sims Numeric, the number of simulations to generate. Default = \code{100}.
#'
#' @return A data table, with the following columns:
#' \enumerate{
#'    \item \code{$POOL}, the population pool ID. \cr
#'    \item \code{$LOCUS}, the locus ID. \cr
#'    \item \code{$BOOT.NUM}, the simulation number. \cr
#'    \item \code{$BOOT.PI}, the simulated pi (Ref allele frequency).
#' }
#'
#' @details The values of \code{PI} and \code{SD} in \code{dat} are used to generate
#' the alpha and beta paramters of a beta distribution, where: \cr
#' \cr
#' \code{alpha = ((1 - Mu) / Var - 1 / Mu) * Mu ^ 2} \cr
#' \code{beta = alpha * (1 / Mu - 1)} \cr
#' \cr
#' Here, values of \code{dat$PI} take on the values of \code{Mu} (the mean) and
#' \code{(dat$SD)^2} take on the values of \code{Var} (the variance). \cr
#' \cr
#' From the resulting beta distribution, \code{num.sims} values are drawn to create
#' a distribution of possible allele frequencies (for each locus) that might exist in the sampled
#' populations (given the associated mean and error estimated by \code{poolne_estim}).
#'
#' @examples
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Get the poolne estimat pi estimates
#' pi.data <- poolne_estim_output(stat='pi', datDir=genomaliciousExtData, lociDir=genomaliciousExtData)
#'
#' # Simulate potential distributions
#' pi.sims <- poolne_estim_boot_pi(pi.data, 100)
#' pi.sims
#'
#' @export
poolne_estim_boot_pi <- function(dat, num.sims=100){

  # BEGIN ...................

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'pbapply')){ require(lib, character.only=TRUE)}

  # Check the class of dat.
  if(!'data.table' %in% class(dat)){
    stop("Argument dat needs to be class 'data.table'.")
  }

  # Check that the correct columns are in dat.
  if(length(which((c('POOL', 'LOCUS', 'PI', 'SD') %in% colnames(dat))==FALSE)) > 0){
    stop("Argument dat needs the columns $POOL, $LOCUS, $PI, and $SD.")
  }

  # --------------------------------------------+
  # Internal function
  # --------------------------------------------+
  beta_est <- function(Mu, Var) {
    # Estimates beta distribution shape params from the mean and variance of
    # a sample.
    #
    # INPUTS:
    #   Mu    (numeric)   The sample mean
    #   Var   (numeric)   The sample variance
    #
    # OUTPUTS:
    #   A list: $alpha = alpha param, $beta = beta param

    # BEGIN ............
    alpha <- ((1 - Mu) / Var - 1 / Mu) * Mu ^ 2
    beta <- alpha * (1 / Mu - 1)
    return(params = list(alpha = alpha, beta = beta))
    # ............ END
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # 1. Get parameters for beta distribution per locus.
  # 2. Use parameters to simulate allele frequencies.
  freq.sims <- dat[, beta_est(PI, SD^2), by=c('POOL','LOCUS')] %>%
    .[, .(BOOT.PI=rbeta(num.sims, shape1=alpha, shape2=beta)), by=c('POOL','LOCUS')] %>%
    .[, BOOT.NUM:=1:.N, by=c('POOL','LOCUS')]

  # Return final dataset
  return(freq.sims)
  # ................... END
}


