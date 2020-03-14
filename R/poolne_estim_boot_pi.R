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
#'    \item \code{$SIM.NUM}, the simulation number. \cr
#'    \item \code{$SIM.PI}, the simulated pi (Ref allele frequency).
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
#' pi.sims <- poolne_estim_sim_pi(pi.data, 100)
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
  # Get index for PI and SD
  index.pi <- which(names(dat)=='PI'); index.sd <- which(names(dat)=='SD')

  # Split the dataset on LOCUS
  dat.spl <- split(dat, dat$LOCUS)

  # This function works on each LOCUS
  cat('Conducting parametric bootstrap on observed frequencies', sep='\n')

  freq.sims <- pblapply(dat.spl, function(W){
    # This function work on each POOL (i.e. the rows). Subset W so it only contains
    # PI (columns = 1) and SD (column = 2)
    Y <- apply(W[,c('POOL','LOCUS','PI','SD')], 1, function(X){
      # Calculate the beta shape params (alpha and beta) from PI and SD (square SD to
      # make it into variance).
      beta.params <- beta_est(as.numeric(X['PI']), as.numeric(X['SD'])^2)
      # Generate a vector of length num.sims of pi values from a beta distribution
      # using the calculated shape params
      pi.sims <- rbeta(num.sims, beta.params$alpha, beta.params$beta)
      # Return the vector of simulated pi values.
      return(data.table(POOL=X['POOL'], LOCUS=X['LOCUS']
         , SIM.NUM=1:num.sims, SIM.PI=pi.sims))
    })

    # Return locus simulations
    return(do.call('rbind', Y))
  })

  # Return final dataset
  return(do.call('rbind', freq.sims))
  # ................... END
}
