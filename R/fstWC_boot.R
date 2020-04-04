#' Bootstrapped Weir & Cockerham's FST
#'
#' Generate a vector of bootstrapped Weir and Cockerham's FST from observed
#' genotype or allele freuqency data.
#'
#' @param dat Matrix: Allele frequencies for populations, or biallelic genotypes
#' of individuals scored as integer counts of the Alt allele (0, 1, 2).
#' Populations or individuals are in rows, loci are in columns.
#' Column names are loci IDs and must exactly match those in argument \code{samp_size}.

#' @param input_type Character: One of two possible values: 'genos', calcualte
#' variance components from genotype matrix, or 'freqs', calculate variance
#' components from an allele frequency matrix.
#'
#' @param pop_id  Charater: A vector of population IDs. Default is NULL and
#' is only required if \code{input_type=='genos'}. Must be the same length
#' as \code{nrow(dat)} and the order of values must match the order of rows in \code{dat}.
#'
#' @param samp_size Matrix: The sample size for each locus in each population.
#' Default is \code{NULL} and is required if \code{input_type=='freqs'}.
#' Allows for different sample sizes at each locus, for example, if there
#' is missing data. Populations in rows, loci in columns. Rows must be
#' in the same order as rows in \code{dat}. Column names are loci IDs
#' and must all occur in \code{dat}.
#'
#' @param boots Integer: The number of bootstrap replicates. Default = 100.
#'
#' @return Returns a vector of bootstrapped multilocus FST values.
#'
#' @references
#' Weir, Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evol. \cr
#'
#' @export
fstWC_boot <- function(dat, input_type, pop_id=NULL, samp_size=NULL, boots=100){

  require(pbapply)

  # Number of loci
  n <- ncol(dat)

  # For each ith bootstrap...
  cat('Performing FST bootstrap calculations', '\n')
  fst_boot <- unlist(pblapply(1:boots, function(i){
    # Get indices of the loci to boostrap resample
    index <- sample(1:n, size=n, replace=TRUE)

    # Get bootstrapped variance components
    if(input_type=='genos'){
      bootVarcomps <- fstWC_varcomps(dat=dat[, index], input_type=input_type, pop_id=pop_id)
    } else if(input_type=='freqs'){
      bootVarcomps <- fstWC_varcomps(dat=dat[, index], input_type=input_type, samp_size=samp_size[, index])
    }

    # The mean permuted FST
    fst <- sum(bootVarcomps$NUMER) / sum(bootVarcomps$DENOM)
    return(fst)
  }))

  return(fst_boot)
}
