#' Perform the Patterson et al. (2006) normalisation to a genotype matrix

#' @param dat Matrix: Counts of Ref allele. E.g. the C(i,j) table described
#' in Patterson et al. (2006).
#'
#' @return Returns a matrix with the same deminsions as \code{dat}, but with
#' genotypes normalised
#'
#' @export
normalise_patterson <- function(dat){
  # Iterate over each j locus, normalise genotypes
  M <- apply(dat, 2, function(j){
    # Locus j correction factor
    u <- sum(j)/length(j)
    # The underlying allele frequency
    p <- u/2
    # The normalised genotypes
    corrected_genos <- j - u
    drift_effect <- sqrt(p * (1-p))
    return(corrected_genos/drift_effect)
  })

  return(M)
}
