#' Permutated Weir & Cockerham's FST
#'
#' Generate a vector of permuted Weir and Cockerham's FST from observed
#' genotype data.
#'
#' @param genoMat Matrix: A genotype matrix with samples in rows, loci in
#' columns, and Alt allele counts in cells (biallelic only: 0, 1, 2).
#'
#' @param pop_id  Charater: A vector of population IDs. Must be the same length
#' as \code{nrow(genoMat)} and the order of values must match the order
#' of rows in \code{genoMat}.
#'
#' @param perms Integer: The number of permutations to conduct. Default = 100.
#'
#' @return Returns a vector of permuted multilocus FST values.
#'
#' @references
#' Weir, Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evol. \cr
#'
#' @export
fstWC_perm <- function(genoMat, pop_id, perms=100){

  require(data.table); require(pbapply)

  # For each ith permutation...
  cat('Performing FST permutation calculations', '\n')
  fst_perm <- unlist(pblapply(1:perms, function(i){
    # Get randomly ordered indices to grab from genoMat
    index <- sample(1:nrow(genoMat), replace=FALSE)

    # Calculate the variance components on a randomly
    # shuffled genotype matrix, but keep population IDs
    # in the observed order.
    permVarcomp <- fstWC_varcomps(genoMat[index,], input_type='genos', pop_id=pop_id)

    # Get permuted multilocus FST
    fst <- sum(permVarcomp$NUMER) / sum(permVarcomp$DENOM)
    return(fst)
  }))

  return(fst_perm)
}
