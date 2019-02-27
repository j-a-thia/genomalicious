#' Proability of locus overlap between two groups
#'
#' This test returns a probability that sets of "significant" loci
#' identified in two groups do not overlap by chance. For example,
#' if outlier tests are carried out in two different groups, the overlap
#' between these groups can be assessed to determine if the number of outliers
#' shared is significant or could be observed due to random chance.
#'
#' @param sigA Character/Integer: A vector of loci in group A.
#'
#' @param sigB Character/Integer: A vector of loci in group B.
#'
#' @param lociA Character/Integer: A vector of loci analysed in group A.
#'
#' @param lociB Character/Integer: A vector of loci analysed in group B.
#'
#' @param perm Integer: The number of permutations to run for significance testing.
#'
#' @details The hypothesis being tested is that the number of observed loci overlapping
#' between two groups is greater than that observed due to random chance. Permutations
#' are run where loci are randomly drawn from each groups set of loci at a size equivalent
#' to that deemed statistically significant. In other words, if a test determined 2 loci out
#' of 100 were significant, a null draw would randomly sample 2 loci from the set of 100.
#'
#' @return The value returned is the proportion of null permutations that were greater or
#' equal to the observed value (i.e. the empirical number of significant loci shared between groups).
#'
#' @example
#' #' locus_overlap_2groups(sigA=c(1,5), sigB=5, lociA=1:100, lociB=1:100, perm=1000)
#'
#'@export
locus_overlap_2groups <- function(sigA, sigB, lociA, lociB, perms){
  # The observed number of shared loci between group A and B
  shareObs <- length(intersect(sigA, sigB))

  # Generate null distribution of shared loci
  shareNull <- lapply(1:perms, function(i){
    # Take loci sets from both groups and draw the same number of loci
    # as observed under selection
    randA <- sample(x=lociA, size=length(sigA), replace=FALSE)
    randB <- sample(x=lociB, size=length(sigB), replace=FALSE)
    # Compare the intersect
    return(length(intersect(randA, randB)))
  })
  shareNull <- unlist(shareNull)

  # Probability that null permutations produced equal or higher numbers
  # of overlap between two groups
  sum(shareNull >= shareObs)/perms
}


