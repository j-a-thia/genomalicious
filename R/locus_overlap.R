#' Proability of random null overlap among groups for outlier loci
#'
#' This functions calculates the null probability of observing
#' overlap in outlier loci between 2+ groups with respect to a scenario
#' whereby loci are drawn at random from a set of loci.
#'
#' This non-parametric sytle test is not a strict population genetic statistical
#' method per se, but provides an informal way to compare parameters of an
#' outlier analysis (number of groups, number of analysed loci,
#' and the number of observed outlier loci) to random chance.
#'
#' @param lociList List: A list of loci analysed for each group. Each index is a group,
#' containing a vector loci names.
#'
#' @param lociSig List: A list of significant loci for each group. Each index is a group,
#' containing a vector of loci names.
#'
#' @param perms Integer: The number of permutations to run for significance testing. Default = 1000
#'
#' @details
#' It is very important that the order of indices in \code{lociList} and \code{lociSig}
#' are in the same, because each index is expected to correpond to a group.
#'
#' The hypothesis being tested is that the number of observed loci overlapping
#' among groups is greater than that observed due to random chance.
#'
#' It is assumed that the amount of overlap is across all groups. E.g. if only 2 groups
#' are provided, then a locus must be significant in both groups. If 3 groups are provided,
#' then a locus must be significant in all 3 groups.
#'
#' Permutations are run where loci are randomly drawn from each group's set of loci at a size equivalent
#' to that deemed statistically significant in empircal tests. In other words, if a test
#' determined 2 loci out of 100 were significant, a null draw would randomly sample 2 loci from the set of 100.
#'
#' @return The value returned is the proportion of null permutations that were greater
#' than the observed value (i.e. the empirical number of significant loci shared between groups).
#'
#' @examples
#' grps <- list(A=c(1:100), B=c(20:100), C=c(20:100))
#'
#' sig <- list(A=c(8, 10, 16, 67, 68, 69, 88, 89, 90), B=c(22, 50, 51, 56, 57, 88, 95, 96), C=c(20, 21, 23, 44, 60, 70, 75, 88, 100))
#'
#' locus_overlap(lociList=grps, lociSig=sig, perms=1000)
#'
#' @export
locus_overlap <- function(lociList, lociSig, perms=1000){
  # Groups are not named
  if(is.null(names(lociList))){ names(lociList) <- paste0('group', 1:length(lociList)) }
  if(is.null(names(lociSig))){ names(lociSig) <- paste0('group', 1:length(lociSig)) }

  # Number of gorups
  numGrps <- length(lociSig)

  # Group names
  grpNames <- names(lociSig)

  # Counts of significant loci
  sigCount <- table(unlist(lociSig))

  # The observed number of loci shared among all groups
  obsShares <- length(sigCount[sigCount==numGrps])

  # Number of significant loci per group
  numSigs <- lapply(lociSig, length)

  # Permutations
  permDat <- lapply(1:perms, function(i){
    sharePerm <- unlist(lapply(grpNames, function(n){
      sample(x=lociList[[n]], size=numSigs[[n]], replace=FALSE)
    }))

    sharePerm <- table(sharePerm)

    return(sum(sharePerm==numGrps))
  })

  permDat <- unlist(permDat)

  # p-value: Numbers of perms with shares > observed, divided by number of perms
  p <- sum(permDat > obsShares)/perms
  return(p)
}


