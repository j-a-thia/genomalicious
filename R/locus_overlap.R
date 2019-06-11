#' Proability of locus overlap between two groups
#'
#' This test returns a probability that sets of "significant" loci
#' identified in 2+ groups do not overlap by chance. For example,
#' if outlier tests are carried out in two different groups, the overlap
#' between these groups can be assessed to determine if the number of outliers
#' shared is greater than what might be expected from random chance.
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
#' @return The value returned is the proportion of null permutations that were greater or
#' equal to the observed value (i.e. the empirical number of significant loci shared between groups).
#'
#' @examples
#' grps <- list(A=c(1:101), B=c(20:100, 101:120), C=c(1:50, 70:120))
#'
#' sig <- list(A=c(1,8,88), B=c(88, 101), C=c(88, 102, 118, 120))
#' locus_overlap_2groups(sigA=c(1,5), sigB=5, lociA=1:100, lociB=1:100, perms=1000)
#'
#' @export
locus_overlap <- function(lociList, lociSig, perms=1000){
  # Name groups is not named
  if(is.null(names(lociList))){ names(lociList) <- paste0('group', length(lociList)) }
  if(is.null(names(lociSig))){ names(lociSig) <- paste0('group', length(lociSig)) }

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

  # p-value: Numbers of perms with shares >= observed, divided by number of perms
  p <- sum(permDat >= obsShares)/perms
  return(p)
}


