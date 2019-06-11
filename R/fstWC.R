#' Calculate Weir and Cockerham's FST
#'
#' Takes a matrix of biallelic allele frequencies, and a matric of sample sizes,
#' and calculates Weir and Cockerham's FST, i.e. theta (Weir & Cockerham, 1984)
#'
#' @param freqMat Matrix: Ref allele counts. Rows = populations,
#' columns = loci; make sure both are named. Row names used to label output FST matrix.
#'
#' @param sampMat Matrix: Number of sampled individuals. Rows = populations,
#' columns = loci.
#'
#' @param pairs Logical: Should pairwise FSTwc be calculated (TRUE) or the
#' among population FSTwc (FALSE)?
#'
#' @param dist Logical: Should a a distance matrix of FST be returned as well?
#' Will only run if \code{pairs==TRUE}.
#'
#' @return A list with two indices. \code{$fst.mean} contains the mean FST. In pairwise analyses,
#' this is a data table with population pairs. \code{$fst.locus} contains the locus-specific
#' FST values. If analysis was not pairwise, this is a vector, but for pairwise analyses,
#' this is a data table with population pairs. Additionally, if a pairwise analysis
#' was run and a distance matrix requested (\code{pairs==TRUE} and \code{dist==TRUE}),
#' the output can be found in \code{$fst.dist}.
#'
#' @references
#' Weir, Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evol. \cr
#' Weir, Hill (2002) Estimating F-statistics. Annu. Rev. Genet
#'
#' @examples
#' data(genomalicious_Freqs)
#' freqMat <- genomalicious_Freqs
#' sampMat <- matrix(rep(30, 32), nrow=4, ncol=8)
#' rownames(sampMat) <- paste0('Pop', 1:4)
#' colnames(sampMat) <- colnames(freqMat); rownames(sampMat) <- rownames(freqMat)
#'
#' fstWC(freqMat, sampMat, pairs=FALSE)
#' fstWC(freqMat, sampMat, pairs=TRUE)
#' fstWC(freqMat, sampMat, pairs=TRUE, dist=TRUE)
#'
#' @export
fstWC <- function(freqMat, sampMat, pairs=FALSE, dist=FALSE){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table)

  # Check dimensions
  if(sum(dim(freqMat) == dim(sampMat))!=2){
    stop('The dimensions of arguments freqMat and sampMat are not equivalent.')
  }

  # Check class
  if(sum(c(class(freqMat), class(sampMat))=='matrix')!=2){
    stop('Arguments freqMat and sampMat must both be matrices.')
  }

  # Check all samples and loci are present in both matrices
  if(sum(rownames(freqMat) %in% rownames(sampMat))!=nrow(sampMat)){
    stop('Make sure all row names in freqMat are also in sampMat.')
  }
  if(sum(colnames(freqMat) %in% colnames(sampMat))!=ncol(sampMat)){
    stop('Make sure all column names in freqMat are also in sampMat.')
  }

  if(dist==TRUE & pairs==FALSE){
    warning('Argument dist set to TRUE, but pairs set to FALSE: Will
    only calcualte a distance matrix for pairwise analyses.')
  }

  # Make sure sampMat and freqMat are in the same order row-wise
  sampMat <- sampMat[rownames(freqMat), ]

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  if(pairs==FALSE){
    # Calculate variance per locus
    lociVar <- fstWC_varcomps(freqMat, sampleMat)
    # Theta across loci
    fst.mean <- sum(lociVar$NUMER) / sum(lociVar$DENOM)
    # Theta per locus
    fst.locus <- lociVar$NUMER / lociVar$DENOM
    names(fst.locus) <- lociVar$LOCUS

    # Return among population analyses
    return(list(fst.mean=fst.mean, fst.locus=fst.locus))

  } else if(pairs==TRUE){
    pairCombos <- combn(x=rownames(freqMat), m=2)

    # For the Xth pair
    pairLs <- apply(pairCombos, 2, function(X){
      # Calculate variance per locus
      lociVar <- fstWC_varcomps(freqMat[X,], sampleMat[X,])
      # Theta across loci
      thetaMean <- data.table(POP1=X[1], POP2=X[2]
                    , FST=sum(lociVar$NUMER) / sum(lociVar$DENOM))
      # Theta per locus
      thetaLocus <- data.table(POP1=X[1], POP2=X[2]
                              , LOCUS=lociVar$LOCUS
                              , FST=lociVar$NUMER / lociVar$DENOM)

      # Return list
      return(list(mean=thetaMean, bylocus=thetaLocus))
    })

    # Create a list of objects to return for pairwise analyses
    outLs <- list(fst.mean=do.call('rbind', lapply(pairLs, function(X){ X$mean}))
                    , fst.locus=do.call('rbind', lapply(pairLs, function(X){ X$bylocus}))
    )

    # If an FST distance matrix was requested:
    if(dist==TRUE){
      # Empty matrix
      fst.dist <- matrix(0, nrow=nrow(freqMat), ncol=nrow(freqMat)
                         , dimnames=list(rownames(freqMat), rownames(freqMat)))
      # Get the pairwise info and fill the distance matrix
      for(j in 1:ncol(pairCombos)){
        pop1 <- pairCombos[1,j]
        pop2 <- pairCombos[2,j]
        fst.dist[pop2, pop1] <- outLs$fst.mean[POP1==pop1 & POP2==pop2]$FST
      }
      # Add to the output list
      outLs$fst.dist <- as.dist(fst.dist, diag=TRUE)
    }

    # Return pairwise analyses
    return(outLs)

  }
}
