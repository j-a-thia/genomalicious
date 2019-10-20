#' Calculate Weir and Cockerham's FST from allele frequencies
#'
#' Takes a matrix of biallelic allele frequencies, and a matrix of sample sizes,
#' and calculates Weir and Cockerham's FST, i.e. theta (Weir & Cockerham, 1984).
#'
#' @param freqMat Matrix: Ref allele counts. Rows = populations,
#' columns = loci; make sure both are named. Row names used to label output FST matrix.
#'
#' @param sampMat Matrix: Number of sampled individuals. Rows = populations,
#' columns = loci.
#'
#' @param doPairs Logical: Should pairwise FSTwc be calculated (TRUE) or the
#' among population FSTwc (FALSE)?
#'
#' @param doDist Logical: Should a a distance matrix of FST be returned as well?
#' Will only run if \code{doPairs==TRUE}.
#'
#' @param perLocus Logical: Should the per locus FST be returned? Default is FALSE.
#'
#' @return A list with up to three different indices.
#' \enumerate{
#'    \item \code{$fst.mean}: contains the mean FST across loci. For among populations,
#'          this is a vector, whereas for pairwise analyses it is a data table
#'          with a \code{$POP1}, \code{$POP2}, and \code{$FST} column.
#'    \item \code{$fst.locus}: contains the locus-specific FST values
#'          if \code{perLocus==TRUE}. For amonong populations, this is a vector,
#'          but for pairwise analyses, this is a data table with population pairs.
#'    \item \code{$fst.dist}: will contain a distance matrix of FST values
#'          if \code{doPairs==TRUE} and \code{doDist==TRUE}.
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
#' fstWC_freqs(freqMat, sampMat, doPairs=FALSE)
#' fstWC_freqs(freqMat, sampMat, doPairs=TRUE)
#' fstWC_freqs(freqMat, sampMat, doPairs=TRUE, doDist=TRUE)
#'
#' @export
fstWC_freqs <- function(freqMat, sampMat, doPairs=FALSE, doDist=FALSE, perLocus=FALSE){
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

  if(doDist==TRUE & doPairs==FALSE){
    warning('Argument doDist set to TRUE, but doPairs set to FALSE: Will
    only calcualte a distance matrix for pairwise analyses.')
  }

  # Make sure sampMat and freqMat are in the same order row-wise
  sampMat <- sampMat[rownames(freqMat), ]

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  if(doPairs==FALSE){
    fst_out <- list()

    # Calculate variance per locus
    lociVar <- fstWC_varcomps(freqMat, sampMat)

    # Theta across loci
    fst_out$fst.mean <- sum(lociVar$NUMER) / sum(lociVar$DENOM)

    # Theta per locus
    if(perLocus==TRUE){
      fst.locus <- lociVar$NUMER / lociVar$DENOM
      names(fst.locus) <- lociVar$LOCUS
      fst_out$fst.locus <- fst.locus
    }

    # Return among population analyses
    return(fst_out)

  } else if(doPairs==TRUE){
    fst_out <- list()

    pairCombos <- combn(x=rownames(freqMat), m=2)

    # For the Xth pair
    pairLs <- apply(pairCombos, 2, function(X){
      # Calculate variance per locus
      lociVar <- fstWC_varcomps(freqMat[X,], sampMat[X,])
      # Theta across loci
      thetaMean <- data.table(POP1=X[1], POP2=X[2]
                    , FST=sum(lociVar$NUMER) / sum(lociVar$DENOM))

      # Theta per locus
      if(perLocus==TRUE){
        thetaLocus <- data.table(POP1=X[1], POP2=X[2]
                                 , LOCUS=lociVar$LOCUS
                                 , FST=lociVar$NUMER / lociVar$DENOM)
        return(list(mean=thetaMean, bylocus=thetaLocus))
      } else{ return(list(mean=thetaMean)) }
    })

    # Create a list of objects to return for pairwise analyses
    if(perLocus==TRUE){
      outLs <- list(fst.mean=do.call('rbind', lapply(pairLs, function(X){ X$mean}))
                    , fst.locus=do.call('rbind', lapply(pairLs, function(X){ X$bylocus})))
    } else{
      outLs <- list(fst.mean=do.call('rbind', lapply(pairLs, function(X){ X$mean})))
    }

    # If an FST distance matrix was requested:
    if(doDist==TRUE){
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
