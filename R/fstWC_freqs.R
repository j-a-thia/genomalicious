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
#' @param doPairs Logical: Should pairwise FST be calculated? Default = \code{FALSE},
#' which calculates the among population FST.
#'
#' @param bootCI Logical: Should bootstrap confidence intervals be estimated?
#' Default = \code{FALSE}.
#'
#' @param boots Integer: The number of bootstrap replicates. Default = 100.
#'
#' @param doDist Logical: Should a a distance matrix of FST be returned?
#' Default = \code{FALSE}. Only applied when \code{doPairs==TRUE}.
#'
#' @param perLocus Logical: Should the per locus FST be returned?
#' Default = \code{FALSE}.
#'
#' @return Returns a list, the contents vary depending on argument choices. \cr
#' When \code{doPairs==FALSE}, the analysis is among populations.
#' \enumerate{
#'    \item \code{$multilocus.mean}: the multilocus FST.
#'    \item \code{$multilocus.ci}: a vector of bootstrap confidence intervals.
#'    \item \code{$perlocus}: a data table of the FST for each locus.
#' }
#' When \code{doPairs==TRUE}, the analysis is between population pairs.
#' \enumerate{
#'    \item \code{$multilocus.mean}: a data table of multilocus FST.
#'    \item \code{$multilocus.perms}: a data table of permuted FST values.
#'    \item \code{$multilocus.pval}: a data table of p-values from permutation tests.
#'    \item \code{$multilocus.ci}: a data table of bootstrap confidence intervals.
#'    \item \code{$multilocus.dist}: a distance matrix of multilocus FST values.
#'    \item \code{$perlocus}: a data table of the FST for each locus.
#' }
#'
#' @references
#' Weir, Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evol. \cr
#' Weir, Hill (2002) Estimating F-statistics. Annu. Rev. Genet
#'
#' @examples
#' data(data_Freqs)
#' freqMat <- data_Freqs
#' sampMat <- matrix(rep(30, 32), nrow=4, ncol=8)
#' rownames(sampMat) <- paste0('Pop', 1:4)
#' colnames(sampMat) <- colnames(freqMat); rownames(sampMat) <- rownames(freqMat)
#'
#' fstWC_freqs(freqMat=freqMat, sampMat=sampMat, doPairs=FALSE)
#' fstWC_freqs(freqMat=freqMat, sampMat=sampMat, doPairs=TRUE)
#' fstWC_freqs(freqMat=freqMat, sampMat=sampMat, doPairs=TRUE, doDist=TRUE)
#'
#' @export
fstWC_freqs <- function(freqMat, sampMat, doPairs=FALSE, bootCI=FALSE, boots=100, doDist=FALSE, perLocus=FALSE){

  # --------------------------------------------+
  # Assertions and environment
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

  # Output list
  fstList <- list()

  # --------------------------------------------+
  # Code: FST among all populations
  # --------------------------------------------+
  # Calculate among populations
  if(doPairs==FALSE){
    # Calculate variance per locus
    lociVar <- fstWC_varcomps(dat=freqMat, input_type='freqs', samp_size=sampMat)

    # Theta across loci
    fstList$multilocus.mean <- sum(lociVar$NUMER) / sum(lociVar$DENOM)

    # Bootstrap theta
    if(bootCI==TRUE){
      fst.boot <- fstWC_boot(dat=freqMat, samp_size=sampMat, input_type='freqs')
      fstList$multilocus.boot <- fst.boot
      fstList$multilocus.ci <- quantile(fst.boot, c(0.025, 0.975))
    }

    # Theta per locus
    if(perLocus==TRUE){
      fst.locus <- lociVar$NUMER / lociVar$DENOM
      names(fst.locus) <- lociVar$LOCUS
      fstList$perlocus <- fst.locus
    }

    # Return among population analyses
    return(fstList)

  }

  # --------------------------------------------+
  # Code: FST between population pairs
  # --------------------------------------------+
  if(doPairs==TRUE){
    # Population pairs
    pairCombos <- combn(x=rownames(freqMat), m=2)

    # For the Xth pair
    pairList <- apply(pairCombos, 2, function(X){
      # Output for each Xth pair
      Xlist <- list()

      # Pair ID
      pair_id <- paste(X, collapse='/')

      # Calculate variance per locus
      lociVar <- fstWC_varcomps(dat=freqMat[X,], input_type='freqs', samp_size=sampMat[X,])

      # Theta across loci
      Xlist$mean <- data.table(POP1=X[1], POP2=X[2]
                    , FST=sum(lociVar$NUMER) / sum(lociVar$DENOM))

      # Bootstrap theta
      if(bootCI==TRUE){
        fst.boot <- fstWC_boot(dat=freqMat[X,], input_type='freqs', samp_size=sampMat[X,], boots=boots)
        fst.ci <- quantile(fst.boot, c(0.025, 0.975))
        Xlist$boots <- data.table(PAIR=pair_id, POP1=X[1], POP2=X[2], BOOT=1:boots, FST=fst.boot)
        Xlist$ci <- data.table(PAIR=pair_id, POP1=X[1], POP2=X[2], CI025=fst.ci[1], CI975=fst.ci[2])
      }

      # Theta per locus
      if(perLocus==TRUE){
        Xlist$perlocus <- data.table(POP1=X[1], POP2=X[2]
                                 , LOCUS=lociVar$LOCUS
                                 , FST=lociVar$NUMER / lociVar$DENOM)
      }

      # Return the data for Xth pair
      return(Xlist)
    })

    # Combine multilocus estimates
    fstList$multilocus.mean <- do.call('rbind', lapply(pairList, function(xx){ xx$mean }))

    # Combine bootstrap replicates
    if(bootCI==TRUE){
      fstList$multilocus.boot <- do.call('rbind', lapply(pairList, function(xx){ xx$boot }))
      fstList$multilocus.ci <- do.call('rbind', lapply(pairList, function(xx){ xx$ci }))
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
        fst.dist[pop2, pop1] <- fstList$multilocus.mean[POP1==pop1 & POP2==pop2]$FST
      }
      # Add to the output list
      fstList$multilocus.dist <- as.dist(fst.dist, diag=TRUE)
    }

    # Combine the per locus estimates
    if(perLocus==TRUE){
      fstList$perlocus <- do.call('rbind', lapply(pairList, function(xx){ xx$perlocus }))
    }

    # Return pairwise analyses
    return(fstList)
  }
}
