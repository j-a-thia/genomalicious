#' Bootstrap Weir & Cockerham's FST
#' 
#' Use a bootstrap resampling procedure to estimate the distribution
#' of values around observed FST value.
#' 
#' @param freqMat Matrix: Ref allele counts. Rows = populations,
#' columns = loci; make sure both are named. Row names used to label output FST matrix.
#'
#' @param sampMat Matrix: Number of sampled individuals. Rows = populations,
#' columns = loci.
#' 
#' @param boots Integer: The number of bootstrap replicates. Default = 100.
#' 
#' @return Returns a vector bootstrapped mean FST values (estimated across
#' bootstrapped resampled loci).
#' 
#' @references
#' Weir, Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evol. \cr
#' Weir, Hill (2002) Estimating F-statistics. Annu. Rev. Genet.
#' 
#' @export
fstWC_boot <- function(freqMat, sampMat, boots=100){
  
  fst_boot <- unlist(lapply(1:boots, function(i){
    # Number of loci
    n <- ncol(freqMat)
    # Get bootstrapped allele frequency matrix
    boot_freqs <- freqMat[,sample(1:n, size=n, replace=TRUE)]
    # The mean permuted FST
    return(fstWC_freqs(boot_freqs, sampMat[, colnames(boot_freqs)])$fst.mean)
  }))
  
  return(fst_boot)
}