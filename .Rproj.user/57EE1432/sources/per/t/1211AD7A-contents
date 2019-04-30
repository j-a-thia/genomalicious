#' Bootstrap FST values from \code{poolfstat}
#'
#' Takes a data table of read counts and bootstraps the estimate of FST
#' calcualted from the function \code{poolfstat::computeFST}. Also requires
#' pool size information.
#'
#' @param dat Pooldata: The main data class object for \code{poolfstat}. Can be created
#' from a data table of read counts using \code{poolfstat_DT}.
#'
#' @param num.boots Integer: The number of bootstrap simulations to run. Default = 100.
#'
#' @param method Character: Either 'Anova' (default) or 'Identity'. Passed to \{method} argument
#' in \code{poolfstat::computeFST()}.
#'
#' @param snp.index List: A list of SNPs to consider. Default = \code{NA}.
#' Passed to \{snp.index} argument in \code{poolfstat::computeFST()}.
#'
#' @return Returns a vector of FST values produced by bootstrapping loci in the
#' original read count dataset, \code{dat}.
#'
#' @examples
#' #' # Load in the pool metadata and reads
#' data(genomaliciousReads)
#' data(genomaliciousInfo)
#'
#' # Subset to keep only Rep1 reads.
#' X <- genomaliciousReads[grep(pattern='Rep1', x=genomaliciousReads$SAMPLE)]
#'
#' # Need to add pool ID.
#' X$POOL <- unlist(lapply(strsplit(X$SAMPLE, '_'), function(X){ return(X[1]) }))
#'
#' # Use poolfstat_DT to compute FST for this dataset and create a pooldata object
#' poolFst <- poolfstat_DT(X, genomaliciousInfo)
#'
#' # Bootstrap FST, using pooldata object from poolFst
#' bootFst <- poolfstat_boot(poolFst$pooldat, 100)
#'
#' @export
poolfstat_boot <- function(dat, num.boots=100, method='Anova', snp.index=NA){
  # Iterate over sims
  fstBoot <- lapply(1:num.boots, function(sim){
    # Make a new copy of data to bootstrap
    bootDat <- dat
    # Rows in data matrices to resample
    bootId <- sample(1:dat@nsnp, replace=TRUE)
    # Resample ref allele counts, read coverage, and SNP info matrices
    bootDat@refallele.readcount <- bootDat@refallele.readcount[bootId, ]
    bootDat@readcoverage <- bootDat@readcoverage[bootId, ]
    bootDat@snp.info <- bootDat@snp.info[bootId, ]
    # Compute FST
    return(computeFST(bootDat, method=method, snp.index=snp.index)$FST)
  })

  fstBoot <- unlist(fstBoot)

  # Return the bootstrapped FST values
  return(fstBoot)
}
