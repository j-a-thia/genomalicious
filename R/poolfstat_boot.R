#' Bootstrap FST values from \code{poolfstat}
#'
#' Takes a data table of read counts and bootstraps the estimate of FST
#' calcualted from the function \code{poolfstat::computeFST}. Also requires
#' pool size information.
#'
#' @param dat Data table: Contains read counts, e.g. like that been
#' produced by the function \code{vcf2DT}. Must contain all the following columns:
#' \code{$CHROM}, \code{$POS}, \code{$REF}, \code{$ALT}, \code{$LOCUS}, \code{$POOL},
#' \code{$DP}, \code{$RO}.
#'
#' @param pool.info Data table: Contains the sample sample sizes (number of diploids) for
#' for each unique pool listed in \code{dat$POOL}. Requires two columns: \nr
#' \enumerate{
#'    \item \code{$POOL} The pools listed in \code{$dat$POOL}.
#'    \item \code{$INDS} The number of diploid individuals for the pools.
#' }
#'
#' @param num.sims Integer: The number of bootstrap simulations to run. Default = 100.
#'
#' @return Returns a vector of FST values produced by bootstrapping loci in the
#' original read count dataset, \code{dat}.
#'
#' @examples
#' #' # Load in the pool metadata and reads
#' data(genomaliciousInfo)
#' data(genomaliciousReads)
#'
#' # Subset to keep only Rep1 reads.
#' X <- genomaliciousReads[grep(pattern='Rep1', x=genomaliciousReads$SAMPLE)]
#'
#' # Need to add pool ID.
#' X$POOL <- unlist(lapply(strsplit(X$SAMPLE, '_'), function(X){ return(X[1]) }))
#'
#' # Bootstrap FST
#' bootFST <- boot_poolfstat(X, genomaliciousInfo, 100)
#'
#' @export
boot_poolfstat <- function(dat, pool.info, num.sims=100){
  # Get the ID and number of unique loci
  idLoci <- unique(dat$LOCUS)
  numLoci <- length(idLoci)

  # Create a vector to hold bootstrapped FST
  bsFst <- numeric()

  # Iterate for num.sims
  for(i in 1:num.sims){
    # Create the bootstrapped loci names
    bsLoci <- data.table(LOCUS.BS=paste0('Locus', 1:numLoci)
                         , LOCUS.OG=sample(x=idLoci, replace=TRUE)
    )

    # Create the bootstrapped dataset
    bsDat <- apply(bsLoci, 1, function(L){
      data.table(LOCUS=L['LOCUS.BS']
                 , dat[LOCUS==L['LOCUS.OG']
                       , c('CHROM', 'POS', 'REF', 'ALT', 'POOL', 'DP', 'RO')]
      )
    })
    bsDat <- do.call('rbind', bsDat)

    # Calculate FST from the bootstrapped dataset
    X <- poolfstat_fromDT(bsDat, pool.info)
    bsFst <- c(bsFst, X$Fst$FST)
  }

  # Return the bootstrapped FST values
  return(bsFst)
}
