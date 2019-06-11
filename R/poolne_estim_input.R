#' Generate input files for \code{poolne_estim}
#'
#' This function prepares biallelic SNP data for analysis by the programme \code{poolne_estim}.
#'
#' @param dat Data table: The biallelic SNP data. Requires all of the following columns: \cr
#'              \enumerate{
#'                \item \code{$POOL} = The population pool ID. \cr
#'                \item \code{$SAMPLE} = The sample replicate ID; e.g. 1, 2, or 3, if three replicates were made. \cr
#'                \item \code{$CHROM} = The chromosome ID. \cr
#'                \item \code{$LOCUS} = The locus ID. \cr
#'                \item \code{$RO} = The number of Ref read counts at the locus. \cr
#'                \item \code{$DP} = The depth (total reads observed) at the locus.
#'              }
#'
#' @param pool.info Data table: The population pool metadata. Requires all of the following columns: \cr
#'              \enumerate{
#'                \item \code{$POOL} = The population pool ID. \cr
#'                \item \code{$INDS} = The number of diploid individuals in each pooled library.
#'              }
#'
#' @param runID Character: Appended inside filename to indicate a particular run or sim identifier.
#'
#' @param mcmc.sampling A vector MCMC sampling parameters for poolne_estim. Can be customised, but the
#' can be left alone (these are the defaults in \code{poolne_estim}). Must have 6 values: \cr
#'             \enumerate{
#'                \item The number of values to sample from the posterior distribution (e.g. 1000).
#'                \item The second Thinning rate (e.g. 25).
#'                \item The burn-in period length (e.g. 5000).
#'                \item The maximal number of pilot runs (e.g. 20) to adjust the parameters of the proposal distributions.
#'                \item The pilot run length (e.g. 500).
#'                \item A positive interger for random seed generation.
#'             }
#'
#' @return Four files are written to the working directory: \cr
#'             \enumerate{
#'               \item \code{[runID]_[poolID]_Ref.txt} = Counts of Ref alleles; loci in rows, sample replicates in columns.
#'               \item \code{[runID]_[poolID]_Dep.txt} = Depth of reads; loci in rows, sample replicates in columns.
#'               \item \code{[runID]_[poolID]_Loci.txt} = A list of loci in the order they appear in the count TXT files.
#'               \item \code{[runID]_[poolID]_Input.txt} = The params input file for \code{poolne_estim}.
#'             }
#'
#' @examples
#' # Load in the pool metadata and reads
#' data(genomalicious_PoolInfo)
#' data(genomalicious_PoolReads)
#'
#' # Have a look at the data: the samples are populations
#' # ('Pop') and replicate library preps ('Rep').
#' genomalicious_PoolReads$SAMPLE
#'
#' # You need to make sure there is a column that contains
#' # the pool ID. Split the SAMPLE column and return the first value:
#' X <- genomalicious_PoolReads
#' X$POOL <- unlist(lapply(strsplit(X$SAMPLE, '_'), function(X){ return(X[1]) }))
#'
#' # Check
#' X
#'
#' # Now make inputs
#' poolne_estim_input(dat=X, pool.info=genomalicious_PoolInfo, runID='genomalicious')
#'
#' @export
poolne_estim_input <- function(dat, pool.info, runID, mcmc.sampling=c(1000, 25, 5000, 20, 50, sample(as.integer(1000:9999),1))) {
  # BEGIN ............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table')){ require(lib, character.only=TRUE) }

  # Check class of dat.
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Check class of pool.info.
  if(!'data.table' %in% class(pool.info)){ stop("Argument pool.info isn't a data table") }

  # Test for the necessary columns in dat.
  if(sum(c("POOL", "SAMPLE", "CHROM", "LOCUS", "DP", "RO") %in% colnames(dat))!= 6){
    stop("Argument dat needs the columns $POOL, $SAMPLE, $CHROM, $LOCUS, $DP, and $RO.")
  }

  # Test for the necessary columns in pool.info.
  if(sum(c('POOL', 'INDS') %in% colnames(pool.info))!=2){
    stop("Argument pool.info needs the columns $POOL and $INDS.")
  }

  # Check length of mcmc.sampling.
  if(length(mcmc.sampling) != 6){ stop("Argument mcmc.sampling needs a length of 6 values.") }


  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Create a list split by $POOL.
  # Each item in the list is a data table for each pool.
  pool.ls <- split(dat[,c('POOL','SAMPLE','CHROM','LOCUS','RO','DP')], dat$POOL)

  # Create a list of unique pool names
  pool.x <- unique(dat$POOL)

  # Extract a Xth individual pool as a data table
  out <- lapply(as.list(pool.x), function(X){
    pool.ind <- pool.ls[[X]]

    # Make separate datables, which have $LOCUS in rows, sequencing replicates in the columns ($SAMPLE)
    # and the reference read/total read counts ($RO or $DP) in the cells.
    reads.ref <- spread(data=pool.ind[,c('CHROM','LOCUS','RO','SAMPLE')], key=SAMPLE, value=RO)
    reads.dep <- spread(data=pool.ind[,c('CHROM','LOCUS','DP','SAMPLE')], key=SAMPLE, value=DP)

    # Write reads.ref and reads.dep to file, also write a list of loci in order of appearance (rows)
    # and an input file for poolene_estim.

    # Make a bunch of params for output files.
    param.reffile <- paste(runID, X, 'Ref.txt',sep='_')
    param.depfile <- paste(runID, X, 'Dep.txt',sep='_')
    param.locifile <- paste(runID, X, 'Loci.txt',sep='_')
    param.inpfile <- paste(runID, X, 'Input.txt',sep='_')
    param.numloci <- nrow(reads.ref)
    param.seqreps <- ncol(reads.ref[,-c('CHROM','LOCUS')])
    param.numinds <- pool.info$INDS[match(X, pool.info$POOL)]

    # Write read counts and loci
    fwrite(reads.ref[,-c('CHROM','LOCUS')], param.reffile, col.names=FALSE)
    fwrite(reads.dep[,-c('CHROM','LOCUS')], param.depfile, col.names=FALSE)
    fwrite(reads.ref[,c('CHROM','LOCUS')], param.locifile, sep='\t', col.names=TRUE)

    # Create a list for the poolne_estim input file.
    input.ls <-list(c(param.numloci, param.seqreps, param.reffile, param.depfile)
                    , rep(param.numinds, param.seqreps)
                    , as.integer(mcmc.sampling)
                    , c(0.3, 0.3, 100, 1.25, 0.25, 0.40, 0.7)
    )

    # Create a blank file
    file.create(param.inpfile)

    # Write input.ls line by line
    lapply(input.ls, function(Y){
      # First convert each vector to a matrix with a single row, then append to out file.
      Y <- matrix(Y, nrow=1, ncol=length(Y))
      write(x=Y, file=param.inpfile, ncolumns=ncol(Y), append=T)
    })

    return(paste0('Written ', X))
  })

  for(i in unlist(out)){ cat(i, sep='\n')}
  # ............ END
}

