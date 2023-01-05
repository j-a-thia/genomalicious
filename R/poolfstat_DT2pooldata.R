#' Convert a data table of read counts into a pooldata object
#'
#' Used to takes a data table of pool-seq read counts and creates an object of class
#' "pooldata", as per the \code{poolfstat} package (Hivert et al. 2018).
#'
#' @param dat Data table: A long-format data table of read counts for pool-seq
#' samples of biallelic SNP loci to be converted into a pooldata object.
#' Must contain all the following columns:
#' \enumerate{
#'    \item \code{$CHROM} The chromosome (contig) ID.
#'    \item \code{$POS} The variant position on the chromosome.
#'    \item \code{$REF} The reference allele.
#'    \item \code{$ALT} The alternate allele.
#'    \item \code{$LOCUS} The locus ID.
#'    \item \code{$POOL} The pool ID.
#'    \item \code{$AO} The number of reads supporting the alternate allele.
#'    \item \code{$RO} The number of reads supporting the reference allele.
#' }
#' Alternatively, a pooldata object to be converted into a data table.
#' See param \code{flip} with respect to specifying the direction of the conversion.
#'
#' @param flip Logical: Should the function be reversed? Default is \code{FALSE}.
#' When \code{flip==FALSE}, a data table is converted into a pooldata object.
#' When \code{flip==TRUE}, a pooldata object is converted into a data table.
#'
#' @param poolInfo Data table: Contains the sample sample sizes (number of diploids) for
#' for each unique pool listed in \code{dat$POOL}. Requires two columns:
#' \enumerate{
#'    \item \code{$POOL} The pools listed in \code{dat$POOL}.
#'    \item \code{$INDS} The number of diploid individuals for the pools.
#' }
#' You only need to specify \code{poolInfo} when \code{flip==FALSE}, that is,
#' transforming a data table itno a pooldata object.
#'
#' @return When \code{flip==FALSE}, returns a pooldata object as per
#' the \code{poolfstat} package. Alternatively, when \code{flip==TRUE}, returns
#' a data table.
#'
#' @references Hivert et al. (2018) Genetics. DOI: 10.1534/genetics.118.300900
#'
#' @examples
#' library(genomalicious)
#'
#' # Load in the pool metadata and reads
#' data(data_PoolInfo)
#' data(data_PoolFreqs)
#'
#' # Pool info
#' data_PoolInfo
#'
#' # Pool reads in $DP, $AO, and $RO
#' data_PoolFreqs[, c('POOL','DP','AO','RO')]
#'
#' # Make pooldata object
#' pooldataObj <-poolfstat_DT2pooldata(
#'    data_PoolFreqs,
#'    flip=FALSE,
#'    poolInfo=data_PoolInfo
#' )
#'
#' class(pooldataObj)
#'
#' pooldataObj
#'
#' # And go back to data table
#' pooldataTab <- poolfstat_DT2pooldata(pooldataObj, flip=TRUE)
#'
#' head(pooldataTab)
#'
#' @export

poolfstat_DT2pooldata <- function(dat, flip=FALSE, poolInfo){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(i in c('tidyr', 'data.table', 'poolfstat')){ require(i, character.only=TRUE); rm(i)}

  # Check that slip is specified correctly
  if(is.logical(flip)==FALSE){
    stop('Argument `flip` must be logical. See ?poolfstat_DT2pooldata.')
  }

  # Specific to flip==FALSE
  if(flip==FALSE){
    dat <- as.data.table(dat)
    dat[, POOL:=as.character(POOL)]
    poolInfo[, POOL:=as.character(POOL)]

    if(sum(c('CHROM', 'POS', 'REF', 'ALT', 'POOL', 'AO', 'RO') %in% colnames(dat)) != 7){
      stop('Argument dat needs the columns $CHROM, $POS, $REF, $ALT, $POOL, $AO, and $RO.')
    }

    if(sum(c('POOL', 'INDS') %in% colnames(poolInfo)) != 2){
      stop('Argument poolInfo needs the columns $POOL and $INDS.')
    }

    if(sum(unique(dat$POOL) %in% poolInfo$POOL) != length(unique(dat$POOL))){
      stop('The pools in argument dat are not all present in argument poolInfo.')
    }
  }

  # Specific to flip==TRUE
  if(flip==TRUE){
    if(!'pooldata' %in% class(dat)){
      stop(
        'Argument `flip` is TRUE, but argument `dat` is not a pooldata object.
        See ?poolfstat_DT2pooldata.'
      )
    }
  }

  # --------------------------------------------+
  # Code: data table to pooldata
  # --------------------------------------------+
  if(flip==FALSE){
    dat[, LOCUS:=paste(CHROM, POS, sep='_')]
    dat$DP <- dat$AO + dat$RO

    setorderv(dat, cols=c('LOCUS', 'POOL'))
    setorder(poolInfo, POOL)

    # Depth count as wide format data frame
    dpDF <- spread(data=dat[, c('LOCUS', 'POOL', 'DP')], key=POOL, value=DP) %>%
      as.data.frame() %>%
      column_to_rownames(., 'LOCUS')

    # Reference allele counts as wide format data frame
    roDF <- spread(data=dat[, c('LOCUS', 'POOL', 'RO')], key=POOL, value=RO) %>%
      as.data.frame() %>%
      column_to_rownames(., 'LOCUS')

    # Locus information
    lociDF <- dat[, c('CHROM', 'POS', 'REF', 'ALT')] %>%
      unique() %>%
      .[, LOCUS:=paste(CHROM, POS, sep='_')] %>%
      as.data.frame() %>%
      column_to_rownames(., 'LOCUS') %>%
      setnames(., new=c("Chromosome", "Position", "RefAllele", "AltAllele"))

    # Sort so all data frames are ordered the same
    roDF <- roDF[rownames(dpDF),]
    lociDF <- lociDF[rownames(dpDF),]

    # Pool sizes
    pool.sizes <- poolInfo$INDS*2
    names(pool.sizes) <- poolInfo$POOL

    # Pool names
    pool.names <- poolInfo$POOL
    names(pool.names) <- poolInfo$POOL

    # Compile pooldata object
    X <- new("pooldata")
    X@npools <- length(unique(dat$POOL))
    X@nsnp <- nrow(lociDF)
    X@refallele.readcount <- as.matrix(roDF)
    X@readcoverage <- as.matrix(dpDF)
    X@snp.info <- lociDF
    X@poolsizes <- pool.sizes[colnames(dpDF)]
    X@poolnames <- pool.names[colnames(dpDF)]

    # Output
    return(X)
  }

  # --------------------------------------------+
  # Code: pooldata to data table
  # --------------------------------------------+
  if(flip==TRUE){
    # SNP info
    snp.tab <- dat@snp.info %>%
      as.data.table %>%
      setnames(., new=c('CHROM', 'POS', 'REF', 'ALT')) %>%
      .[, LOCUS:=paste0(CHROM, '_', POS)]

    # Sample sizes
    samp.tab <- data.table(SAMPLE=dat@poolnames, N=dat@poolsizes)

    # Refereance allele counts as matrix
    ref.mat <- dat@refallele.readcount
    rownames(ref.mat) <- snp.tab$LOCUS
    colnames(ref.mat) <- dat@poolnames

    # Total depth as matrix
    dp.mat <- dat@readcoverage
    rownames(dp.mat) <- snp.tab$LOCUS
    colnames(dp.mat) <- dat@poolnames

    # Reference and depth counts as tables
    ro.tab <- ref.mat %>%
      as.data.frame %>%
      rownames_to_column(., 'LOCUS') %>%
      as.data.table %>%
      melt(., id.vars='LOCUS', variable.name='POOL', value.name='RO')

    dp.tab <- dp.mat %>%
      as.data.frame %>%
      rownames_to_column(., 'LOCUS') %>%
      as.data.table %>%
      melt(., id.vars='LOCUS', variable.name='POOL', value.name='DP')

    # Combine
    left_join(dp.tab, ro.tab) %>%
      as.data.table %>%
      .[, AO:=DP-RO] %>%
      left_join(., snp.tab) %>%
      return()
  }
}
