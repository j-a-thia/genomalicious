#' Calculate FST with \code{poolfstat} from a data table of read counts
#'
#' Takes a data table of read counts and creates an object of class
#' \code{poolfstat}. The FST for the pools in the data table is calculated using
#' the function \code{poolfstat::computeFST}. Also requires pool size information.
#'
#' @param dat Data table: Contains read counts, e.g. like that been
#' produced by the function \code{vcf2DT}. Must contain all the following columns:
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
#'
#' @param pool.info Data table: Contains the sample sample sizes (number of diploids) for
#' for each unique pool listed in \code{dat$POOL}. Requires two columns:
#' \enumerate{
#'    \item \code{$POOL} The pools listed in \code{dat$POOL}.
#'    \item \code{$INDS} The number of diploid individuals for the pools.
#' }
#'
#' @return Returns a \code{pooldata} object as per the \code{poolfstat} package.
#'
#' @examples
#' # Load in the pool metadata and reads
#' data(data_PoolInfo)
#' data(data_PoolReads)
#'
#' # Pool info
#' data_PoolInfo
#'
#' # Replicate samples for each pool
#' data_PoolReads
#'
#' # Subset to keep only Rep1 reads.
#' subReads <- data_PoolReads[grep(pattern='Rep1', x=data_PoolReads$SAMPLE)]
#'
#' # Need to add pool ID: the ID used to match $POOL in data_PoolInfo
#' subReads$POOL <- unlist(lapply(strsplit(subReads$SAMPLE, '_'), function(X){ return(X[1]) }))
#'
#' # Make pooldata object
#' pooldataObj <-poolfstat_DT2pooldata(subReads, data_PoolInfo)
#'
#' # And go back to data table
#' pooldataTab <- poolfstat_pooldata2DT(pooldataObj)
#'
#' @export

poolfstat_DT2pooldata <- function(dat, pool.info){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(i in c('tidyr', 'data.table', 'poolfstat')){ require(i, character.only=TRUE); rm(i)}

  if(sum(c('CHROM', 'POS', 'REF', 'ALT', 'POOL', 'AO', 'RO') %in% colnames(dat)) != 7){
    stop('Argument dat needs the columns $CHROM, $POS, $REF, $ALT, $POOL, $AO, and $RO.')
  }

  if(sum(c('POOL', 'INDS') %in% colnames(pool.info)) != 2){
    stop('Argument pool.info needs the columns $POOL and $INDS.')
  }

  if(sum(unique(dat$POOL) %in% pool.info$POOL) != length(unique(dat$POOL))){
    stop('The pools in argument dat are not all present in argument pool.info.')
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  dat[, LOCUS:=paste(CHROM, POS, sep='_')]
  dat$DP <- dat$AO + dat$RO

  setorderv(dat, cols=c('LOCUS', 'POOL'))
  setorder(pool.info, POOL)

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
  pool.sizes <- pool.info$INDS*2
  names(pool.sizes) <- pool.info$POOL

  # Pool names
  pool.names <- pool.info$POOL
  names(pool.names) <- pool.info$POOL

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
