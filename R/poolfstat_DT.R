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
#' @param method Character: Either 'Anova' (default) or 'Identity'. Passed to \code{method} argument
#' in \code{poolfstat::computeFST}.
#'
#' @param snp.index List: A list of SNPs to consider. Default = \code{NA}.
#' Passed to \code{snp.index} argument in \code{poolfstat::computeFST}.
#'
#' @return Returns a list with two indices: \code{$Fst} is the calculated FST among the
#' pools using a function call of \code{poolfstat::computeFST}, whereas \code{$pooldat} is the
#' \code{poolfstat} object used to generate said FST values.
#'
#' @examples
#' # Load in the pool metadata and reads
#' data(data_PoolInfo)
#' data(data_PoolReads)
#'
#' # Subset to keep only Rep1 reads.
#' X <- data_PoolReads[grep(pattern='Rep1', x=data_PoolReads$SAMPLE)]
#'
#' # Need to add pool ID.
#' X$POOL <- unlist(lapply(strsplit(X$SAMPLE, '_'), function(X){ return(X[1]) }))
#'
#' # Calculate FST using poolfstat
#' Y <- poolfstat_DT(X, data_PoolInfo)
#'
#' # Output is a list
#' class(Y)
#'
#' # Outout from poolfstat::computeFST
#' Y$Fst
#'
#' # The pooldata class object, generated from data table of pooled reads
#' class(Y$pooldat)
#' Y$pooldat
#'
#'@export
poolfstat_DT <- function(dat, pool.info, method='Anova', snp.index=NA){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(i in c('tidyr', 'data.table', 'poolfstat')){ require(i, character.only=TRUE); rm(i)}

  if(sum(c('CHROM', 'POS', 'REF', 'ALT', 'LOCUS', 'POOL', 'AO', 'RO') %in% colnames(dat)) != 8){
    stop('Argument dat needs the columns $CHROM, $POS, $REF, $ALT, $LOCUS, $POOL, $AO, and $RO.')
  }

  if(sum(c('POOL', 'INDS') %in% colnames(pool.info)) != 2){
    stop('Argument pool.info needs the columns $POOL and $INDS.')
  }

  if(sum(unique(dat$POOL) %in% pool.info$POOL) != length(unique(dat$POOL))){
    stop('The pools in argument dat are not all present in argument pool.info.')
  }

  setorderv(dat, cols=c('LOCUS', 'POOL'))
  setorder(pool.info, POOL)

  dat$DP <- dat$AO + dat$RO

  dpMat <- spread(data=dat[, c('LOCUS', 'POOL', 'DP')], key=POOL, value=DP)
  roMat <- spread(data=dat[, c('LOCUS', 'POOL', 'RO')], key=POOL, value=RO)
  loci <- dat[, c('CHROM', 'POS', 'REF', 'ALT')]; loci <- loci[!duplicated(loci),]

  X <- new("pooldata")
  X@npools <- length(unique(dat$POOL))
  X@nsnp <- nrow(loci)
  X@refallele.readcount <- as.matrix(roMat[, !'LOCUS'], rownames=roMat$LOCUS)
  X@readcoverage <- as.matrix(dpMat[, !'LOCUS'], rownames=dpMat$LOCUS)
  X@snp.info <- as.matrix(loci)
  X@poolsizes <- pool.info[which(pool.info$POOL %in% colnames(dpMat[, !'LOCUS']))]$INDS * 2
  X@poolnames <- pool.info[which(pool.info$POOL %in% colnames(dpMat[, !'LOCUS']))]$POOL

  return(list(Fst=computeFST(X, method=method, snp.index=snp.index), pooldat=X))

  rm(X)
}
