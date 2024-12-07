#' Convert a \code{pooldata} object to a data table of read counts
#'
#' Takes a \code{pooldata} class object (from \code{poolfstat} package) and
#' converts into a long-format data.table.
#'
#' @param pooldata Pool data: an \code{pooldata} object, such as that obtained
#' from importing VCF with \code{poolfstat::vcf2pooldata}.
#'
#' @return Returns a long-format data table with the columns:
#' \enumerate{
#'    \item \code{$CHROM} The chromosome (contig) ID.
#'    \item \code{$POS} The variant position on the chromosome.
#'    \item \code{$REF} The reference allele.
#'    \item \code{$ALT} The alternate allele.
#'    \item \code{$LOCUS} The locus ID ([chrom]_[pos])
#'    \item \code{$POOL} The pool ID.
#'    \item \code{$DP} The total number of reads.
#'    \item \code{$AO} The number of reads supporting the alternate allele.
#'    \item \code{$RO} The number of reads supporting the reference allele.
#' }
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
#' pooldataObj <-poolfstat_DT2pooldata(data_PoolFreqs, data_PoolInfo)
#'
#' # And go back to data table
#' pooldataTab <- poolfstat_pooldata2DT(pooldataObj)
#'
#'@export

poolfstat_pooldata2DT <- function(pooldata){
  require(poolfstat); require(data.table); require(tidyverse)

  # SNP info
  snp.tab <- pooldata@snp.info %>%
    as.data.table %>%
    setnames(., new=c('CHROM', 'POS', 'REF', 'ALT')) %>%
    .[, LOCUS:=paste0(CHROM, '_', POS)]

  # Sample sizes
  samp.tab <- data.table(SAMPLE=pooldata@poolnames, N=pooldata@poolsizes)

  # Refereance allele counts as matrix
  ref.mat <- pooldata@refallele.readcount
  rownames(ref.mat) <- snp.tab$LOCUS
  colnames(ref.mat) <- pooldata@poolnames

  # Total depth as matrix
  dp.mat <- pooldata@readcoverage
  rownames(dp.mat) <- snp.tab$LOCUS
  colnames(dp.mat) <- pooldata@poolnames

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

