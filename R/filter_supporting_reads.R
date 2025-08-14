#' Filter samples by the minimum supporting reads for alleles.
#'
#' This function can almost be seen like a minor allele frequency or count filter
#' at the level of a the sample (instead of the whole dataset). It will mark
#' a sample as having insufficient supporting reads for the allele with lower coverage
#' if they are below a certain threshold. This might be useful, for example, when
#' using pooled allele frequencies, or when genotyped individuals are sequenced at
#' low-to-moderate coverage.
#'
#' @param dat Data.table: Contains the information of samples, loci, the total depth
#' of coverage, and the read count of the alterante allele. The reference allele read count
#' is assumed to be 1 - alternate allele read count.
#' Must contain the columns:
#' \enumerate{
#'    \item The sample ID (see param \code{sampCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The total read depth (see param \code{dpCol}).
#'    \item The alternate allele read counts (see param \code{aoCol}).
#' }
#'
#' @param sampCol Character: The column with the sample information.
#' Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param dpCol Character: The column with the total read depth information.
#' Default = \code{'DP'}.
#'
#' @param aoCol Character: The column with the alternate allele read count information.
#' Default = \code{'AO'}.
#'
#' @param suppReads Integer: The minimum number of supporting reads for the allele that
#' is least well covered by reads within a sample.
#'
#' @details Note, this sample will only evaluate sites for each there are reads
#' supporting both alleles. It will not evaluate sites that only have reads for the
#' reference alleles, or only have reads for the alternate allele.
#'
#' @returns Returns a data.table with the columns \code{$SAMPLE} and \code{$LOCUS},
#' the sample and locus information, and \code{KEEP}, a logical column with TRUE or FALSE
#' indicating whether a sample + locus observation should be kept based on uncertainty
#' in the supporting reads. Note, all samples + loci observations are returned, such that
#' they will match \code{dat}. This facilitates merging of the original data and results.
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#'
#' # Take a look at the read distribution for alternate alleles
#' hist(data_Genos$AO, xlab='Alt allele read counts', main='')
#'
#' # Let's find those sample + loci observations where there are not
#' # at least 5 reads supporting each allele
#' suppTest <- filter_supporting_reads(data_Genos, suppReads=5)
#'
#' head(suppTest)
#'
#' suppTest[KEEP==FALSE]
#'
#' # You could use this information to filter loci. For example, removing
#' # a locus if any sample does not meet the supporting read threshold for
#' # both alleles.
#' uniq_bad_loci <- unique(suppTest[KEEP==FALSE]$LOCUS)
#'
#' uniq_bad_loci
#'
#' data_Genos[!LOCUS %in% uniq_bad_loci]
#'
#' @export


filter_supporting_reads <- function(dat, sampCol='SAMPLE', locusCol='LOCUS', dpCol='DP', aoCol='AO', suppReads=3){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table); require(tidyverse)

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){
    stop("Argument `dat` isn't a data table. See ?filter_space_loci.")
  }

  # Check for correct columns
  if(sum(c(sampCol, locusCol, dpCol, aoCol) %in% colnames(dat))!=4){
    stop("Not all specified columns (`sampCol`, `locusCol`, `dpCol`, `aoCol`) are in data.table dat. See ?filter_space_loci.")
  }

  # Check that the supporting reads is a positive value
  if(suppReads<1){
    stop('Argument `stepSize` must be >=1. See ?filter_space_loci.')
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Rename columns
  dat <- dat %>%
    copy %>%
    setnames(., c(sampCol, locusCol, dpCol, aoCol), c('SAMPLE','LOCUS','DP','AO')) %>%
    .[, RO:=DP-AO]

  # Assign "zygosity"
  dat[, ZYG:=if_else(DP==RO | DP==AO, 'Hom', 'Het')]

  # Assign "keep" status based on zygosity and supporting read threshold
  dat[ZYG=='Hom', KEEP:=TRUE]
  dat[ZYG=='Het', KEEP:=if_else(RO>=suppReads & AO>=suppReads, TRUE, FALSE)]

  # Out
  return(dat[, c('SAMPLE','LOCUS','KEEP')])
}
