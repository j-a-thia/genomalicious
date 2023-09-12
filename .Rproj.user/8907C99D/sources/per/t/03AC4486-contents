#' Generate an data.table of allelic read depth
#'
#' A data.table of samples and read observations in RO (Ref) and AO (Alt) columns
#' is converted into another long-format data.table of counts for each allele.
#' This stacks Ref and Alt alleles into a single column and can be used on loci
#' with >2 alleles.
#'
#' @param dat Data.table: The input data of variants, samples, and read counts.
#' Requires the following columns:
#' \enumerate{
#'   \item The chromosome ID (see param \code{chromCol}).
#'   \item The variant position (see param \code{posCol}).
#'   \item The locus ID (see param \code{locusCol}).
#'   \item The sampled individuals (see param \code{sampCol}).
#'   \item The reference allele (see param \code{refCol}).
#'   \item The alternate allele(s) (see param \code{altCol}).
#'   \item The reference read observations (see param \code{roCol}).
#'   \item The alternate read observations (see param \code{aoCol}).
#' }
#'
#' @param chromCol The column with chromosome information. Default is \code{'CHROM'}.
#'
#' @param posCol The column with variant position information. Default is \code{'POS'}.
#'
#' @param locusCol The column with locus ID information. Default is \code{'LOCUS'}.
#'
#' @param sampCol The column with sample ID information. Default is \code{'SAMPLE'}.
#'
#' @param refCol The column with Ref allele information. Default is \code{'REF'}.
#'
#' @param altCol The column with Alt information. For loci with more than one
#' Alt allele, each allele's nucleotides should be separated by a comma. Default is \code{'ALT'}.
#'
#' @param roCol The column with Ref read observations. Default is \code{'RO'}.
#'
#' @param aoCol The column with Alt read observations. For loci with more than one
#' Alt allele, each allele's reads should be separated by a comma. Default is \code{'CHROM'}.
#'
#' @returns Returns a long format data.table with the columns:
#' \enumerate{
#'    \item \code{$CHROM}: character, the chromosome ID.
#'    \item \code{$POS}: integer, the variant position ID.
#'    \item \code{$LOCUS}: character, the locus ID.
#'    \item \code{$SAMPLE}: character, the sample ID.
#'    \item \code{$ALLELE}: character, the allele nucleotides.
#'    \item \code{$IS.REF}: logical, is the allele the Ref allele?
#'    \item \code{$READS}: integer, the observed reads for the allele.
#' }
#'
#' @examples
#' library(genomalicious)
#' data("data_Genos")
#'
#' readsPerAllele <- alleleReadsDT(data_Genos)
#' readsPerAllele
#'
#'@export

alleleReadsDT <- function(
    dat, chromCol='CHROM', posCol='POS', locusCol='LOCUS',
    sampCol='SAMPLE', refCol='REF', altCol='ALT', roCol='RO', aoCol='AO'
){
  require(data.table); require(tidyverse)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   INTERNAL FUNCTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  FUN_allele_reads <- function(loc, alle, samp, reads, is.ref){
    require(data.table)
    data.table(
      LOCUS=loc,
      SAMPLE=samp,
      ALLELE=str_split(alle, ',')[[1]],
      IS.REF=is.ref,
      READS=as.integer(str_split(reads, ',')[[1]])
    )
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   CHECKS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Check columns
  keyCols <- c(chromCol,posCol,locusCol,refCol,altCol,sampCol,roCol,aoCol)
  if(sum(keyCols %in% colnames(dat))!=8){
    stop('Argument `dat` must have all the correctly specified columns. See ?alleleCountsDT.')
  }

  # Rename columns
  dat <- setnames(
      copy(dat),
      c(chromCol,posCol,locusCol,refCol,altCol,sampCol,roCol,aoCol),
      c('CHROM','POS','LOCUS','REF','ALT','SAMPLE','RO','AO')
    )

  # Add in LOCUS column if it doesn't exist
  if(!'LOCUS' %in% colnames(dat)){
    dat[, LOCUS:=paste0(CHROM, '_', POS)]
  }

  # Add in a row index
  dat[, INDEX:=1:.N]

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   EXECUTE   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Unique variants
  uniqVars <- dat[, c('CHROM','POS','LOCUS','REF','ALT')] %>% unique

  # Allele counts
  alleReads <- rbind(
    dat[, FUN_allele_reads(loc=LOCUS, samp=SAMPLE, alle=REF, is.ref=TRUE, reads=RO), by=INDEX],
    dat[, FUN_allele_reads(loc=LOCUS, samp=SAMPLE, alle=ALT, is.ref=FALSE, reads=AO), by=INDEX]
  ) %>%
    left_join(., uniqVars[, c('LOCUS','CHROM','POS')]) %>%
    .[, c('CHROM','POS','LOCUS','SAMPLE','ALLELE','IS.REF','READS')]

  # Output
  return(alleReads)
}

