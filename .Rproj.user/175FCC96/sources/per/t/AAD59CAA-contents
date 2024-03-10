#' Convert a long-format genotype data table into polyRAD's RADdata class
#'
#' For a data table with samples and loci as rows (long-format data table),
#' create a RADdata object as required for the R package, polyRAD (Clark et al 2019).
#'
#' @param dat Data.table: A long-format data.table with columns specified in
#' \code{chromCol}, \code{posCol}, \code{refCol}, \code{altCol}, \code{locusCol},
#' and in \code{sampCol}.
#'
#' @param chromCol Character: Column with chromosome (or contig) information. Default = \code{'CHROM'}.
#'
#' @param posCol Character: Column with position information. Default = \code{'POS'}.
#'
#' @param refCol Character: Column with reference allele information. Default = \code{'REF'}
#'
#' @param altCol Character: Column with alternate allele information. Default = \code{'ALT'}
#'
#' @param roCol Character: Column with reference observed read counts information. Default = \code{'RO'}
#'
#' @param aoCol Character: Column with alternate observed read information. Default = \code{'AO'}
#'
#' @param sampCol Character: Column with sample information. Default = \code{'SAMPLE'}.
#'
#' @references Clark et al (2019). polyRAD: Genotype calling with uncertainty from sequencing data in polyploids and diploids. G3. doi.org/10.1534/g3.118.200913
#'
#' @return Returns a RADdata object class.
#'
#' @examples
#' data(data_Genos)
#'
#' genoRADdata <- polyrad_DT2RADdata(data_Genos)
#'
#' genoRADdata
#'
#' @export
#'
polyrad_DT2RADdata <- function(
    dat, chromCol='CHROM', posCol='POS', refCol='REF', altCol='ALT',
    roCol='RO', aoCol='AO', sampCol='SAMPLE'){

  # -------------------------------------------------+
  # ASSERTIONS AND CHECKS
  # -------------------------------------------------+
  library(data.table); library(polyRAD); library(tidyverse)

  dat <- copy(dat)

  # Check that all columns are present
  check_cols <- c(chromCol, posCol, refCol, altCol, roCol, aoCol, sampCol) %in% colnames(dat)

  if(sum(check_cols)!=7){
    stop('Argument `dat` does not have all required columns. See ?polyrad_DT2RADdata')
  }

  # -------------------------------------------------+
  # BEGIN
  # -------------------------------------------------+
  # Rename columns
  dat <- dat %>%
    setnames(
      .,
      old=c(chromCol, posCol, refCol, altCol, roCol, aoCol, sampCol),
      new=c('CHROM', 'POS', 'REF', 'ALT', 'RO', 'AO', 'SAMPLE')
    )

  # Add in locus
  dat[, LOCUS:=paste0(CHROM, '_', POS)]

  # Get the unique loci
  uniqLociDt <- copy(dat[, c('LOCUS', 'CHROM', 'POS', 'REF', 'ALT')]) %>% unique()
  setorder(uniqLociDt, 'CHROM', 'POS')

  # Create allele names and number
  uniqLociDt[, VAR.REF:=paste(LOCUS, REF, sep='_')]
  uniqLociDt[, VAR.ALT:=paste(LOCUS, ALT, sep='_')]
  uniqLociDt[, VAR.NUM:=1:nrow(uniqLociDt)]

  # Get the unique samples
  uniqSamps <- sort(unique(dat$SAMPLE))

  # Create a matrix of allele read counts for each sample.
  # Each locus is represented by two columns: one for the Ref
  # allele, and another for the Alt allele counts.
  # Therefore, for n samples and m loci, the dimensions of
  # the matrix are n x 2m.
  # Create unique alleles names for each locus
  dat[, VAR.REF:=paste(LOCUS, REF, sep='_')]
  dat[, VAR.ALT:=paste(LOCUS, ALT, sep='_')]

  # Get the Ref and Alt observed read counts per allele.
  # Loci will be in the same order in each of these wide data tables.
  ro_wDT <- dcast(dat, SAMPLE ~ VAR.REF, value.var = 'RO')
  ao_wDT <- dcast(dat, SAMPLE ~ VAR.ALT, value.var = 'AO')

  # Convert to matrices
  roMat <- ro_wDT %>%
    column_to_rownames(var='SAMPLE') %>%
    as.matrix()
  aoMat <- ao_wDT %>%
    column_to_rownames(var='SAMPLE') %>%
    as.matrix()

  # Replace NA with 0
  roMat[is.na(roMat)] <- 0
  aoMat[is.na(aoMat)] <- 0

  # Make sure rownames in roMat and aoMat match
  aoMat <- aoMat[match(rownames(roMat), rownames(aoMat)),]

  # Make sure alleles are ordered by loci
  roMat <- roMat[, uniqLociDt$VAR.REF]
  aoMat <- aoMat[, uniqLociDt$VAR.ALT]

  # Combine into single matrix and transform to integer class
  depthMat <- apply(cbind(roMat, aoMat), 2, as.integer)
  rownames(depthMat) <- rownames(roMat)

  # Genetic info vectors
  chrVec <- uniqLociDt$CHROM
  posVec <- uniqLociDt$POS
  locVec <- uniqLociDt$LOCUS
  nucVec <- c(uniqLociDt$REF, uniqLociDt$ALT)

  # Build RADdata object
  rad_data <- RADdata(
    alleleDepth=depthMat,

    alleles2loc=c(1:ncol(roMat), 1:ncol(aoMat)),

    locTable=data.frame(Chr=chrVec, Pos=posVec, row.names=locVec),

    possiblePloidies=list(2L),

    contamRate=0.001,

    alleleNucleotides=nucVec,

    taxaPloidy=2
  )

  return(rad_data)
}



