#' Convert a long-format data table into polyRAD's RADdata object class.
#'
#' Used to generate the RADdata class object from the \code{polyRAD} package
#' described in Clark et al. (2019).
#'
#' @param dat Data.table: A data table of read counts for genotypes loci.
#' Expects the columns:
#' \enumerate{
#'    \item The chromosome (contig) ID (see param \code{chromCol}).
#'    \item The position information (see param \code{posCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The sample ID (see param \code{sampCol}).
#'    \item The reference allele nucleotides (see param \code{refCol}).
#'    \item The reference allele read counts (see param \code{roCol}).
#'    \item The alternate allele nucleotides (see param \code{altCol}).
#'    \item The alternate allele read counts (see param \code{aoCol}).
#'    (see param \code{freqCol})
#' }
#'
#' @param chromCol Character: The column name with the chromosome information.
#' Default = \code{'CHROM'}.
#'
#' @param posCol Character: The column name with the position information.
#' Default = \code{'POS'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default = \code{'SAMPLE'}. Only needed when \code{type=='genos'}.
#'
#' @param refCol Character: The column with the reference allele nucleotide information.
#' Default = \code{'REF'}.
#'
#' @param roCol Character: The column with the reference allele read count information.
#' Default = \code{'RO'}.
#'
#' @param altCol Character: The column with the alternate allele nucleotide information.
#' Default = \code{'ALT'}.
#'
#' @param aoCol Character: The column with the alternate allele read count information.
#' Default = \code{'AO'}.
#'
#' @param possPloidy List: A list of integers or numerics that represent the unique
#' ploidy values in the dataset. Default is \code{list(2)}, i.e., all samples have
#' a ploidy of 2. A list of \code{list(2, 4)} would represent possible ploidies of
#' 2 and 4.
#'
#' @param sampPloidy Integer/Numeric: Either a single value or a named vector of values
#' for each sample. Default is a single value, 2, i.e., all samples have a ploidy
#' of 2. A vector \code{c('Ind1'=2, 'Ind2'=4, 'Ind3'=2)}, for example, is a vector
#' of ploidies for three individuals, with ploidy values of 2, 4, and 2 for individuals 1, 2
#' and 3, respectively.
#'
#' @param contamRate Numeric: The contamination rate. Default is 0.001.
#'
#' @references Clark et al. (2019). G3. DOI: 10.1534/g3.118.200913
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#'
#' # Using a single ploidy
#' RD1 <- DT2RADdata(data_Genos, sampPloidy=2)
#'
#' # Using a vector of ploidies.
#' samps_uniq <- unique(data_Genos$SAMPLE)
#' samps_ploid <- rep(2, length(samps_uniq))
#' names(samps_ploid) <- samps_uniq
#'
#' samps_ploid
#'
#' RD2 <- DT2RADdata(data_Genos, sampPloidy=samps_ploid)
#'
#' Specifying multiple ploidies.
#' samps_ploid[20:40] <- 4
#'
#' samps_ploid %>% table
#'
#' RD3 <- DT2RADdata(data_Genos, possPloidy=list(2, 4), sampPloidy=samps_ploid)
#'
#' @export

DT2RADdata <- function(
    dat, chromCol='CHROM', posCol='POS', locusCol='LOCUS', sampCol='SAMPLE',
    refCol='REF', roCol='RO', altCol='ALT', aoCol='AO',
    possPloidy=list(2L), sampPloidy=2L, contamRate=0.001){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  library(data.table); library(polyRAD); library(tidyverse)

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){
    stop("Argument `dat` isn't a data table. See ?DT2RADdata.")
  }

  # Check columns are in dat
  if(sum(c(chromCol, posCol, locusCol, sampCol, refCol, roCol, altCol, aoCol) %in% colnames(dat))!=8){
    stop("Not all specified columns (`chromCol`, `posCol`, `locusCol`, `sampCol`, `refCol`, `roCol`, `altCol`, `aoCol`) are in data.table dat. See ?DT2RADdata.")
  }

  # Check that possible ploidies is a list
  if(!'list' %in% class(possPloidy)){
    stop('Argument `possPloidy` must be a "list" class object of integers. See ?DT2RADdata.')
  }

  # Check that sample ploidies are specified as in integer/numeric.
  if(!class(sampPloidy) %in% c('numeric','integer')){
    stop('Argument `sampPloidy` value must be an integer/numeric class. See ?DT2RADdata.')
  }

  # Check that sample ploidies are specified correctly as a named vector
  if(length(sampPloidy)>1){
    # Names present
    if(is.null(names(sampPloidy))){
      stop('Argument `sampPloidy` is specified as a vector, but sample names were not detected. See ?DT2RADdata.')
    }
    # All names match samples in dat
    uniq_samps <- dat[[sampCol]] %>% unique
    n_samps <- length(uniq_samps)
    if(sum(uniq_samps %in% names(sampPloidy)) != n_samps){
      stop('The samples in argument `dat` must all appear in the names in the argument `sampPloidy` if `sampPloidy` is specified as a vector. See ?DT2RADdata')
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Rename columns
  dat <- dat %>%
    copy %>%
    setnames(
      .,
      old=c(chromCol,posCol,locusCol,sampCol,refCol,roCol,altCol,aoCol),
      new=c('CHROM','POS','LOCUS','SAMPLE','REF','RO','ALT','AO')
    )

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

    possiblePloidies=possPloidy,

    contamRate=contamRate,

    alleleNucleotides=nucVec,

    taxaPloidy=sampPloidy
  )

  return(rad_data)
  # ........ END
}



