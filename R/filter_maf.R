#' Filter loci by minor allele frequency (MAF)
#'
#' Parses a data table of genotypes/allele frequencies and returns a list of
#' loci that conform to a desired MAF threshold.
#'
#' @param dat A data table of genotypes or allele frequencies. Ggenotypes are
#' recorded either as '/' separated alleles (0/0, 0/1 1/1), or as counts of the
#' Alt allele (0, 1, 2). If allele frequencies, can be either the Ref or Alt
#' allele, so long as it is consistent across samples, populations, loci, etc.
#' Eexpects the columns:
#' \enumerate{
#'    \item The sample ID (see param \code{sampCol}), for genotype datasets only.
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The population ID (see param \code{popCol}).
#'    \item The genotypes (see param \code{genoCol}), or the allele frequencies
#'    (see param \code{freqCol})
#' }
#'
#' @param maf Numeric: The minor allele frequency. E.g. 0.05 will filter for 5%, which will remove
#' a locus if its frequency is < 0.05 or > 0.95. Default is 0.05, and the value
#' must be <=0.5.
#'
#' @param type Character: Is \code{dat} a data table of genotypes (\code{'genos'})
#' or a data table of allele frequencies (\code{'freqs'})? Default = \code{'genos'}.
#'
#' @param method Character: The method by which MAF filtering is performed.
#' One of \code{'mean'}, or \code{'any_pop'}. Default = \code{'mean'}.
#' For \code{'mean'}, the mean MAF across populations is calculated and used to
#' assess the MAF threshold at each locus. For \code{'any_pop'}, if any population
#' has a MAF less than the threshold at a locus, then that locus will be removed.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default = \code{'SAMPLE'}. Only needed when \code{type=='genos'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param popCol Character: The column name with population information.
#' Default = \code{'POP'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default = \code{'GT'}. Only needed when \code{type=='genos'}.
#'
#' @param freqCol Character: The column name with the allele frequency information.
#' Default = \code{'freqCol'}. Only needed when \code{type=='freqs'}.
#'
#' @return Returns a character vector of locus names in \code{dat[[locusCol]]} that conform
#' to the MAF threshold (>= value of \code{maf}).
#'
#' @examples
#' # LONG TABLE OF GENOTYPES
#' data(data_Genos)
#'
#' # Filter for MAF=0.20
#' loci.genos <- filter_maf(data_Genos, maf=0.20, type='genos')
#'
#' data_Genos[LOCUS %in% dt.loci]
#'
#' # LONG TABLE OF ALLELE FREQUENCIES
#' freqs_4pops <- data_Genos %>%
#'    .[, .(FREQ=sum(GT)/(length(GT)*2)), by=c('LOCUS','POP')]
#'
#' loci.freqs <- filter_maf(freqs_4pops, maf=0.20, type='freqs')
#'
#' @export
filter_maf <- function(
    dat, maf=0.05, type='genos', method='mean',
    sampCol='SAMPLE', locusCol='LOCUS', popCol='POP', genoCol='GT', freqCol='FREQ'
  ){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- c('data.table','tidyverse')
  for(L in libs){ require(L, character.only=TRUE)}

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){
    stop("Argument `dat` isn't a data table. See ?filter_maf.")
  }

  # Check columns for data.table
  if(type=='genos'){
    if(sum(c(sampCol, locusCol, genoCol, popCol) %in% colnames(dat))!=4){
      stop("Not all specified columns (`sampCol`, `locusCol`, `genoCol`, `popCol`) are in data.table dat. See ?filter_maf.")
    }
  } else if (type=='freqs'){
    if(sum(c(locusCol, freqCol, popCol) %in% colnames(dat))!=3){
      stop("Not all specified columns (`locusCol`, `freqCol`, `popCol`) are in data.table dat. See ?filter_maf.")
    }
  }

  # Check that the MAF is <=0.5
  if(maf>0.5 | maf<0){
    stop("Argument `maf` needs to be a numeric between 0 and 0.5. See ?filter_maf.")
  }

  # Check that the method has been specified correctly.
  if(!method %in% c('mean','any_pop')){
    stop("Argument `method` must be one of 'mean' or 'any_pop'. See ?filter_maf.")
  }

  # Convert genotypes to integers if needed
  if(type=='genos'){
    if(class(dat[[genoCol]])=='character'){
      dat[[genoCol]] <- genoscore_converter(dat[[genoCol]])
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Specify the min and max MAF
  minF <- maf
  maxF <- 1 - maf

  # If the input if a data.table of individuals and genotypes.
  if(type=='genos'){
    # Reassign column names
    colReass <- match(c(sampCol, locusCol, genoCol), colnames(dat))
    colnames(dat)[colReass] <- c('SAMPLE', 'LOCUS', 'GT')

    # Frequencies
    datF <- dat[, .(FREQ=sum(GT)/(length(GT)*2)), by=c('POP','LOCUS')]
  } else if(type=='freqs'){
    datF <- dat %>% copy()
  }

  # If method=='mean'
  if(method=='mean'){
    datMaf <- datF[, .(FREQ.MEAN=mean(FREQ)), by=LOCUS] %>%
      .[, MAF:=if_else(FREQ.MEAN > 0.5, 1-FREQ.MEAN, FREQ.MEAN)]

    good.loci <- datMaf[MAF>=maf]$LOCUS
  }

  # If method=='any_pop'
  if(method=='any_pop'){
    num.pops <- datF$POP %>% unique() %>% length()

    good.loci <- datF[FREQ>=minF & FREQ<=maxF] %>%
      .[, .(NUM.POPS=length(unique(POP))), by=LOCUS] %>%
      .[NUM.POPS==num.pops] %>%
      .[['LOCUS']]
  }

  # Output
  return(good.loci)

  # ........ END
}
