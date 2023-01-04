#' Calculate heterozygosity from genotypes or allele frequencies
#'
#' This function takes in a long format data table of genotypes or allele
#' frequencies (as read counts) and calculates heterozygosity. Heterozygosity
#' is calculated as either the SNPwise heterozygosity (at each SNP, excluding
#' monomorphic sites), or as the genomic heterozygosity (the SNPwise
#' heterozygosity standardised for total sites assayed, including monomorphic
#' sites). Expects biallelic SNP loci. Please see Details for further explanations.
#'
#' @param snpData Data table: Genotypes for individuals or read counts for
#' populations in long-format.
#'
#' @param type Character: One of \code{'genos'} or \code{'freqs'}, to calculate
#' FST on genotype or allele frequency data, respectively.
#'
#' @param snpwise Logical: Should SNPwise heterozygosity be calculated?
#' Default is FALSE.
#'
#' @param chromData Data table: "Chromosome" size information for calculating
#' genomic heterozygosity. Default is NULL. Required when \code{snpwise==FALSE}.
#'
#' @param popData Data table: Population sampling information for calculating
#' genomic heterozygosity from allele frequencies. Default is NULL.
#' Required when \code{type=='freqs'}.
#'
#' @param byPop Logical: Should heterozgyosity be estimated for each
#' population? Default is TRUE. Note, if \code{byPop} is TRUE, then
#' \code{bySamp} must be FALSE. Required when \code{type=='genos'}.
#'
#' @param bySamp Logical: Should heterozygosity be estimated for each sampled
#' individual? Default is FALSE. Note, if \code{bySamp} is TRUE, then
#' \code{byPop} must be FALSE. Required when \code{type=='genos'}.
#'
#' @param chromCol Character: The column name with the "chromosome" information.
#' Default is \code{'CHROM'}. Must be in input data \code{snpData} and
#' \code{chromData}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default is \code{'LOCUS'}. Must be in input data \code{snpData}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default is \code{'GT'}. Must be in input data \code{snpData}, but only when
#' \code{type=='genos'}.
#'
#' @param dpCol Character: The column name with total read depth counts.
#' Default is \code{'DP}. Must be in input data \code{snpData}, but only when
#' \code{type=='freqs'}.
#'
#' @param aoCol Character: The column name with alternate allele read counts.
#' Default is \code{'AO'}. Must be in input data \code{snpData}, but only when
#' \code{type=='freqs'}.
#'
#' @param popCol Character: The column name with the population information.
#' Default is \code{'POP'}. Must be in input data \code{snpData} and \code{popData}.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default is \code{'SAMPLE'}. Must be in input data \code{snpData}, but only
#' when \code{type=='genos'}.
#'
#' @param widthCol Character: The column name with the "chromosome" width as
#' the number of nucleotide sites (base pairs). Default is 'WIDTH'.
#' Must be in input data \code{chromData}.
#'
#' @param indsCol Character: The column name with information on the number of
#' pooled individuals. Default is \code{'INDS'}. Must be in input data \code{popData}.
#'
#' @details The genomic heterozygosity (also known as the autosomal
#' heterozygosity) has been demonstrated to be the more accurate and robust
#' measure of heterozygosity (Schmidt et al. 2021). SNPwise heterozygosity can
#' suffer from sampling biases (filtering, missing data, sample size, etc).
#' Note, estimates of genomic heterozygosity are orders of magnitude less than
#' SNPwise heterozygosity (so do not be alarmed if you see very small values!).
#'
#' It is important that you DO NOT filter your data for unlinked loci when
#' calculating genomic heterozygosity (\code{snpwise==FALSE}). This is because
#' genomic heterozygosity requires measures of the number of polymorphic sites
#' relative to monomorphic sites. Contrastingly, if you choose to calculate
#' the SNPwise heterozygosity (\code{snpwise==TRUE}), then it is advisable that
#' you DO filter your data for unlinked loci, otherwise multiple loci
#' within the same genomic region will contribute to estimates (pseudoreplication).
#'
#' This function allows calculation of heterozygosity from genotype or allele
#' frequencies using the argument, \code{type}. When \code{type=='genos'}, the
#' columns \code{popCol}, \code{sampCol}, \code{chromCol}, \code{locusCol},
#' \code{genoCol} are required in the input data table, \code{snpData}.
#' When \code{type=='freqs'}, the columns \code{popCol}, \code{chromCol},
#' \code{locusCol}, \code{dpCol} and \code{aoCol} are required in \code{snpData}.
#' Heterozygosity from allele frequencies is calculated from the ratio of
#' alternate allele reads counts relative to the total read depth at a SNP locus.
#'
#' For \code{type=='genos'}, we have the option of calculating heterozygosity
#' at the level of the population or the level of the sampled individuals.
#' This is controlled using a pair of arguments, \code{byPop} and \code{bySamp}.
#' These logical arguments are directly opposed, if one is TRUE, the other has
#' to be FALSE. Population level heterozygosity is calculated when
#' \code{byPop==TRUE} and \code{bySamp==FALSE}. Conversely, sample level
#' heterozygosity is calculated when \code{byPop==FALSE} and \code{bySamp==TRUE}.
#' Heterozygosity is first calculated for each "chromosome" within individuals.
#' Estimates are then averaged across "chromosomes" for each individual, and
#' then averaged again across individuals for each population.
#'
#' For \code{type=='freqs'}, genomic heterozygosity (\code{snpwise==FALSE}),
#' is calculated using a modified version of the method described in
#' Ferretti et al. (2013). The SNPwise heterozygosity (\code{snpwise==TRUE}),
#' is calculated simply as the expected heterozygosity (2pq). Heterozygosity
#' is first calculated for each "chromosome", then averaged across "chromosomes",
#' for each population.
#'
#' When genomic heterozygosity is desired (snpwise==FALSE), then we need to
#' to provide the data table \code{chromData}. This data table as the columns
#' \code{chromCol} and \code{widthCol}. Each row is a "chromosome" (or RAD contig)
#' that was assayed for polymorphisms. The information in \code{chromData} is
#' used to standardise the polymorphisms in \code{snpData} relative to the size
#' of the genomic region they came from.
#'
#' When \code{type=='freqs'} and \code{snpwise==FALSE} (genomic heterozygosity
#' from allele frequencies) then an additional data table is required,
#' \code{popData}. This data table contains two columns, \code{popCol} and
#' \code{indsCol}. Each row contains information on the population ID and
#' the number of pooled individuals that went into estimating the allele frequencies.
#'
#' @return Returns a data table of heterozygosity estimates.
#'
#' @references
#' Ferretti et al. (2013) Molecular Ecology. DOI: 10.1111/mec.12522
#' Schmidt et al. (2021) Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.13659
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#' data(data_PoolFreqs)
#' data(data_PoolInfo)
#'
#' # Create width information
#' chrom500bp <- data.table(CHROM=unique(data_Genos$CHROM), WIDTH=500)
#'
#' # Calculate genomic heterozygosity for each population
#' het_calc(snpData=data_Genos, chromData=chrom500bp, type='genos')
#'
#' # Calculate genomic heterozygosity for each sample (individual)
#' het_calc(
#'    snpData=data_Genos, chromData=chrom500bp,
#'    type='genos', byPop=FALSE, bySamp=TRUE
#' )
#'
#' # Calculate SNPwise heterozygosity for each populations
#' het_calc(snpData=data_Genos, type='genos', snpwise=TRUE)
#'
#' # Calculate SNPwise heterozygosity for each sample (individual)
#' het_calc(
#'    snpData=data_Genos, type='genos', chromData=chrom500bp,
#'    snpwise=TRUE, byPop=FALSE, bySamp=TRUE
#' )
#'
#' # Calculate genomic heterozygosity for pooled population allele frequencies
#' het_calc(
#'    snpData=data_PoolFreqs, type='freqs', snpwise=FALSE,
#'    chromData=chrom500bp, popData=data_PoolInfo,
#'    chromCol='CHROM', locusCol='LOCUS', widthCol='WIDTH',
#'    dpCol='DP', aoCol='AO', popCol='POOL', indsCol='INDS'
#' )
#'
#' # Calculate SNPwise heterozygosity for pooled population allele frequencies
#' het_calc(
#'    snpData=data_PoolFreqs, type='freqs', snpwise=TRUE,
#'    chromData=chrom500bp, popData=data_PoolInfo,
#'    chromCol='CHROM', locusCol='LOCUS', widthCol='WIDTH',
#'    dpCol='DP', aoCol='AO', popCol='POOL', indsCol='INDS'
#' )
#'
#' @export
het_calc <- function(
    snpData, type, snpwise=FALSE, chromData=NULL, popData=NULL, byPop=TRUE, bySamp=FALSE,
    chromCol='CHROM', locusCol='LOCUS', genoCol='GT', dpCol='DP', aoCol='AO',
    popCol='POP', sampCol='SAMPLE', widthCol='WIDTH', indsCol='INDS'
){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'tidyverse')){ require(lib, character.only = TRUE)}

  # Make sure frequencies are not specified without population info
  if(type=='freqs' & is.null(popData)){
    stop('Arguments `type` is "freqs", but no `popData` specified. See ?het_calc.')
  }

  # Make sure the data is in data.table format
  snpData <- as.data.table(snpData)

  if(snpwise==FALSE){
    chromData <- as.data.table(chromData)
  }

  if(!is.null(popData)){
    popData <- as.data.table(popData)
  }

  # Check that type is specified correctly
  if(!type %in% c('genos','freqs')){
    stop('Argument `type` must be either "genos" or "freqs". See ?het_calc.')
  }

  # Check columns
  if(type=='genos'){
    gcol.check <- sum(c(chromCol,locusCol,genoCol,sampCol,popCol) %in% colnames(snpData))
    if(gcol.check!=5)
      stop(
        'Argument `snpData` must have columns specified in arguments `chromCol`,
        `locusCol`, `genoCol`, `popCol`, and `sampCol`. See ?het_calc.'
      )
  }

  if(type=='freqs'){
    fcol.check <- sum(c(chromCol,locusCol,aoCol,dpCol,popCol) %in% colnames(snpData))
    if(fcol.check!=5){
      stop(
        'Argument `snpData` must have columns specified in arguments `chromCol`,
        `locusCol`, `aoCol`, `dpCol`, and `popCol`. See ?het_calc.'
      )
    }
  }

  if(snpwise==FALSE){
    ccol.check <- sum(c(chromCol,widthCol) %in% colnames(chromData))
    if(ccol.check!=2){
      stop(
        'Argument `chromData` must have columns specified in arguments `chromCol`
        and `widthCol`. See ?het_calc.'
      )
    }
  }

  # Adjust columns
  if(type=='genos'){
    snpData <- snpData %>%
      copy %>%
      setnames(
        .,
        old=c(chromCol,locusCol,genoCol,popCol,sampCol),
        new=c('CHROM','LOCUS','GT','POP','SAMPLE')
      )
  }

  if(type=='freqs'){
    snpData <- snpData %>%
      copy %>%
      setnames(
        .,
        old=c(chromCol,locusCol,aoCol,dpCol,popCol),
        new=c('CHROM','LOCUS','AO','DP','POP')
      )

    popData <- popData %>%
      copy %>%
      setnames(., old=c(popCol,indsCol), new=c('POP','INDS'))
  }

  if(snpwise==FALSE){
    chromData <- chromData %>%
      copy %>%
      setnames(
        .,
        old=c(chromCol,widthCol),
        new=c('CHROM','WIDTH')
      )
  }

  # Check that logicals are correctly specified for type=='genos'
  if(is.logical(snpwise)==FALSE & type=='genos'){
    stop('Argument `snpwise` must be a logical value. See ?het_calc.')
  }

  if(is.logical(bySamp)==FALSE & type=='genos'){
    stop('Argument `bySamp` must be a logical value. See ?het_calc.')
  }

  if(is.logical(byPop)==FALSE & type=='genos'){
    stop('Argument `byPop` must be a logical value. See ?het_calc.')
  }

  # Check specification of bySamp and byPop for type=='genos'
  if(byPop==FALSE & bySamp==FALSE & type=='genos'){
    stop(
      'Arguments `byPop` and `bySamp` are both FALSE. One must be FALSE and the
    other must be TRUE. See ?het_calc.'
    )
  }

  if(byPop==TRUE & bySamp==TRUE & type=='genos'){
    stop(
      'Arguments `byPop` and `bySamp` are both TRUE. One must be FALSE and the
    other must be TRUE. See ?het_calc.'
    )
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Heterozygosity from genotypes
  if(type=='genos'){
    indHetSnp <- snpData %>%
      .[, .(HET.SNP=sum(GT==1)/length(GT)), by=c('POP','SAMPLE','CHROM')]

    # SNPwise heterozygosity
    if(snpwise==TRUE & bySamp==TRUE & byPop==FALSE){
      outHet <- indHetSnp %>%
        .[, .(HET.SNP=sum(HET.SNP==1)/length(HET.SNP)), by=c('POP','SAMPLE')]
    }

    if(snpwise==TRUE & bySamp==FALSE & byPop==TRUE){
      outHet <- indHetSnp %>%
        .[, sum(HET.SNP==1)/length(HET.SNP), by=c('POP','SAMPLE')] %>%
        .[, .(HET.SNP=mean(V1)), by=POP]
    }

    # Genomic heterozygosity
    if(snpwise==FALSE){
      indHetGenome <- indHetSnp %>%
        left_join(., chromData) %>%
        .[, .(HET.GENOME=HET.SNP/WIDTH), by=c('POP','SAMPLE','CHROM')]
    }

    if(snpwise==FALSE & bySamp==TRUE & byPop==FALSE){
      outHet <- indHetGenome %>%
        .[, .(HET.GENOME=mean(HET.GENOME)), by=c('POP','SAMPLE')]
    }

    if(snpwise==FALSE & bySamp==FALSE & byPop==TRUE){
      outHet <- indHetGenome %>%
        .[, mean(HET.GENOME), by=c('POP','SAMPLE')] %>%
        .[, .(HET.GENOME=mean(V1)), by=POP]
    }
  }

  # Heterozygosity from allele frequencies
  if(type=='freqs'){
    if(snpwise==TRUE){
      outHet <- snpData %>%
        .[, FREQ:=AO/DP] %>%
        .[, 2 * FREQ * (1-FREQ), by=c('POP','LOCUS')] %>%
        .[, .(HET.SNP=mean(V1)), by=POP]
    } else{
      outHet <- snpData %>%
        .[, NUMER:=2*(AO*(DP-AO))] %>%
        .[, DENOM:=DP*(DP-1)] %>%
        .[, FRAC:=NUMER/DENOM] %>%
        .[, .(SUM=FRAC), by=c('POP', 'CHROM')] %>%
        left_join(., chromData) %>%
        left_join(., popData) %>%
        .[, (INDS/(INDS-1)) * (1/WIDTH) * SUM, by=c('POP','CHROM')] %>%
        .[, .(HET.GENOME=mean(V1)), by=POP]
    }
  }

  # Output
  return(outHet)
}


