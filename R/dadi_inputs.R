#' Genertate dadi input from genotype or allele frequency data
#'
#' Creates an input file for the program dadi, described in Gutenkunst et al. (2009).
#' The input is biallelic genotypes or allele frequencies at SNP loci in a
#' long-format data table.
#'
#' @param dat Data table: A long-format data table of biallelic genotypes,
#' coded as '/' separated alleles ('0/0', '0/1', '1/1'), or counts
#' of the Alt alleles (0, 1, 2, respectively). Alternatively, a long-format
#' data table of allele frequencies.
#' Columns required for \strong{both} genotypes and allele frequencies:
#' \enumerate{
#'    \item The population ID (see param \code{popCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The reference allele (see param \code{refCol}).
#'    \item The alternate allele (see param \code{altCol}).
#' }
#' Columns required only for \strong{genotypes}:
#' \enumerate{
#'    \item The sample ID (see param \code{sampCol}).
#'    \item The genotypes (see param \code{genoCol}).
#' }
#' Columns required only for \strong{allele frequencies}:
#' \enumerate{
#'    \item The allele frequencies (see param \code{freqCol}).
#'    \item The number of individuals used to obtain the allele frequency
#'    estimate (see param \code{indsCol}).
#' }
#'
#' @param type Character: One of \code{'genos'} or \code{'freqs'}, to calculate
#' F-statistics from genotype or allele frequency data, respectively.
#'
#' @param sampCol Character: Sample ID. Default = \code{'SAMPLE'}.
#'
#' @param popCol Character: Population ID. Default = \code{'POP'}.
#'
#' @param locusCol Character: Locus ID. Default = \code{'LOCUS'}.
#'
#' @param refCol Character: Reference allele. Default = \code{'REF'}.
#'
#' @param altCol Character: Alternate allele. Default = \code{'ALT'}.
#'
#' @param genoCol Character: The genotype. Default = \code{'GT'}.
#'
#' @param freqCol Character: The reference allele frequency. Default = \code{'FREQ'}.
#'
#' @param indsCol Character: The number of individuals per population pool. Default = \code{'INDS'}.
#'
#' @param freqMethod Character: The method to estimate the SFS from allele
#' frequency data. Either \code{'probs'} or \code{'counts'}.
#' Default = \code{'probs'}. Only applicable when \code{type=='freqs'}.
#' See Details for parameterisation.
#'
#' @param popSub Character: The populations to subset out of \code{popCol}. Default = \code{NULL}.
#'
#' @param popLevels Character: An optional vector of the population IDs used
#' to manually specify the first and second population order. Default = \code{NULL}.
#'
#' @details Because pool-seq provides estimates of allele frequencies,
#' not direct observations of allele counts, we have to infer the SFS from
#' the allele frequencies. This is determined by the argument \code{freqMethod}.
#' \cr\cr
#' When \code{freqMethod=='counts'}, the default, the allele counts are simply rounded to the
#' nearest integer (e.g. 1.5 = 2, and 1.4 = 1), relative to the number of chromosomes.
#' The Ref allele counts are made first, then the Alt allele counts are made.
#' For instance, if 20 diploid individuals were pooled and the Ref allele frequency was 0.82,
#' from the 40 haploid chromosomes, 33 (32.8 rounded up) would be expected to contain the
#' Ref allele, whilst 7 (40 - 33) would be expected to carry the Alt allele. NOTE: if the
#' estimated number of individuals for the Ref allele is < 1 but > 0, this will always be
#' rounded to 1. This method will produce a consistent SFS, but note that extremely low
#' Ref allele frequencies will have a tendency to produce counts of 1.
#' \cr\cr
#' When \code{freqMethod=='probs'}, the allele counts are derived from a binomial draw using
#' R's \code{rbinom()} function. Again, if the Ref allele frequency from pooled diploids was
#' 0.82, then the SFS would be generated from the command call: \code{rbinom(n=1, size=40, prob=0.82)},
#' which would produce a probable number of Ref allele counts, and the Alt allele counts would
#' be 40 minus this number. This method will not produce consistently reproducible SFSs due
#' to the nature of the probabilistic draws. However, it does avoid potentially biasing
#' the SFS from rounding errors when allele frequencies are low.
#'
#' @return Returns a data table in the dadi input format.
#'
#' @references Gutenkunst et al. (2009) Inferring the joint demographic history of multiply populations
#' from multidimensional SNP frequency data. PLoS Genetics: 10, e1000695.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#' data(data_PoolFreqs)
#' data(data_PoolInfo)
#'
#' ### Make the dadi input from genotype data
#' dadi_inputs(dat=data_Genos, type='genos', popSub=c('Pop1', 'Pop2'))
#'
#' ### Make the dadi input from allele frequency data
#' colnames(data_PoolFreqs)
#'
#' # We need to add in the $INDS column to the data, data_PoolFreqs
#' newFreqData <- left_join(data_PoolFreqs, data_PoolInfo)
#' colnames(newFreqData)
#'
#' # Three
#' dadi_inputs(newFreqData, type='freqs', freqMethod='probs', )
#'
#' @export
dadi_inputs <- function(
    dat, type, sampCol='SAMPLE', popCol='POP', locusCol='LOCUS',
    refCol='REF', altCol='ALT', genoCol='GT', freqCol='FREQ', indsCol='INDS',
    freqMethod='probs', popSub=NULL, popLevels=NULL
  ){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'dplyr')){ require(lib, character.only = TRUE)}

  # Reassign names
  if(type=='genos'){
    colReass <- match(c(sampCol, popCol, locusCol, refCol, altCol, genoCol), colnames(dat))
    colnames(dat)[colReass] <- c('SAMPLE', 'POP', 'LOCUS', 'REF', 'ALT', 'GT')
  }

  if(type=='freqs'){
    colReass <- match(c(popCol, locusCol, refCol, altCol, freqCol, indsCol), colnames(dat))
    colnames(dat)[colReass] <- c('POP', 'LOCUS', 'REF', 'ALT', 'FREQ', 'INDS')
  }

  # Sub out the populations if specified
  if(is.null(popSub)==FALSE){
    dat <- dat[POP %in% popSub]
  }

  # Genotype scores
  if(type=='genos'){
    # Get the class of the genotypes
    gtClass <- class(dat$GT)

    # Check that genotypes are characters or counts
    if(!gtClass %in% c('character', 'numeric', 'integer')){
      stop("Check that genotypes are coded as '/' separated characters or as
         counts of the Alt allele.")
    }

    # Convert characters of separated alleles to counts
    if(gtClass=='character'){
      dat$GT <- genoscore_converter(dat$GT)
    }

    # Convert numeric allele counts to integers
    if(gtClass=='numeric'){
      dat$GT <- as.integer(dat$GT)
    }
  }

  # Check that the `freqMethod` argument has been assigned properly for
  # allele frequency estimates
  if((freqMethod %in% c('counts', 'probs'))==FALSE & type=='freqs'){
    stop(
    "Argument `freqMethod` must be 'probs' or 'counts' for allele frequency
    data. See ?dadi_inputs_freqs.")
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # For genotype data
  if(type=='genos'){
    r <- spread(dat[, (length(SAMPLE)*2) - sum(GT), by=c('LOCUS', 'POP', 'REF')]
                , key='POP', value='V1')
    setorder(r, 'LOCUS'); setnames(r, 'REF', 'Allele1')
    
    r <- as.data.table(r)
    
    a <- spread(dat[, sum(GT), by=c('LOCUS', 'POP', 'ALT')], key='POP', value='V1')
    setorder(a, 'LOCUS'); setnames(a, 'ALT', 'Allele2')
    
    a <- as.data.table(a)

    if(!is.null(popLevels)){
      r <- r[, c('LOCUS','Allele1',popLevels),with=FALSE]
      a <- a[, c('LOCUS','Allele2',popLevels),with=FALSE]
    }

    return(
      data.table(REF=paste0('-', r$Allele1, '-')
                 , ALT=paste0('-', a$Allele2, '-')
                 , r[, !'LOCUS']
                 , a[, !'LOCUS']
                 , LOCUS=r$LOCUS)
    )
  }

  # For allele frequencies
  if(type=='freqs'){
    # Simple REF counts
    if(freqMethod=='counts'){
      # Convert frequency into estimated counts of individuals with each allele
      dat$REF.COUNT <- apply(dat[, c('FREQ', 'INDS')], 1, function(X){
        p <- X[['FREQ']]
        inds <- X[['INDS']]
        ref.count <- p * (inds * 2)
        if(p != 0 & ref.count < 1){ ref.count <- 1
        } else if(p != 1 & ref.count < 1){ ref.count <- 1
        } else{ ref.count <- round(ref.count)
        }
        return(ref.count)
      })
    }

    # Probabilistic REF counts
    if(freqMethod=='probs'){
      # Get probabilistic counts of alleles using binomial draws
      dat$REF.COUNT <- apply(dat[, c('FREQ', 'INDS')], 1, function(X){
        rbinom(n=1, size=X['INDS']*2, prob=X['FREQ'])
      })
    }

    # Get the ALT allele counts
    dat[, ALT.COUNT:=(INDS*2)-REF.COUNT]

    # Some manipulations
    r <- spread(dat[,c('LOCUS', 'REF', 'POP', 'REF.COUNT')], key=POP, value=REF.COUNT)
    setorder(r, 'LOCUS'); setnames(r, 'REF', 'Allele1')
    
    r <- as.data.table(r)

    a <- spread(dat[,c('LOCUS', 'ALT', 'POP', 'ALT.COUNT')], key=POP, value=ALT.COUNT)
    setorder(a, 'LOCUS'); setnames(a, 'ALT', 'Allele2')
    
    a <- as.data.table(a)

    if(!is.null(popLevels)){
      r <- r[, c('LOCUS','Allele1',popLevels),with=FALSE]
      a <- a[, c('LOCUS','Allele2',popLevels),with=FALSE]
    }

    # Mash it all together
    data.table(
      REF=paste0('-', r$Allele1, '-')
      , ALT=paste0('-', a$Allele2, '-')
      , r[, !'LOCUS']
      , a[, !'LOCUS']
      , LOCUS=r$LOCUS
    ) %>%
      return()
  }

  # .......... END
}
