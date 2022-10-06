#' Genertate dadi input from pool-seq data
#'
#' Creates an input file for the program dadi, described in Gutenkunst et al. (2009).
#' This function specifically uses allele frequencies to produce inferred
#' allele counts, e.g., from pool-seq data.
#'
#' @param dat Data table: Must contain columns with the following information,
#' \enumerate{
#' \item Population pool ID
#' \item Locus ID
#' \item Reference allele
#' \item Alternate alelle
#' \item Reference allele freuqency
#' \item Number of individuals per population pool
#'             }
#' @param popCol Character: Population pool ID. Default = \code{'POOL'}
#' @param locusCol Character: Locus ID. Default = \code{'LOCUS'}
#' @param refCol Character: Reference allele. Default = \code{'REF'}
#' @param altCol Character: Alternate allele. Default = \code{'ALT'}
#' @param freqCol Character: The reference allele frequency. Default = \code{'FREQ'}.
#' @param indsCol Character: The number of individuals per population pool. Default = \code{'INDS'}.
#' @param poolSub Character: The pools to subset out of \code{popCol}. Default = \code{NULL}.
#' @param methodSFS Character: The method to estimate the SFS, either \code{'counts'} or
#' \code{'probs'}. Default = \code{'counts'}. See Details for parameterisation.
#' @param popLevels Character: An optional vector of the population pool IDs used
#' to manually specify the first and second population order. Default = \code{NULL}.
#'
#' @details Because pool-seq provides estimates of allele frequencies, not direct observations
#' of allele counts, we have to infer the SFS from the allele frequencies. This is determined
#' by the argument \code{methodSFS}.
#' \cr\cr
#' When \code{methodSFS=='counts'}, the default, the allele counts are simply rounded to the
#' nearest integer (e.g. 1.5 = 2, and 1.4 = 1), relative to the number of chromosomes.
#' The Ref allele counts are made first, then the Alt allele counts are made.
#' For instance, if 20 diploid individuals were pooled and the Ref allele frequency was 0.82,
#' from the 40 haploid chromosomes, 33 (32.8 rounded up) would be expected to contain the
#' Ref allele, whilst 7 (40 - 33) would be expected to carry the Alt allele. NOTE: if the
#' estimated number of individuals for the Ref allele is < 1 but > 0, this will always be
#' rounded to 1. This method will produce a consistent SFS, but note that extremely low
#' Ref allele frequencies will have a tendency to produce counts of 1.
#' \cr\cr
#' When \code{methodSFS=='probs'}, the allele counts are derived from a binomial draw using
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
#' data(data_PoolFreqs)
#' data_PoolFreqs
#'
#' # Default allele count estimation
#' dadi_inputs_freqs(
#'    dat=data_PoolFreqs,
#'    popCol='POOL',
#'    locusCol='LOCUS',
#'    refCol='REF',
#'    altCol='ALT',
#'    freqCol='FREQ',
#'    indsCol='INDS',
#'    poolSub=c('Pop1', 'Pop2')
#' )
#'
#' # Using probabilistic allele count estimation
#' dadi_inputs_freqs(
#'    dat=data_PoolFreqs,
#'    popCol='POOL',
#'    locusCol='LOCUS',
#'    refCol='REF',
#'    altCol='ALT',
#'    freqCol='FREQ',
#'    indsCol='INDS',
#'    poolSub=c('Pop1', 'Pop2'),
#'    methodSFS='probs')
#'
#'
#' @export
dadi_inputs_freqs <- function(
  dat
  , popCol='POOL'
  , locusCol='LOCUS'
  , refCol='REF'
  , altCol='ALT'
  , freqCol='FREQ'
  , indsCol='INDS'
  , poolSub=NULL
  , methodSFS='counts'
  , popLevels=NULL){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'dplyr', 'tidyr')){ require(lib, character.only = TRUE)}

  # Check that the `methodSFS` argument has been assigned properly.
  if((methodSFS %in% c('counts', 'probs'))==FALSE){
    stop("Argument `methodSFS` must be 'counts' or 'probs'. See ?dadi_inputs_freqs.")
  }

  # Reassign names
  colReass <- match(c(popCol, locusCol, refCol, altCol, freqCol, indsCol), colnames(dat))
  colnames(dat)[colReass] <- c('POOL', 'LOCUS', 'REF', 'ALT', 'FREQ', 'INDS')

  # Sub out the pools if specified
  if(is.null(poolSub)==FALSE){
    dat <- dat[POOL %in% poolSub]
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Simple REF counts
  if(methodSFS=='counts'){
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

    # Probabilistic REF counts
  } else if(methodSFS=='probs'){
    # Get probabilistic counts of alleles using binomial draws
    dat$REF.COUNT <- apply(dat[, c('FREQ', 'INDS')], 1, function(X){
      rbinom(n=1, size=X['INDS']*2, prob=X['FREQ'])
    })
  }

  # Get the ALT allele counts
  dat[, ALT.COUNT:=(INDS*2)-REF.COUNT]

  # Some manipulations
  r <- spread(dat[,c('LOCUS', 'REF', 'POOL', 'REF.COUNT')], key=POOL, value=REF.COUNT)
  setorder(r, 'LOCUS'); setnames(r, 'REF', 'Allele1')

  a <- spread(dat[,c('LOCUS', 'ALT', 'POOL', 'ALT.COUNT')], key=POOL, value=ALT.COUNT)
  setorder(a, 'LOCUS'); setnames(a, 'ALT', 'Allele2')

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

  # ........... END
}


