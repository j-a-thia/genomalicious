#' Genertate dadi input from individual genotype data
#'
#' Creates an input file for the program dadi, described in Gutenkunst et al. (2009).
#'
#' @param dat Data table: Contains population and genotype information. Genotypes
#' must be coded as '/' separated characters (e.g. '0/0', '0/1', or '1/1') or
#' integers of Alt allele counts (e.g. 0, 1, 2). Must contain the following columns:
#' \enumerate{
#' \item Sample ID (see argument \code{sampCol})
#' \item Population ID
#' \item Locus ID
#' \item Reference allele
#' \item Alternate alelle
#' \item Genotype
#' }
#'
#' @param sampCol Character: Sample ID. Default = \code{'SAMPLE'}.
#' @param popCol Character: Population pool ID. Default = \code{'POP'}.
#' @param locusCol Character: Locus ID. Default = \code{'LOCUS'}.
#' @param refCol Character: Reference allele. Default = \code{'REF'}.
#' @param altCol Character: Alternate allele. Default = \code{'ALT'}.
#' @param genoCol Character: The genotype. Default = \code{'GT'}.
#' @param popSub Character: The populations to subset out of \code{popCol}. Default = \code{NULL}.
#'
#' @return Returns a data table in the dadi input format.
#'
#' @references Gutenkunst et al. (2009) Inferring the joint demographic history of multiply populations
#' from multidimensional SNP frequency data. PLoS Genetics: 10, e1000695.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_4pops)
#' datGt <- copy(data_4pops)
#'
#' # Make the dadi input
#' dadi_inputs_genos(dat=datGt, popSub=c('Pop1', 'Pop2'))
#'
#' @export
dadi_inputs_genos <- function(dat
                            , sampCol='SAMPLE'
                            , popCol='POP'
                            , locusCol='LOCUS'
                            , refCol='REF'
                            , altCol='ALT'
                            , genoCol='GT'
                            , popSub=NULL){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'dplyr')){ require(lib, character.only = TRUE)}

  # Reassign names
  colReass <- match(c(sampCol, popCol, locusCol, refCol, altCol, genoCol), colnames(dat))
  colnames(dat)[colReass] <- c('SAMPLE', 'POP', 'LOCUS', 'REF', 'ALT', 'GT')

  # Sub out the pools if specified
  if(is.null(popSub)==FALSE){
    dat <- dat[POP %in% popSub]
  }

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

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  r <- spread(dat[, (length(SAMPLE)*2) - sum(GT), by=c('LOCUS', 'POP', 'REF')]
          , key='POP', value='V1')
  setorder(r, 'LOCUS'); setnames(r, 'REF', 'Allele1')
  a <- spread(dat[, sum(GT), by=c('LOCUS', 'POP', 'ALT')], key='POP', value='V1')
  setorder(a, 'LOCUS'); setnames(a, 'ALT', 'Allele2')

  return(
  data.table(REF=paste0('-', r$Allele1, '-')
             , ALT=paste0('-', a$Allele2, '-')
             , r[, !'LOCUS']
             , a[, !'LOCUS']
             , LOCUS=r$LOCUS)
  )

  # .......... END
}
