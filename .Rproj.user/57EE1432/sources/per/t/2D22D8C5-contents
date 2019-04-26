#' Genertate DADI input from pool-seq data
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
#' @param poolCol Character: Population pool ID. Default = \code{'POOL'}
#' @param locusCol Character: Locus ID. Default = \code{'LOCUS'}
#' @param refCol Character: Reference allele. Default = \code{'REF'}
#' @param altCol Character: Alternate allele. Default = \code{'ALT'}
#' @param freqCol Character: The reference allele frequency. Default = \code{'FREQ'}.
#' @param indsCol Character: The number of individuals per population pool. Default = \code{'INDS'}.
#' @param poolSub Character: The pools to subset out of \code{poolCol}. Default = \code{NULL}.
#'
#' @return Returns a data tablein the DADI input format.
#'
#' @examples
#' data(genomaliciousPi)
#' genomaliciousPi
#'
#' dadi_inputs_pools(genomaliciousPi, poolCol='POOL', locusCol='LOCUS', refCol='REF', altCol='ALT', freqCol='P', indsCol='INDS')
#'
#'
#' @export
dadi_inputs_pools <- function(dat
                             , poolCol='POOL'
                             , locusCol='LOCUS'
                             , refCol='REF'
                             , altCol='ALT'
                             , freqCol='FREQ'
                             , indsCol='INDS'
                             , poolSub=NULL){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'dplyr')){ require(lib, character.only = TRUE)}

  # Reassign names
  setnames(dat, c(poolCol, locusCol, refCol, altCol, freqCol, indsCol)
           , c('POOL', 'LOCUS', 'REF', 'ALT', 'P', 'INDS'))

  # Sub out the pools if specified
  if(is.null(poolSub)==FALSE){
    dat <- dat[POOL %in% poolSub]
  }

  # Are there more than 3 populations in the POOL column?


  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Convert frequency into estimated counts of individuals with each allele
  dat[, REF.COUNT:=round(P*INDS)]
  dat[, ALT.COUNT:=INDS-REF.COUNT]

  # Some manipulations
  r <- spread(X[,c('LOCUS', 'REF', 'POOL', 'REF.COUNT')], key=POOL, value=REF.COUNT)
  setorder(r, 'LOCUS'); setnames(r, 'REF', 'Allele1')
  a <- spread(X[,c('LOCUS', 'ALT', 'POOL', 'ALT.COUNT')], key=POOL, value=ALT.COUNT)
  setorder(a, 'LOCUS'); setnames(a, 'ALT', 'Allele2')

  # Mash it all together
  return(data.table(REF=paste0('-', r$Allele1, '-')
                    , ALT=paste0('-', a$Allele2, '-')
                    , r[, !'LOCUS']
                    , a[, !'LOCUS']
                    , LOCUS=r$LOCUS
  ))

  # ........... END
}


