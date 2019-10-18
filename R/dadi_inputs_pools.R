#' Genertate dadi input from pool-seq data
#'
#' Creates an input file for the program dadi, described in Gutenkunst et al. (2009).
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
#' @return Returns a data table in the dadi input format. NOTE: The estimated number
#' of indivdiuals carrying an allele is rounded to nearest integer (e.g. 1.5 = 2, and 1.4 = 1),
#' with the exception when the number of individauls is < 1 but > 0, in which it is always rounded to 1.
#'
#' @references Gutenkunst et al. (2009) Inferring the joint demographic history of multiply populations
#' from multidimensional SNP frequency data. PLoS Genetics: 10, e1000695.
#'
#' @examples
#' data(genomalicious_PoolPi)
#' genomalicious_PoolPi
#'
#' dadi_inputs_pools(dat=genomalicious_PoolPi
#'                   , poolCol='POOL'
#'                   , locusCol='LOCUS'
#'                   , refCol='REF'
#'                   , altCol='ALT'
#'                   , freqCol='PI'
#'                   , indsCol='INDS'
#'                   , poolSub=c('Pop1', 'Pop2'))
#'
#'
#' @export
dadi_inputs_pools <- function(dat
                              , poolCol='POOL'
                              , locusCol='LOCUS'
                              , refCol='REF'
                              , altCol='ALT'
                              , freqCol='P'
                              , indsCol='INDS'
                              , poolSub=NULL){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'dplyr')){ require(lib, character.only = TRUE)}

  # Reassign names
  colReass <- match(c(poolCol, locusCol, refCol, altCol, freqCol, indsCol), colnames(dat))
  colnames(dat)[colReass] <- c('POOL', 'LOCUS', 'REF', 'ALT', 'P', 'INDS')

  # Sub out the pools if specified
  if(is.null(poolSub)==FALSE){
    dat <- dat[POOL %in% poolSub]
  }


  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Convert frequency into estimated counts of individuals with each allele
  dat$REF.COUNT <- apply(dat[, c('P', 'INDS')], 1, function(X){
    p <- X[['P']]
    inds <- X[['INDS']]
    ref.count <- p * (inds * 2)
    if(p != 0 & ref.count < 1){ ref.count <- 1
    } else if(p != 1 & ref.count < 1){ ref.count <- 1
    } else{ ref.count <- round(ref.count)
    }
    return(ref.count)
  })

  dat[, ALT.COUNT:=(INDS*2)-REF.COUNT]

  # Some manipulations
  r <- spread(dat[,c('LOCUS', 'REF', 'POOL', 'REF.COUNT')], key=POOL, value=REF.COUNT)
  setorder(r, 'LOCUS'); setnames(r, 'REF', 'Allele1')
  a <- spread(dat[,c('LOCUS', 'ALT', 'POOL', 'ALT.COUNT')], key=POOL, value=ALT.COUNT)
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


