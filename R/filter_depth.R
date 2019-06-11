#' Filter loci by sequencing depth
#'
#' Parses a data table of sequencing read information to determine which loci
#' conform to desired sequencing depth values.
#'
#' @param dat Data table: The sequencing read information, must contain the columns:
#' \enumerate{
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The read depth (see param \code{dpCol}).
#' }
#'
#' @param minDP Integer: The minimum sequencing depth. Loci below this value are flagged
#' as 'bad' loci. Default = \code{10}. Must be specified; set \code{minDP = 0} if there is no
#' minimum number of reads.
#'
#' @param maxDP Integer: The maximum sequencing depth. Loci above this value are flagged
#' as 'bad' loci. Default = \code{NA}.
#'
#' @param locusCol Character: The column name with the locus information. Default = \code{'LOCUS'}.
#'
#' @param dpCol Character: The column with read depth information. Default = \code{'DP'}.
#'
#' @return Returns a character vector of locus names in \code{dat$LOCUS} that conform
#' to the read depth values specified, i.e. 'good' loci.
#'
#' @examples
#' data(genomalicious_PoolReads)
#' datReads <- genomalicious_PoolReads
#'
#' # Look at the sequencing coverage
#' datReads
#'
#' # Exclude loci with coverage < 30 reads
#' min30 <- filter_depth(datReads, minDP=30)
#' min30
#' datReads[LOCUS %in% min30]
#'
#' # Exclude loci with coverage < 30 and > 110 reads
#' min30max110 <- filter_depth(datReads, minDP=30, maxDP=110)
#' min30max110
#' datReads[LOCUS %in% min30max110]
#'
#' # Alternatively, subset data.table to only contain the
#' # bad loci with coverage < 30 reads
#' datReads[!(LOCUS %in% min30)]
#'
#' @export
filter_depth <- function(dat, minDP=10, maxDP=NA, locusCol='LOCUS', dpCol='DP'){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- library(data.table)
  for(L in libs){ require(L, character.only=TRUE)}

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Reassign colnames
  colReass <- match(c(locusCol, dpCol), colnames(dat))
  colnames(dat)[colReass] <- c('LOCUS', 'DP')

  # Turn data into lists of data.tables, where each index is a LOCUS
  dat.spl <- split(dat, dat$LOCUS)

  # If 'maxDP' isn't specified:
  if(is.na(maxDP)){
    badLoci <- lapply(dat.spl
                      , function(L){
                        if(min(L$DP) < minDP){
                          return(data.table(LOCUS=L$LOCUS[1], KEEP='No'))
                        } else if(min(L$DP) >= minDP){
                          return(data.table(LOCUS=L$LOCUS[1], KEEP='Yes'))
                        }
                      }
    )
  }

  # If 'maxDP' is specified:
  if(!is.na(maxDP)){
    badLoci <- lapply(dat.spl
                      , function(L){
                        if(min(L$DP) >= minDP & max(L$DP) <= maxDP){
                          return(data.table(LOCUS=L$LOCUS[1], KEEP='Yes'))
                        } else {
                          return(data.table(LOCUS=L$LOCUS[1], KEEP='No'))
                        }
                      }
    )
  }

  badLoci <- do.call('rbind', badLoci)
  return(badLoci[KEEP=='Yes']$LOCUS)

  # ............ END
}
