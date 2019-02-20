#' Filter loci by sequencing depth
#'
#' Parses a data table of sequencing read information to determine which loci
#' conform to desired sequencing depth values.
#'
#' @param dat Data table: The sequencing read information, must contain the columns:
#' \code{$LOCUS} = the locus ID; and \code{$DP} = the depth of sequencing reads.
#'
#' @param minDP Integer: The minimum sequencing depth. Loci below this value are flagged
#' as 'bad' loci. Default = \code{10}. Must be specified; set \code{minDP = 0} if there is no
#' minimum number of reads.
#'
#' @param maxDP Integer: The maximum sequencing depth. Loci above this value are flagged
#' as 'bad' loci. Default = \code{NA}.
#'
#' @return Returns a character vector of locus names in \code{dat$LOCUS} that conform
#' to the read depth values specified, i.e. 'good' loci.
#'
#' @examples
#' data(pgposerReads)
#'
#' # Look at the sequencing coverage
#' pgposerReads
#'
#' # Exclude loci with coverage < 30 reads
#' min30 <- filter_depth(pgposerReads, minDP=30)
#' min30
#' pgposerReads[LOCUS %in% min30]
#'
#' # Exclude loci with coverage < 30 and > 200 reads
#' min30max200 <- filter_depth(pgposerReads, minDP=30, maxDP=200)
#' min30max200
#' pgposerReads[LOCUS %in% min30max200]
#'
#' # Alternatively, subset data.table to only contain the
#' # bad loci with coverage < 30 reads
#' pgposerReads[!(LOCUS %in% min30)]
#'
#' @export
filter_depth <- function(dat, minDP=10, maxDP=NA){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- library(data.table)
  for(L in libs){ require(L, character.only=TRUE)}

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Test for the necessary columns in dat.
  if(length(which((c('LOCUS', 'DP') %in% colnames(dat))==FALSE)) > 0){
    stop("Argument dat needs the columns $LOCUS and $DP.")
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
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
