#' Filter loci by sequencing depth
#'
#' Parses a data table of genotypes/allele frequencies and returns a list of
#' loci that conform to a desired read depth threshold.
#'
#' @param dat Data table: The sequencing read information, must contain the columns:
#' \enumerate{
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The read depth (see param \code{dpCol}).
#' }
#'
#' @param minDP Integer: The minimum sequencing depth. Loci below this value are flagged
#' as 'bad' loci. Default is \code{NULL}.
#'
#' @param maxDP Integer: The maximum sequencing depth. Loci above this value are flagged
#' as 'bad' loci. Default is \code{NULL}.
#'
#' @param locusCol Character: The column name with the locus information. Default = \code{'LOCUS'}.
#'
#' @param dpCol Character: The column with read depth information. Default = \code{'DP'}.
#'
#' @return Returns a character vector of locus names in \code{dat[[locusCol]]}
#' that conform to the read depth threshold (>= value of \code{minDP} and <= value
#' of \code{maxDP}).
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#'
#' # Exclude loci with coverage < 10 reads
#' min10 <- filter_depth(data_Genos, minDP=10)
#' min10
#' data_Genos[LOCUS %in% min10]$DP %>% summary
#'
#' # Exclude loci with coverage < 10 and > 100 reads
#' min10max100 <- filter_depth(data_Genos, minDP=10, maxDP=100)
#' min10max100
#' data_Genos[LOCUS %in% min10max100]$DP %>% summary
#'
#' # Alternatively, subset data.table to only contain the
#' # bad loci with coverage < 10 reads
#' data_Genos[!(LOCUS %in% min10)]$DP %>% summary
#'
#' @export
filter_depth <- function(dat, minDP=NULL, maxDP=NULL, locusCol='LOCUS', dpCol='DP'){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- library(data.table)
  for(L in libs){ require(L, character.only=TRUE)}

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){
    stop("Argument dat isn't a data table. See ?filter_depth.")
    }

  # Test that at least one of the depth filters is specified
  if(is.null(minDP) & is.null(maxDP)){
    stop(
    'At least one of arguments `minDP` and `maxDP` need to be specified.
    See ?filter_depth.')
  }

  # Make sure both depth filter are specified correctly
  if((!is.null(minDP) & !is.null(maxDP))){
    if((minDP > maxDP)){
      stop('Argument `minDP` must be smaller than `maxDP`. See ?filter_depth.')
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Reassign colnames
  colReass <- match(c(locusCol, dpCol), colnames(dat))
  colnames(dat)[colReass] <- c('LOCUS', 'DP')

  # Table of loci
  locTab <- left_join(
    dat[, .(MIN=min(DP)), by=LOCUS],
    dat[, .(MAX=max(DP)), by=LOCUS]
  )

  # Filter
  if(!is.null(minDP) & is.null(maxDP)){
    good.loci <- locTab[MIN>=minDP]$LOCUS
  } else if(is.null(minDP) & !is.null(maxDP)){
    good.loci <- locTab[MAX<=maxDP]$LOCUS
  } else if(!is.null(minDP) & !is.null(maxDP)){
    good.loci <- locTab[MIN>=minDP & MAX<=maxDP]$LOCUS
  }

  return(good.loci)

  # ............ END
}
