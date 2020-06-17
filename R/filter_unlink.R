#' Filter for unlinked loci
#'
#' Parses a data table of sequence read information to randomly draw one locus
#' per contig.
#'
#' @param dat Data table: The sequencing read information, must contain the columns:
#' \code{CHROM} = the chromosome ID, i.e. contigs; and \code{LOCUS} = the locus ID.
#' \enumerate{
#'    \item The chromosome (or contig) ID (see param \code{chromCol}).
#'    \item The locus ID (see param \code{locusCol}).
#' }
#'
#' @param chromCol Character: The chromosome (or contig) information column. Default = \code{'CHROM'}.
#'
#' @param locusCol Character: The locus information column. Default = \code{'LOCUS'}.
#'
#' @return Returns a character vector of locus names in \code{dat$LOCUS} that are not on the same
#' contig in \code{dat$CHROM}.
#'
#' @examples
#' data(data_4pops)
#'
#' filter_unlink(data_4pops)
#'
#' @export
filter_unlink <- function(dat, chromCol='CHROM', locusCol='LOCUS'){

  # BEGIN ............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- library(data.table)
  for(L in libs){ require(L, character.only=TRUE)}

  # Make sure that dat is a data table.
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Make sure the right columns are in dat.
  if(length(which((c('CHROM', 'LOCUS') %in% colnames(dat))==FALSE)) > 0){
    stop("Argument dat needs the columns $CHROM and $LOCUS.")
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Reassign column names
  colReass <- match(c(chromCol, locusCol), colnames(dat))
  colnames(dat)[colReass] <- c('CHROM', 'LOCUS')

  # Split the data based on contigs, i.e. CHROM.
  dat.spl <- split(dat[,c('CHROM','LOCUS')], dat$CHROM)

  # For each CHROM, sample one LOCUS at random.
  keep.loci <- lapply(dat.spl, function(X){
    if (length(unique(X$LOCUS))>1){ return(sample(x=unique(X$LOCUS), size=1)) }
    if (length(unique(X$LOCUS))==1){ return(unique(X$LOCUS)) }
  }
  )

  # Subset original dataset by loci in 'keep.loci'.
  return(unlist(keep.loci))

  # ............ END
}
