#' Filter loci based on their spacing
#'
#' A function to space loci based on a particular step size.
#' For each chromosome (contig), starting from the first locus, each
#' subsequent locus is assessed with respect to whether it is further
#' than the step size from the previous locus. If a locus is not >= the
#' step size, it will be skipped, and the next locus evaluated.

#' @param dat Data.table: Contains the information on the genomic context of each
#' locus, that is, their position and which chromosome/contig they reside on. You can
#' pass this function a genotype data.table (e.g., as produced from \code{vcf2DT()}),
#' because it will subset only the unique chromosome, position, and locus information.
#' Must contain the columns:
#' \enumerate{
#'    \item The chromosome/contig ID (see param \code{chromCol}).
#'    \item The positional information (see param \code{posCol}).
#'    \item The locus ID (see param \code{locusCol}).
#' }
#'
#' @param chromCol Character: The column name with the chromosome information.
#' Default = \code{'CHROM'}.
#'
#' @param posCol Character: The column name with the position information.
#' Default = \code{'POS'}.
#'
#' @param locusCol Character: The column name with the locus name information.
#' Default = \code{'LOCUS'}.
#'
#' @param stepSize Integer: the size of steps between loci.
#'
#' @returns Returns a vector of loci IDs to be kept.
#'
#' @export

filter_space_loci <- function(dat, chromCol='CHROM', posCol='POS', locusCol='LOCUS', stepSize){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table); require(tidyverse)

  # Test that the data table is the correct class.
  if(!'data.table' %in% class(dat)){
    stop("Argument `dat` isn't a data table. See ?filter_space_loci.")
  }

  # Check for correct columns
  if(sum(c(chromCol, posCol, locusCol) %in% colnames(dat))!=3){
    stop("Not all specified columns (`chromCol`, `posCol`, `locusCol`) are in data.table dat. See ?filter_space_loci.")
  }

  # Check that the step size is a positive value
  if(stepSize<1){
    stop('Argument `stepSize` must be >=1. See ?filter_space_loci.')
  }

  # --------------------------------------------+
  # Internal function
  # --------------------------------------------+

  FUN_sample_positions <- function(pos.vec, stepSize) {
    # Sort positions so we always move left → right
    pos <- sort(pos.vec)
    n <- length(pos)

    # Handle trivial cases early
    if (n <= 1) {
      return(pos) # zero or one value, nothing to filter
    }
    if (stepSize <= 0) {
      return(pos) # no spacing constraint, return everything
    }

    # Greedy single-pass selection
    keep <- logical(n) # preallocate
    keep[1] <- TRUE # always keep the first
    acc <- 0 # distance since last kept

    for (i in 2:n) {
      acc <- acc + (pos[i] - pos[i - 1])
      if (acc >= stepSize) {
        keep[i] <- TRUE
        acc <- 0 # reset counter
      }
    }

    pos[keep]
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Rename columns
  dat <- dat %>%
    copy %>%
    setnames(., c(chromCol, posCol, locusCol), c('CHROM','POS','LOCUS'))

  # Get unique loci
  D.uniq.loci <- dat[, c('CHROM','POS','LOCUS')] %>% unique()

  # Loci to keep
  D.keep.loci <- D.uniq.loci[, .(POS=unique(POS)), by=CHROM] %>%
    .[, .(POS=FUN_sample_positions(POS, stepSize)), by=CHROM] %>%
    merge.data.table(., D.uniq.loci, by=c('CHROM', 'POS'))

  # Output
  return(D.keep.loci$LOCUS)
  # ... END
}
