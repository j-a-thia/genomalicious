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
  # Code
  # --------------------------------------------+

  # Rename columns
  dat <- dat %>%
    copy %>%
    setnames(., c(chromCol, posCol, locusCol), c('CHROM','POS','LOCUS'))

  # Get unique loci
  D.uniq.loci <- dat[, c('CHROM','POS','LOCUS')] %>% unique()

  # Iterate over chromosomes
  D.keep.loci <- unique(D.uniq.loci$CHROM) %>%
    lapply(., function(chrom){
      # Subset the chormosome
      D.chr.sub.loci <- D.uniq.loci[CHROM==chrom] %>%
        setorder(., POS)

      # Keep stepping past the smallest position measured.
      # Start with the first locus, then go from there.
      # WHILE loop will only step if stepping past the most recent
      # position will take you past the extent of the dataset.
      small.pos <- D.chr.sub.loci$POS[1]
      max.pos <- max(D.chr.sub.loci$POS)
      keep.pos <- D.chr.sub.loci$POS[1]

      check.loop <- 0
      while(check.loop==0){
        x <- D.chr.sub.loci[POS>(small.pos+stepSize)][1,]$POS
        small.pos <- x
        keep.pos <- c(keep.pos, x)
        if( (small.pos+stepSize) > max.pos ){
          check.loop <- 1
        }
      }

      # Output
      D.chr.sub.loci[POS %in% keep.pos]
    }) %>%
    # Combine all chromosomes
    do.call('rbind', .)

  # Output
  return(D.keep.loci$LOCUS)
  # ... END
}
