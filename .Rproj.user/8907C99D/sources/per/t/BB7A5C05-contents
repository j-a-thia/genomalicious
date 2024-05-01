#' Generate a data table of unique alleles per locus
#'
#' Takes a data.table of variants in long-format and converts to a wide-format
#' data.table of unique alleles. Can be used for multiallelic datasets.
#'
#' @param dat Data.table: Long-format data table of variants, e.g., as read in
#' with \code{genomalicious::vcf2DT}.
#'
#' @param chromCol Character: The column with chromosome information. Default is "CHROM".
#'
#' @param posCol Character: The column with position information. Default is "POS".
#'
#' @param locusCol Character: The column with locus ID information. Default is "LOCUS".
#'
#' @param refCol Character: The column with reference allele information. Default is "REF".
#'
#' @param altCol Character: The column with alternate allele information. Default is "ALT".
#'
#' @returns Returns a wide-format data table with columns \code{$CHROM}, \code{$POS},
#' and \code{$LOCUS}, and multiple \code{$ALLELE.X} columns, where 'X' is the
#' i-th allele at a locus. The maximum number of allele columns will depend on
#' whichever locus has the most alleles. Columns with NA indicate absence of an
#' i-th allele, otherwise, the allele sequence will fill the cell. The reference
#' allele is \code{$ALLELE.0}, with subsequent columns representing alternate
#' alleles, e.g., \code{$ALLELE.1}, \code{$ALLELE.2}, etc.
#'
#' @examples
#'
#' library(genomalicious)
#'
#' data(data_Genos)
#'
#' allele_uniq_DT(data_Genos)
#'
#' @export
allele_uniq_DT <- function(
    dat, chromCol='CHROM', posCol='POS', locusCol='LOCUS',
    refCol='REF', altCol='ALT'
){
  # ----------------------------------------+
  # ASSERTIONS + REFORMAT
  # ----------------------------------------+

  require(data.table); require(tidyverse)

  # Column check
  check.cols <- c(chromCol, posCol, locusCol, refCol, altCol) %in% colnames(dat)
  if(sum(check.cols)!=5){
    stop('At least one of `chromCol`, `posCol`, `locusCol`, `refCol` or `altCol` is not in `dat`. See ?allele_DT.')
  }

  # Reassign
  dat <- dat %>% copy %>% as.data.table()
  setnames(
    dat,
    old=c(chromCol, posCol, locusCol, refCol, altCol),
    new=c('CHROM','POS','LOCUS','REF','ALT')
  )

  # ----------------------------------------+
  # MAIN CODE
  # ----------------------------------------+

  # Unique loci
  uniq.loc.tab <- dat[, c('CHROM','POS','LOCUS','REF','ALT')] %>% unique

  # Temporary allele table
  tmp.alle <- uniq.loc.tab[, tstrsplit(ALT,',')] %>%
    cbind(uniq.loc.tab[,'REF'], .)

  # Number of alleles
  colnames(tmp.alle) <- paste0('ALLELE.', 0:(ncol(tmp.alle)-1) )

  # Add back in
  uniq.loc.tab <- cbind(uniq.loc.tab, tmp.alle)

  uniq.loc.tab <- uniq.loc.tab[, !c('CHROM','POS','REF','ALT')] %>%
    melt(., id.vars = 'LOCUS', variable.name = 'VARIANT', value.name='ALLELE') %>%
    .[, .(NUM.ALLELES=sum(!is.na(ALLELE))), by=LOCUS] %>%
    data.table::merge.data.table(x=., y=uniq.loc.tab, by='LOCUS')

  # Output
  return(uniq.loc.tab)
}
