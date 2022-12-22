#' Convert a genotype data table to an allele data frame
#'
#' A long-format data table of genotypes is converted into a wide-format
#' data table of alleles, where each locus is represented by two columns:
#' one for the first allele, the other for the second allele.
#'
#' @param dat Data.table: Contains samples and loci in rows, with a
#' separate column endcoding genotypes, where each allele is separated
#' by a '/'.
#'
#' @param sampCol Character: Column with sample information. Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: Column with locus information. Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: Column with genotype information. Default = \code{'GT'}.
#'
#' @return A data frame in wide format, with a sample column \code{$SAMPLE} and
#' 2 * number of loci additional columns containing alleles for each locus.
#' For example, \code{$locus1_1} and \code{$locus1_2} will contain the first and
#' second allele for "locus1".
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#'
#' data_Genos
#'
#' alle_dt <- genos2indAlleles(data_Genos)
#'
#' # The first three loci, with count of allele 1 and 2.
#' alle_dt[1:4, 1:7]
#'
#' @export
genos2indAlleles <- function(dat, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT'){

  require(data.table); require(tidyverse)

  if(sum(c(sampCol, locusCol, genoCol) %in% colnames(dat))!=3){
    stop('One of arguments sampCol, locusCol, genoCol is not in dat. See ?genoDT2alleleDF')
  }

  if(class(dat$GT) %in% c('integer', 'numeric')){
    dat$GT <- genoscore_converter(dat$GT)
  }

  # Unique loci
  uniq_loci <- unique(dat$LOCUS)

  # Convert data table to wide matrix of genotypes
  mat <- DT2Mat_genos(dat, 'SAMPLE', 'LOCUS', 'GT')

  # Create the allele data frame.
  # Iterate over loci.
  alleDF <- lapply(uniq_loci, function(loc){
    # Get genotypes for a locus
    gen <- mat[,loc]
    # Split the genotypes by '/', unlist, convert to integer,
    # then convert to 2 column matrix
    alle <- strsplit(gen, split='/') %>%
      unlist() %>%
      as.integer() %>%
      matrix(., ncol=2, byrow=TRUE)
    # Colnames
    colnames(alle) <- paste(loc, 1:2, sep='_')
    # Return
    return(alle)
  }) %>%
    do.call('cbind', .) %>%
    as.data.table() %>%
    data.table(SAMPLE=rownames(mat), .)

  return(alleDF)
}
