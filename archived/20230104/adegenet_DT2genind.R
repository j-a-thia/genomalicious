#' Convert a long genotype data table into a \code{genind} object
#'
#' For a data table with samples and loci as rows (long-format data table),
#' create a \code{genind} object, as per the \code{adegenet} package
#' (Jombart 2008 Bioinformatics). Does not have to be biallelic.
#'
#' @param dat Data.table: Contains samples and loci in rows, with a
#' separate column endcoding genotypes, where each allele is separated
#' by a '/'. Assumes missing genotyes are coded as './.'.
#'
#' @param sampCol Character: Column with sample information. Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: Column with locus information. Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: Column with genotype information. Default = \code{'GT'}.
#'
#' @param popCol Character: Column with population information. Optional.
#'
#' @param n Integer: The ploidy level. Default = 2.
#'
#' @return A \code{genind} object, with the slot \code{pop} slot filled if argument
#' \code{popCol} is specified.
#'
#' @references
#' Jombart (2008) adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#'
#' data_Genos
#'
#' adegenet_DT2genind(data_Genos)
#' adegenet_DT2genind(data_Genos, popCol='POP')
#'
#' @export
adegenet_DT2genind <- function(dat, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', popCol=NULL, n=2L){

  for(lib in c('data.table', 'adegenet')){require(lib, character.only=TRUE)}

  # Convert long format data table into a wide matrix
  mat <- DT2Mat_genos(dat, sampCol=sampCol, locusCol=locusCol, genoCol=genoCol)

  # Create the genind object
  genInd <- df2genind(X=mat, sep='/', ind.names=rownames(mat), loc.names=colnames(mat), ploidy=n, NA.char='./.')

  if(is.null(popCol)==FALSE){
    # The population and samples
    popDat <- unique(dat[, c(sampCol, popCol), with=FALSE])
    # Assign
    genInd@pop <- as.factor(popDat[[popCol]][match(rownames(mat), popDat$SAMPLE)])
  }

  return(genInd)
}
