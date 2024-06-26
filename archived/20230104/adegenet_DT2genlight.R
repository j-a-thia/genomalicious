#' Convert a long genotype data table into a \code{genlight} object
#'
#' For a data table with samples and loci as rows (long-format data table),
#' create a \code{genlight} object, as per the \code{adegenet} package
#' (Jombart 2011 Bioinformatics). Must be biallelic.
#'
#' @param dat Data.table: Contains samples and loci in rows, with a
#' separate column endcoding genotypes, where each allele is separated
#' by a '/', or as a count of ref alleles (0, 1 or 2). Assumes biallelic loci.
#' Code missing data as './.' for characters, and \code{NA} for Ref allele counts.
#'
#' @param sampCol Character: Column with sample information. Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: Column with locus information. Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: Column with genotype information. Default = \code{'GT'}.
#'
#' @param popCol Character: Column with population information. Optional.
#'
#' @return A \code{genind} object, with the slot \code{pop} slot filled if argumnet
#' \code{popCol} is specified.
#'
#' @references
#' Jombart (2008) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#'
#' data_Genos
#'
#' adegenet_DT2genlight(data_Genos)
#' adegenet_DT2genlight(data_Genos, popCol='POP')
#'
#' @export
adegenet_DT2genlight <- function(dat, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', popCol=NULL){

  for(lib in c('data.table', 'adegenet')){require(lib, character.only=TRUE)}

  # What is the genotype scoring method? Convert to counts if needed.
  if(class(dat[[genoCol]])=='character'){
    dat[[genoCol]] <- genoscore_converter(dat[[genoCol]])
  }

  # Create a genlight object
  genoMat <- DT2Mat_genos(dat, 'SAMPLE', 'LOCUS', 'GT')
  genLight <- new('genlight', as.list(as.data.frame(t(genoMat))))

  # Add in locus names
  genLight@loc.names <- colnames(genoMat)

  if(is.null(popCol)==FALSE){
    # The population and samples
    popDat <- unique(dat[, c(sampCol, popCol), with=FALSE])
    genLight@pop <- as.factor(popDat[[popCol]][match(genLight@ind.names, popDat[[sampCol]])])
  }

  return(genLight)
}
