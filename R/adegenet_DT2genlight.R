#' Convert a long genotype data table into a \code{genlight} object
#' 
#' For a data table with samples and loci as rows (long-format data table),
#' create a \code{genlight} object, as per the \code{adegenet} package
#' (Jombart 2011 Bioinformatics). Must be biallelic.
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
#' @param popCol Character: Column with population information. Optional.
#' 
#' @return A \code{genind} object, with the slot \code{pop} slot filled if argumnet
#' \code{popCol} is specified.
#' 
#' @references 
#' Jombart (2008) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics.
#' 
#' @examples 
#' data(genomalicious_4pops)
#' 
#' genomalicious_4pops
#' 
#' adegenet_DT2genlight(genomalicious_4pops)
#' adegenet_DT2genlight(genomalicious_4pops, popCol='POP')
#' 
#' @export
adegenet_DT2genlight <- function(dat, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', popCol=NULL){
  
  for(lib in c('data.table', 'adegenet')){require(lib, character.only=TRUE)}
  
  # What is the genotype scoring method?
  if(class(dat[[genoCol]])=='character'){
    gt.score <- 'sep'
  } else if(class(dat[[genoCol]])=='integer'){
    gt.score <- 'counts'
  }
  
  # Split data by samples into a list
  datspl <- split(dat, dat[[sampCol]])
  
  # Iterate through samples, reorder loci, and get genotypes
  if(gt.score=='sep'){
    alleles <- lapply(datspl, function(samp){
      return(genoscore_converter(samp[order(LOCUS)]$GT))
    })
  } else{
    alleles <- lapply(datspl, function(samp){
      return(samp[order(LOCUS)]$GT)})
  }

  # Convert to genlight object
  genLight <- new('genlight', alleles)
  
  if(is.null(popCol)==FALSE){
    # The population and samples
    popDat <- unique(dat[, c(sampCol, popCol), with=FALSE])
    genLight@pop <- as.factor(popDat[[popCol]][match(names(alleles), popDat[[sampCol]])])
  }
  
  return(genLight)
}
  