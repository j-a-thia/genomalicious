#' Convert a long-format data table of genotypes into a genind/genlight object
#'
#' For a data table with samples and loci as rows (long-format data table),
#' create a genind/genlight object, as per the \code{adegenet} package
#' (Jombart 2011 Bioinformatics). Data must be biallelic SNPs.
#'
#' @param dat Data.table: Contains samples and loci in rows, with a
#' separate column endcoding genotypes, where each allele is separated
#' by a '/', or as a count of ref alleles (0, 1 or 2). Assumes biallelic loci.
#' Code missing data as './.' for characters, and \code{NA} for Ref allele counts.
#'
#' @param genX Character: One of \code{'genind'} or \code{'genlight'}.
#'
#' @param sampCol Character: Column with sample information. Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: Column with locus information. Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: Column with genotype information. Default = \code{'GT'}.
#'
#' @param popCol Character: Column with population information. Optional. Default is NULL.
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
#' # Genind without and with populations
#' adegenet_DT2genX(data_Genos, genX='genind')
#' adegenet_DT2genX(data_Genos, genX='genind', popCol='POP')
#'
#' # Genlight without and with populations
#' adegenet_DT2genX(data_Genos, genX='genlight')
#' adegenet_DT2genX(data_Genos, genX='genlight', popCol='POP')
#'
#' @export
adegenet_DT2genX <- function(dat, genX, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', popCol=NULL){

  for(lib in c('data.table', 'adegenet')){require(lib, character.only=TRUE)}

  # Create a genind object
  if(genX=='genind'){
    # What is the genotype scoring method? Convert to characters if needed.
    if(class(dat[[genoCol]])=='integer' | class(dat[[genoCol]])=='numeric'){
      dat[[genoCol]] <- genoscore_converter(dat[[genoCol]])
    }

    # Convert long format data table into a wide matrix
    mat <- DT2Mat_genos(dat, sampCol=sampCol, locusCol=locusCol, genoCol=genoCol)

    # Create the genind object
    genInd <- df2genind(X=mat, sep='/', ind.names=rownames(mat), loc.names=colnames(mat), ploidy=2, NA.char='./.')

    # Add in populations?
    if(is.null(popCol)==FALSE){
      # The population and samples
      popDat <- unique(dat[, c(sampCol, popCol), with=FALSE])
      # Assign
      genInd@pop <- as.factor(popDat[[popCol]][match(rownames(mat), popDat$SAMPLE)])
    }

    # Output
    return(genInd)
  }

  # Create a genlight object
  if(genX=='genlight'){
    # What is the genotype scoring method? Convert to counts if needed.
    if(class(dat[[genoCol]])=='character'){
      dat[[genoCol]] <- genoscore_converter(dat[[genoCol]])
    }

    genoMat <- DT2Mat_genos(dat, sampCol=sampCol, locusCol=locusCol, genoCol=genoCol)
    genLight <- new('genlight', as.list(as.data.frame(t(genoMat))))

    # Add in locus names
    genLight@loc.names <- colnames(genoMat)

    # Add in populations?
    if(is.null(popCol)==FALSE){
      # The population and samples
      popDat <- unique(dat[, c(sampCol, popCol), with=FALSE])
      genLight@pop <- as.factor(popDat[[popCol]][match(genLight@ind.names, popDat[[sampCol]])])
    }

    # Output
    return(genLight)
  }
}
