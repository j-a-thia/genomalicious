#' Conduct a PCA on individual genotypes stored in a data table
#'
#' @param dat Data table: A long data table, e.g. like that imported from
#' \code{vcf2DT}.
#'
#' @param scaling Character: How should the data (loci) be scaled? Set to
#' \code{'center'} (default) to scale to mean = 0 and variance = 1. Set to
#' \code{'patterson'} to use the Patteron et al. (2006) normalisation.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default is \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default is \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default is \code{'GT'}.


pca_DTinds <- function(dat, scaling='center', sampCol='SAMPLE'
                       , locusCol='LOCUS', genoCol='GT'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  if(!scaling %in% c('center', 'patterson')){
    stop('Argument scaling is invalid')
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Convert to a genotype matrix
  genoMat <- DT2Mat_genos(dat, sampCol=sampCol, locusCol=locusCol, genoCol=genoCol)

  # Scale
  if(scaling=='center'){
    genoMat <- apply(genoMat, 2, scale, center=TRUE, scale=TRUE)
  }

  if(scaling=='patterson'){
    genoMat <- normalise_patterson(genoMat)
  }

  # The PCA
  pca <- prcomp(genoMat, scale=FALSE, center=FALSE)
  return(PCA)
}
