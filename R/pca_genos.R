#' Perform a PCA (principal components analysis) on individual genotypes
#'
#' Takes a long-format data table of genotypes and conducts a PCA using R's
#' \code{prcomp()} function. Different options for scaling the genotypes
#' pre-PCA are available.
#'
#' @param dat Data table: A long-fomate data table. Genotypes can be coded as
#' '/' separated characters (e.g. '0/0', '0/1', '1/1'), or integers of Alt
#' allele counts (e.g. 0, 1, 2). Must contain the following columns,
#' \enumerate{
#'   \item The sampled individuals (see param \code{sampCol}).
#'   \item The locus ID (see param \code{locusCol}).
#'   \item The genotype column (see param \code{genoCol}).
#' }
#' Optionally, a population ID column can also be included (see param \code{popCol}).
#'
#' @param scaling Character: How should the data (loci) be scaled?
#' Default is \code{'covar'} to scale to mean = 0, but variance is not
#' adjusted, i.e. PCA on a covariance matrix. Set to \code{'corr'}
#' to scale to mean = 0 and variance = 1, i.e. PCA on a
#' correlation matrix. Set to \code{'patterson'} to use the
#' Patteron et al. (2006) normalisation. Set to \code{'none'} to
#' if you do not want to do any scaling before PCA.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default is \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default is \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default is \code{'GT'}.
#'
#' @param popCol Character: An optional argument. The column name with the
#' population information. Default is \code{NULL}. If specified, population
#' membership is stored in the returned object.
#'
#' @return Returns a \code{prcomp} object. If argument \code{popCols} was specified,
#' and additional index of \code{$pops} is also also present.
#'
#' @references
#' Patterson et al. (2006) Population structure and eigenanalysis. PLOS Genetics.
#'
#' @examples
#' library(genomalicious)
#'
#' # Data
#' data(data_Genos)
#' data_Genos
#'
#' # Conduct the PCA with Patterson et al.'s (2006) normalisation, and
#' # population specified
#' pca <- pca_genos(dat=data_Genos, scaling='patterson', popCol='POP')
#'
#' # Plot the PCA
#' pca_plot(pca)
#'
#' @export

pca_genos <- function(dat, scaling='covar', sampCol='SAMPLE'
                       , locusCol='LOCUS', genoCol='GT', popCol=NULL){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'dplyr')){ require(lib, character.only = TRUE)}

  # Check that scaling is specified
  if(!scaling %in% c('covar', 'corr', 'patterson', 'none')){
    stop('Argument `scaling`` is invalid. See: ?pca_genos')
  }

  # Get the class of the genotypes
  gtClass <- class(dat[[genoCol]])

  # Check that genotypes are characters or counts
  if(!gtClass %in% c('character', 'numeric', 'integer')){
    stop("Check that genotypes are coded as '/' separated characters or as
         counts of the Alt allele. See: ?pca_genos")
  }

  # Convert characters of separated alleles to counts
  if(gtClass=='character'){
    dat[[genoCol]] <- genoscore_converter(dat[[genoCol]])
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Convert to a genotype matrix
  genoMat <- DT2Mat_genos(dat
                          , sampCol=sampCol
                          , locusCol=locusCol
                          , genoCol=genoCol)

  # Get individuals in rows
  sampRows <- rownames(genoMat)

  # The PCA
  if(scaling=='covar'){
    pca <- prcomp(genoMat, center=TRUE, scale=FALSE)
  } else if(scaling=='corr'){
    pca <- prcomp(genoMat, center=TRUE, scale=TRUE)
  } else if(scaling=='patterson'){
    pca <- prcomp(normalise_patterson(genoMat), center=FALSE, scale=FALSE)
    pca$scale <- 'Patterson et al. (2006)'
  }
  rownames(genoMat) <- sampRows

  # Was the population ID column specified?
  if(is.null(popCol)==FALSE){
    # Get population memberships for each sample
    popMems <- dat[, c(popCol, sampCol), with=FALSE]
    popMems <- unique(popMems)
    # Match individuals in rows from genoMat with popMems
    pca$pops <- popMems[[popCol]][match(sampRows, popMems[[sampCol]])]
  }

  return(pca)

  # ........... END
}
