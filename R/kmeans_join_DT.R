#' Join K-means clustering results into a reference (long) \code{data.table}
#'
#' Takes the results of K-means clustering performed by \code{kmeans} and
#' joins them into a reference data table including the same samples.
#'
#' @param refTab Data.table: the reference data.table with a column of sample
#' names, see param \code{sampCol}.
#'
#' @param sampCol Character: a column in \code{refTab} where sample names are
#' stored. These sample names must match those in \code{names(kFit$cluster)},
#' see param \code{kFit}. Default is \code{SAMPLE}.
#'
#' @param kFit Kmeans: an object produced by the function \code{kmeans}.
#' The cluster predictions in \code{kFit$cluster} will be aligned against samples
#' in \code{refTab[[sampCol]]}. It is therefore necessary that the names in
#' this object match those expected in \code{refTab}.
#'
#' @details It is important that all the samples in \code{refTab} match the names
#' of samples in \code{kFit}. The function will throw a warning if this is not
#' the case, but will force the join regardless.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#' data_Genos
#'
#' # Perform PCA.
#' PCA <- pca_genos(data_Genos)
#'
#' # Get a table of PC axes
#' pcTable <- pca2DT(PCA)
#'
#' # K-means clustering
#' K4 <- kmeans(PCA$x[,1:4], centers=4, nstart=10)
#'
#' K4$cluster
#'
#' # Add in inferred populations.
#' pcTable2 <- kmeans_join_DT(pcTable, 'SAMPLE', K4)
#' pcTable2
#'
#' # Plot
#' (ggplot(pcTable2, aes(x=PC1, y=PC2, colour=CLUSTER))
#'    + geom_point()
#' )
#'
#' @export

kmeans_join_DT <- function(refTab, sampCol='SAMPLE', kFit){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table); require(tidyverse)

  # Convert reference to a data.table if not already
  if(!'data.table' %in% class(refTab)){
    refTab <- as.data.table(refTab)
  }

  # Check that sample column is present in reference
  if(!sampCol %in% colnames(refTab)){
    stop('Argument `sampCol` must be a column in `refTab`. See ?kmeans_join_DT.')
  }

  # Check that the samples in the reference are same in the kmeans object
  if(sum(!refTab[[sampCol]] %in% names(kFit$cluster))>0){
    warning('Not all samples in `refTab[[sampCol]]` are present in `names(kFit$cluster)`.
            See ?kmeans_join_DT.')
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  data.table(X=names(kFit$cluster), Y=as.character(kFit$cluster)) %>%
    setnames(., old=c('X','Y'), new=c(sampCol,'CLUSTER')) %>%
    left_join(refTab, .) %>%
    return()
}
