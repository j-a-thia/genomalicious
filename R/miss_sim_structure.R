#' Simulate missing data with a specified strcuture
#'
#' Takes a long format data table of genotypes that
#' has no missing data (clean) and simulates a missing data
#' analogous to a pathcy matrix used as a guide.
#'
#' @param dat_clean Data table: Must be in long format and contains
#' no missing data (i.e. it is "clean"). Genotypes can be coded as
#' '/' separated characters or a integers. Must contain the following
#' columns:
#' \enumerate{
#'    \item The sampled individual ID (see param \code{sampCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The genotype (see param \code{genoCol}).
#' }
#'
#' @param mat_patchy Matrix: Samples in rows, loci in columns.
#' This is meant to simulate a genotype matrix with missing data
#' (i.e. it is "patchy"), but the the actual contents of the matrix
#' can be anything, as long as missing values are coded by \code{NA}.
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
#' @details The number of samples and loci in the guide matrix, \code{mat_patchy}
#' does not have to be equal \code{dat_clean}. In the simulations, the rows and
#' columns of \code{mat_patchy} are randomly sampled (with replacement) to create
#' a bootstrap simulated  missing data structure of the same dimensions as the
#' observed clean data. This simulated missing data structure is then used to
#' modify the observed clean data. \cr
#' \cr
#' NOTE: This method is blind to population, so will simulate
#' missing data across individuals. Therefore, if missing data exhibits
#' population structuring, it would be necessary to subset the clean dataset
#' by population and run specific simulations per population, then combine
#' all the simulated subsets.
#'
#' @values A data table with the same dimensions as \code{dat_clean}, but
#' genotypes will be knocked out to reproduce a data structure anologous to
#' that specified in \code{mat_patchy}. If genotypes were initially coded
#' as characters, the missing genotypes will be coded as './.'; if they were
#' initially coded as integers, they will be coded as \code{NA}.
#'
#' @examples
#' library(data.table)
#' data(data_PatchyGTs)
#' data(data_4pops)
#'
#' # Take a look at the patchy dataset
#' data_PatchyGTs
#'
#' # Simulate missing data structure
#' patchy4pops <- miss_sim_structure(
#'     dat_clean=data_4pops
#'     , mat_patchy=data_PatchyGTs
#'     , sampCol='SAMPLE'
#'     , locusCol='LOCUS'
#'     , genoCol='GT'
#' )
#'
#' # Take a look at the output.
#' # Histograms are ordered with samples in the first
#' # row and loci in the second row. The first column
#' # is the observed clean data, the second column
#' # is the simulated missing data, and the thir column
#' # is the missing data in the patchy guide matrix.
#' par(mfrow=c(2,3))
#' hist(data_4pops[, sum(GT=='./.')/length(GT), by=SAMPLE]$V1
#'     , 100, xlim=c(0,1), main='Obs clean: Samples', xlab='% missing')
#'
#' hist(patchy4pops[, sum(GT=='./.')/length(GT), by=SAMPLE]$V1
#'     , 100, xlim=c(0,1), main='Sim patchy: Samples', xlab='% missing')
#'
#' hist(
#'     unlist(
#'         apply(data_PatchyGTs
#'            , 1
#'            , function(i){sum(is.na(i))/length(i)}))
#'     , 100, xlim=c(0,1), main='Guide patchy: Samples', xlab='% missing')
#'
#' hist(data_4pops[, sum(GT=='./.')/length(GT), by=LOCUS]$V1
#'     , 100, xlim=c(0,1), main='Obs clean: Loci', xlab='% missing')
#'
#' hist(patchy4pops[, sum(GT=='./.')/length(GT), by=LOCUS]$V1
#'     , 100, xlim=c(0,1), main='Sim patchy: Loci', xlab='% missing')
#'
#' hist(
#'     unlist(
#'         apply(data_PatchyGTs
#'            , 2
#'            , function(i){sum(is.na(i))/length(i)}))
#'     , 100, xlim=c(0,1), main='Guide patchy: Loci', xlab='% missing')
#' par(mfrow=c(1,1))
#'
#' @export
#'
miss_sim_structure <- function(dat_clean, mat_patchy
                               , sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'tidyr', 'plyr')){ require(lib, character.only=TRUE)}

  # Genotype column class
  gt_class <- class(dat_clean[[genoCol]])

  if(gt_class=='numeric'){
    dat_clean[[genoCol]] <- as.integer(dat_clean[[genoCol]])
  }

  if(gt_class!='character'){
    dat_clean[[genoCol]] <- genoscore_converter(dat_clean[[genoCol]])
  }

  # Check the column arguments are specified
  if(sum(length(c(sampCol, locusCol, genoCol) %in% colnames(dat_clean)))!=3){
    stop('Required columns not in `dat_clean`: see ?miss_sim_structure')
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Sampling parameters
  nloci <- length(unique(dat_clean[[locusCol]]))
  nsamps <- length(unique(dat_clean[[sampCol]]))

  # Make a bootstrap missing dataset to have the
  # same number of loci and samples as in observed dataset, `dat_clean`
  mat_boot <- mat_patchy[
    sample(x=1:nrow(mat_patchy), size=nsamps, replace=TRUE)
    , sample(x=1:ncol(mat_patchy), size=nloci, replace=TRUE)
    ]

  # The simulated missing dataset:
  # Start with the original, convert to a wide format matrix, then
  # add in missing values that match positions in `mat_boot`
  mat_sim <- DT2Mat_genos(
    dat_clean
    , sampCol='SAMPLE'
    , locusCol='LOCUS'
    , genoCol='GT')

  mat_sim[which(is.na(mat_boot))] <- './.'

  # Convert the simulated dataset back into a long format data table.
  # Then subset to only include sample/loci combinations with missing data.
  ko_genos <- DT2Mat_genos(
    mat_sim
    , sampCol=sampCol
    , locusCol=locusCol
    , genoCol=genoCol
    , flip=TRUE)[GT=='./.']

  # Make a copy of the clean dataset
  dat_sim <- copy(dat_clean)

  # Unique identifiers for each sample/locus combination in the
  # simulated and observed datasets
  ko_genos$SAMPLE.LOCUS <- paste(ko_genos$SAMPLE, ko_genos$LOCUS, sep='/')
  dat_sim$SAMPLE.LOCUS <- paste(dat_sim[[sampCol]], dat_sim[[locusCol]], sep='/')

  # Get the index of the missing sample/loci combinations in the
  # observed dataset, then knock 'em out.
  ko_indx <- match(ko_genos$SAMPLE.LOCUS, dat_sim$SAMPLE.LOCUS)
  dat_sim[ko_indx, genoCol] <- './.'

  # Convert genotypes back to integer if that was the original specification
  if(gt_class %in% c('numeric', 'integer')){
    dat_clean[[genoCol]] <- genoscore_converter(dat_clean[[genoCol]])
  }

  # Return
  return(dat_sim[, !'SAMPLE.LOCUS'])
}
