#' Simulate missing data with a specified strcuture
#'
#' Takes a long format data table of genotypes that
#' has no missing data (clean) and simulates missing data following
#' specified per locus proportions.
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
#' @param miss_loci Numeric: Vector of proportions. Each value represents the
#' proportion of missing data at a locus.
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
#' @details The number of values in \code{miss_loci} does not have to be equal to
#' \code{dat_clean}. In the simulations, the values in \code{miss_loci}
#' are randomly sampled (with replacement) to create a bootstrap simulated
#' missing data structure of the same dimensions (loci and samples) as the
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
#' library(tidyverse)
#' data(data_4pops)
#'
#' # Create a vector of proportions
#' prop_miss <- rbeta(1000, 0.5, 10)
#' hist(prop_miss, main='Input missing proportions')
#'
#' # Simulate missing data structure
#' patchy4pops <- miss_sim_structure(
#'     dat_clean=data_4pops
#'     , miss_loci=prop_miss
#'     , sampCol='SAMPLE'
#'     , locusCol='LOCUS'
#'     , genoCol='GT'
#' )
#'
#' patchy4pops[, sum(GT=='./.')/length(GT), by=LOCUS]$V1 %>%
#' hist(main='Simulated missing distribution')
#'
#' @export
#'
miss_sim_structure <- function(dat_clean, miss_loci
                               , sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'tidyverse')){ require(lib, character.only=TRUE)}

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

  # Unique loci and samples
  uniqloci <- unique(dat_clean[[locusCol]])
  uniqsamps <- unique(dat_clean[[sampCol]])

  # Sampling parameters
  nloci <- length(uniqloci)
  nsamps <- length(uniqsamps)

  # Make a bootstrap missing dataset to have the
  # same number of loci and samples as in observed dataset, `dat_clean`
  mat_boot <- round(sample(miss_loci, nloci, TRUE) * nsamps) %>%
    as.integer() %>%
    lapply(., function(ii){
      c(rep('X', nsamps-ii), rep(NA, ii)) %>%
        sample(., replace=FALSE)
    }) %>%
    do.call('cbind', .)

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
