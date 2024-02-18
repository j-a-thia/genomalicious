#' Annotate pairwise matrix of simulated family members
#'
#' Takes individuals generated from \code{family_sim_data} and annotates
#' familial relationships across all pairs of individuals. Takes either a pairwise
#' matrix of individuals or a character vector.
#'
#' @param simFam Matrix or Character: Either a pairwise matrix or a character vector of
#' individual sample IDs. If a pairwise matrix, the matrix must be symmetrical with sample
#' IDs in row and column names: the diagonal compares individuals to themselves,
#' and off-diagonal is comparisons between individuals; must set \code{type=="matrix"}.
#' If a character vector of sample IDs, must set \code{type=="character"}.
#'
#' @param type Character: Either "matrix" or "id_vec".
#'
#' @returns Returns a long-format data.table with the columns:
#' \enumerate{
#'    \item \code{$SIM}, the simulation number for pairs of simulated individuals,
#'    or 'NA' for the pairs of observed individuals.
#'    \item \code{$SAMPLE1}, the sample ID for the first individual.
#'    \item \code{$SAMPLE2}, the sample ID for the second individual.
#'    \item \code{$FAMILY}, the familial relationship for simulated individuals.
#'    \item \code{$Y}, the value in the cell contents from the input matrix, if
#'    \code{type=='matrix'}.
#' }
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#'
#' # Subset Pop1 genotypes
#' genosPop1 <- data_Genos[POP=='Pop1', c('SAMPLE', 'LOCUS', 'GT')]
#'
#' # Get the allele frequencies for Pop1
#' freqsPop1 <- genosPop1[, .(FREQ=sum(GT)/(length(GT)*2)), by=LOCUS]
#'
#' # Simulate 100 families
#' simFamily <- family_sim_data(
#'    freqData=freqsPop1,
#'    locusCol='LOCUS',
#'    freqCol='FREQ',
#'    numSims=100
#' )
#'
#' ### THE SIMULATED GENOMIC RELATIONSHIPS MATRIX
#' library(AGHmatrix)
#'
#' # Convert simulated families into a genotype matrix
#' simGenosMat <- DT2Mat_genos(simFamily)
#'
#' # Calculate the GRM
#' simGRM <- Gmatrix(simGenosMat, method='Yang', ploidy=2)
#'
#' ### THE FAMILAL ANNOTATIONS
#' # From a pairwise matrix (the GRM)
#' family_sim_annotate(simFam = simGRM, type = 'matrix')
#'
#' # From a vector of unique individual IDs
#' family_sim_annotate(simFam = rownames(simGRM), type = 'id_vec')
#'
#' @export

family_sim_annotate <- function(simFam, type){

  if(type=='matrix'){
    if(!'matrix' %in% class(simFam)){
      stop('Argument `simFam` must be of class "matrix" if argument `type` == "matrix". See ?family_sim_annotate.')
    }
  }

  if(type=='id_vec'){
    if(!'character' %in% class(simFam)){
      stop('Argument `simFam` must be of class "character" if argument `type` == "id_vec". See ?family_sim_annotate.')
    }
  }

  require(data.table); require(tidyverse)

  # Sort the individual sample IDs based on input type
  if(type=='matrix'){
    annotDT <- simFam %>%
      as.matrix() %>%
      pairwiseMat2DT(., X1='SAMPLE1', X2='SAMPLE2', Y='Y') %>%
      .[SAMPLE1!=SAMPLE2]
  } else if(type=='id_vec'){
    annotDT <- combn(sort(unique(simFam)), 2) %>%
      t() %>%
      as.data.table %>%
      setnames(., new=c('SAMPLE1', 'SAMPLE2'))
  }

  # Simulation ID
  annotDT[, SIM1:=sub('_.*', '', SAMPLE1)]
  annotDT[, SIM2:=sub('_.*', '', SAMPLE2)]

  # All individuals from different simulations are unrelated
  annotDT[SIM1!=SIM2, FAMILY:='Unrelated']

  # The simulated unrelated pair within simulations
  annotDT[
    SIM1 == SIM2 &
      grepl('_unrel_',SAMPLE1)==TRUE &
      grepl('_unrel_',SAMPLE2)==TRUE,
    FAMILY:='Unrelated'
  ]

  # The simulated sibling pair within simulations
  annotDT[
    SIM1 == SIM2 &
      grepl('_sib_',SAMPLE1)==TRUE &
      grepl('_sib_',SAMPLE2)==TRUE,
    FAMILY:='Siblings'
  ]

  # The simulated half-sibling pair within simulations
  annotDT[
    SIM1 == SIM2 &
      grepl('_halfsib_',SAMPLE1)==TRUE &
      grepl('_halfsib_',SAMPLE2)==TRUE,
    FAMILY:='Half-siblings'
  ]

  # The simulated cousins pair within simulations
  annotDT[
    SIM1 == SIM2 &
      grepl('_cous_',SAMPLE1)==TRUE &
      grepl('_cous_',SAMPLE2)==TRUE,
    FAMILY:='Cousins'
  ]

  # The simulated half-cousins pair within simulations
  annotDT[
    SIM1 == SIM2 &
      grepl('_halfcous_',SAMPLE1)==TRUE &
      grepl('_halfcous_',SAMPLE2)==TRUE,
    FAMILY:='Half-cousins'
  ]

  # All other pairs within simulations are unrelated
  annotDT[SIM1 == SIM2 & is.na(FAMILY), FAMILY:='Unrelated']

  # Add in a simulation column
  annotDT <- annotDT[, SIM:=paste0(SIM1, '|', SIM2)] %>%
    .[, !c('SIM1','SIM2')]

  # Output
  if(type=='matrix'){
    return(annotDT[, c('SIM','SAMPLE1','SAMPLE2','FAMILY','Y')])
  } else if(type=='id_vec'){
    return(annotDT[, c('SIM','SAMPLE1','SAMPLE2','FAMILY')])
  }
}
