#' The additive genetic relatedness matrix for simulated individuals
#'
#' Takes individuals generated from \code{family_sim_data} and returns an
#' additive genetic relatedness matrix
#'
#' @param simFam Character: A vector of simulated indivdiual sample IDs
#'
#' @returns Returns a symmetrical matrix of additive genetic relationships
#' between all pairs of individuals.
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
#' ### THE SIMULATED ADDITIVE GENETIC RELATEDNESS MATRIX
#' simARM <- family_sim_Amatrix(unique(simFamily$SAMPLE))
#'
#' simARM[1:8, 1:8]
#'
#' @export

family_sim_Amatrix <- function(simFam, type){

  require(data.table); require(tidyverse)

  Amat <- family_sim_annotate(simFam, type='id_vec') %>%
    .[FAMILY=='Unrelated', A:=0] %>%
    .[FAMILY=='Half-cousins', A:=0.0625] %>%
    .[FAMILY=='Cousins', A:=0.125] %>%
    .[FAMILY=='Half-siblings', A:=0.25] %>%
    .[FAMILY=='Siblings', A:=0.5] %>%
    pairwiseMat2DT(., flip=TRUE, X1='SAMPLE1', X2='SAMPLE2', Y='A', diagAdd=TRUE, diagVal=1)

  return(Amat)
}
