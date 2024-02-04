#' Simulate families of individuals from population allele frequencies
#'
#' This function produces families of siblings, half-siblings, cousins and
#' half-cousins from observed population allele frequencies. These simulations can
#' then be used to test the power to discern between individuals with different
#' levels of relatedness. Assumes diploid genotypes.
#'
#' The output can be used to generate a simulated genetic relationship matrix (GRM)
#' that can be compared against an observed GRM using the function,
#' \code{family_sim_compare}. This can provide a graphical comparison
#' between the observed and expected (simulated) distribution of relatedness
#' values that you might expect for different familial relationships, given the
#' number of loci and their allele frequencies.
#'
#' @param freqData Data.table: Population allele frequencies.
#'
#' @param locusCol Character: Column name in \code{freqData} with the locus IDs.
#' Default is 'LOCUS'.
#'
#' @param freqCol Character: Column name in \code{freqData} with the allele frequencies.
#' Default is 'FREQ'.
#'
#' @param numSims Integer: The number of simulated individuals for each family
#' relationship. Default is 100.
#'
#' @details In this function, a 'simulation' comprises a draw of 10 individuals:
#' 2 are unrelated to each other and all other individuals in the simulation,
#' 2 are siblings, 2 are half-siblings, 2 are cousins, and 2 are half-cousins.
#'
#' Each simulation takes the naming convention 'S[sim number]', for example, S1
#' S2 ... S10 if 10 simulations are requested. The naming convention for each
#' simulated individual is 'S[sim number]_[relationship]', with the 'relationship'
#' tag being one of the following:
#' \itemize{
#'  \item 'unrel_1' and 'unrel_2' are two completely unrelated individuals.
#'  \item 'sib_1' and 'sib_2' are full siblings.
#'  \item 'halfsib_1' are 'halfsib_2' are half siblings.
#'  \item 'cous_1' and 'cous_2' are cousins.
#'  \item 'halfcous_1' and 'halfcous_2' are half cousins.
#' }
#' Note, that within the sth simulation, all pairs are unrelated to each other as
#' well. For example, in simulation S1, sib_1 and sib_2 are completely unrelated
#' to halfsib_1 and halfsib_2. This means the number of unrelated pairs within
#' an sth simulation is 45 pairs, 4 of these are related (siblings, half-siblings,
#' cousins, and half-cousins), and 41 are unrelated. Individuals from different
#' simulation are completely unrelated. For example, all individuals from simulation
#' S1 are completely unrelated to all individuals from simulation S2.
#'
#' @returns Returns a data.table with the following columns:
#' \enumerate{
#'    \item \code{$SIM}, the simulation number.
#'    \item \code{$LOCUS}, the locus ID.
#'    \item \code{$SAMPLE}, the sample ID.
#'    \item \code{$GT}, the diploid genotype as counts of the Alt alleles.
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
#' # Create some siblings in Pop1 from two sets of parents
#' parentList <- list(c('Ind1_1','Ind1_2'), c('Ind1_3','Ind1_4'))
#' genosSibs <- lapply(1:2, function(i){
#'   parents <- parentList[[i]]
#'
#'   child <- paste(sub('Ind1_', '', parents), collapse='.')
#'
#'   gamete1 <- genosPop1[SAMPLE == parents[1]] %>%
#'     .[, .(GAMETE=rbinom(n=2,size=1,prob=GT/2)), by=c('SAMPLE','LOCUS')] %>%
#'     .[, SIB:=1:2, by=LOCUS]
#'
#'   gamete2 <- genosPop1[SAMPLE == parents[2]] %>%
#'     .[, .(GAMETE=rbinom(n=2,size=1,prob=GT/2)), by=c('SAMPLE','LOCUS')] %>%
#'     .[, SIB:=1:2, by=LOCUS]
#'
#'   rbind(gamete1, gamete2) %>%
#'     .[, .(GT=sum(GAMETE)), by=c('LOCUS','SIB')] %>%
#'     .[, SAMPLE:=paste0('Child_',child,'_',SIB)]
#' }) %>%
#'   do.call('rbind', .)
#'
#' ### THE OBSERVED GENOMIC RELATIONSHIPS MATRIX
#' library(AGHmatrix)
#'
#' # Combine the population samples and the created siblings
#' # into a single genotype matrix
#' obsGenosMat <- rbind(genosPop1, genosSibs[, c('SAMPLE','LOCUS','GT')]) %>%
#'   DT2Mat_genos()
#'
#' # Calculate the GRM
#' obsGRM <- Gmatrix(obsGenosMat, method='Yang', ploidy=2)
#'
#' ### THE SIMULATED GENOMIC RELATIONSHIPS MATRIX
#' # Convert simulated families into a genotype matrix
#' simGenosMat <- DT2Mat_genos(simFamily)
#'
#' # Calculate the GRM
#' simGRM <- Gmatrix(simGenosMat, method='Yang', ploidy=2)
#'
#' ### COMPARE THE OBSERVED AND SIMULATED
#' relComp <- family_sim_compare(
#'    simGRM=simGRM,
#'    obsGRM=obsGRM,
#'    look='classic'
#' )
#'
#' # The data
#' relComp$data
#'
#' # Simulated dataset
#' relComp$data[!is.na(SIM)]
#'
#' # The observed dataset
#' relComp$data[is.na(SIM)]
#'
#' # Plot of relatedness values. Dashed lines denote relatedness
#' # values of 0, 0.0625, 0.125, 0.25, and 0.5, which are the theoretical
#' # expectations for unrelated individuals, half-cousins, cousins, half-siblings,
#' # and siblings, respectively.
#' # You will note a large variance are the expected values, which
#' # is not surprising for this very small SNP dataset (200 loci).
#' relComp$plot
#'
#' # Take a look at the "known" relationships in the observed dataset using the
#' # offspring we created.
#' # Note, siblings and parent-offspring pairs have a theoretical
#' # relatedness of 0.5. But you will probably find the "observed"
#' # relatedness might be lower.
#' relComp$data[SAMPLE1=='Child_1.2_1' & SAMPLE2%in%c('Child_1.2_2','Ind1_1','Ind1_2')]
#' relComp$data[SAMPLE1=='Child_3.4_1' & SAMPLE2%in%c('Child_3.4_2','Ind1_3','Ind1_4')]
#'
#' # Now take a look at the simulated distribution: the possible values for
#' # siblings given this dataset are quite wide.
#' relComp$data[FAMILY=='Half-siblings']$RELATE %>% summary()
#' relComp$data[FAMILY=='Siblings']$RELATE %>% summary()
#'
#' @export

family_sim_data <- function(freqData, locusCol='LOCUS', freqCol='FREQ', numSims=100L){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('tidyr', 'data.table')){require(lib, character.only=TRUE)}

  # Check if data.table/can be converted to data.table
  freqData <- as.data.table(freqData)
  if(!'data.table' %in% class(freqData)){
    stop('Argument `freqData` must be a data.table class. See ?sim_family.')
  }

  # Check that columns are specified correctly
  column.check <- sum(c(locusCol,freqCol) %in% colnames(freqData))
  if(column.check!=2){
    stop(
      'Arguments `locusCol` and `freqCol` must be column names in the argument
    freqData. See ?sim_family.')
  }

  # Reassign column names.
  freqData <- freqData %>%
    copy %>%
    setnames(., old=c(locusCol,freqCol), new=c('LOCUS','FREQ'))

  # Check that numSims is >0 and is an integer
  numSims <- as.integer(numSims)
  if(!numSims>0){
    stop('Argument numSims must be an integer >0. See ?sim_family.')
  }

  # --------------------------------------------+
  # Internal functions
  # --------------------------------------------+
  FUN_draw_genos <- function(D, ploidy){
    # Function to create genotypes assuming random binomial draw at
    # each locus. Assumes loci are unlinked.

    # D = data.table of population allele freqs, with columns $LOCUS and $FREQ.
    # ploidy = integer of number of allles per genotype

    D[, .(FREQ=rbinom(n=1, size=ploidy, prob=FREQ)/ploidy), by=LOCUS] %>%
      data.table(PLOIDY=ploidy, .)
  }

  FUN_make_diploid_offspring <- function(D1, D2){
    # Function to create diploid offspring from the allele frequencies of two
    # diploid individuals.

    # D1 and D2 = parental allele frequencies as a data.table for the 1st and 2nd
    # parent, respectively. Contain $LOCUS and $FREQ. In this case, allele freqs
    # represent the relative dosage of a diploid genotype, so 0 = 0/0, 0.5 = 0/1,
    # and 1 = 1/1.

    left_join(
      FUN_draw_genos(D1, ploidy=1) %>%
        setnames(., old='FREQ', new='GAMETE1') %>%
        .[, c('LOCUS','GAMETE1')],
      FUN_draw_genos(D2, ploidy=1) %>%
        setnames(., old='FREQ', new='GAMETE2')%>%
        .[, c('LOCUS','GAMETE2')],
      by='LOCUS'
    ) %>%
      as.data.table %>%
      .[, FREQ:=(GAMETE1 + GAMETE2)/2] %>%
      .[, c('LOCUS','FREQ')] %>%
      data.table(PLOIDY=2, .)
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  simFamily <- lapply(1:numSims, function(sim){
    # Note, 'g[x]' indicates generation in each group.

    # Create sibling pair
    sib.g1.1 <- FUN_draw_genos(freqData, ploidy=2)
    sib.g1.2 <- FUN_draw_genos(freqData, ploidy=2)

    sib.g2.1 <- FUN_make_diploid_offspring(sib.g1.1, sib.g1.2)
    sib.g2.2 <- FUN_make_diploid_offspring(sib.g1.1, sib.g1.2)

    # Create half sibling pair
    hsib.g1.1 <- FUN_draw_genos(freqData, ploidy=2)
    hsib.g1.2 <- FUN_draw_genos(freqData, ploidy=2)
    hsib.g1.3 <- FUN_draw_genos(freqData, ploidy=2)

    hsib.g2.1 <- FUN_make_diploid_offspring(hsib.g1.1, hsib.g1.2)
    hsib.g2.2 <- FUN_make_diploid_offspring(hsib.g1.1, hsib.g1.3)

    # Create cousin pair
    cuz.g1.1 <- FUN_draw_genos(freqData, ploidy=2)
    cuz.g1.2 <- FUN_draw_genos(freqData, ploidy=2)

    cuz.g2.1 <- FUN_draw_genos(freqData, ploidy=2)
    cuz.g2.2 <- FUN_draw_genos(freqData, ploidy=2)
    cuz.g2.3 <- FUN_make_diploid_offspring(cuz.g1.1, cuz.g1.2)
    cuz.g2.4 <- FUN_make_diploid_offspring(cuz.g1.1, cuz.g1.2)

    cuz.g3.g1 <- FUN_make_diploid_offspring(cuz.g2.3, cuz.g2.1)
    cuz.g3.g2 <- FUN_make_diploid_offspring(cuz.g2.4, cuz.g2.2)

    # Create half-cousin pair
    hcuz.g1.1 <- FUN_draw_genos(freqData, ploidy=2)
    hcuz.g1.2 <- FUN_draw_genos(freqData, ploidy=2)
    hcuz.g1.3 <- FUN_draw_genos(freqData, ploidy=2)

    hcuz.g2.1 <- FUN_draw_genos(freqData, ploidy=2)
    hcuz.g2.2 <- FUN_draw_genos(freqData, ploidy=2)
    hcuz.g2.3 <- FUN_make_diploid_offspring(hcuz.g1.1, hcuz.g1.2)
    hcuz.g2.4 <- FUN_make_diploid_offspring(hcuz.g1.1, hcuz.g1.3)

    hcuz.g3.1 <- FUN_make_diploid_offspring(hcuz.g2.3, hcuz.g2.1)
    hcuz.g3.2 <- FUN_make_diploid_offspring(hcuz.g2.4, hcuz.g2.2)

    # Random completely unrelated individuals
    unrel.g1.1 <- FUN_draw_genos(freqData, ploidy=2)
    unrel.g1.2 <- FUN_draw_genos(freqData, ploidy=2)

    rbind(
      unrel.g1.1 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_unrel_1')),
      unrel.g1.2 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_unrel_2')),
      sib.g2.1 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_sib_1')),
      sib.g2.2 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_sib_2')),
      hsib.g2.1 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_halfsib_1')),
      hsib.g2.2 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_halfsib_2')),
      cuz.g3.g1 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_cous_1')),
      cuz.g3.g2 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_cous_2')),
      hcuz.g3.1 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_halfcous_1')),
      hcuz.g3.2 %>% data.table(SIM=sim, SAMPLE=paste0('S',sim,'_halfcous_2'))
    )
  }) %>%
    do.call('rbind', .) %>%
    .[, GT:=FREQ*PLOIDY]

  # Output
  return(simFamily[, c('SIM','SAMPLE','LOCUS','GT')])
}
