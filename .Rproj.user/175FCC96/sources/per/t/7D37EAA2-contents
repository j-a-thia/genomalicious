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
#' @param returnParents Logical: Should the parental genotypes also be returned?
#' Default is FALSE.
#'
#' @param returnPedigree Logical: Shoudl the pedgree also be returned?
#' Default is FALSE.
#'
#' @details In this function, a 'simulation' comprises a draw of 10 individuals:
#' 2 are unrelated to each other and all other individuals in the simulation,
#' 2 are a parent-offspring pair, 2 are siblings, 2 are half-siblings,
#' 2 are cousins, and 2 are half-cousins.
#'
#' The function always returns the genotypes for 'focal pairs' from a set of \code{numSim}
#' simulations. Each simulation takes the naming convention 'S[sim number]',
#' for example, S1, S2 ... S10 if 10 simulations are requested. The naming
#' convention for each simulated individual is 'S[sim number]_[relationship]',
#' with the 'relationship' tag being one of the following:
#' \itemize{
#'  \item 'unrel_ind_1' and 'unrel_ind_2' are two completely unrelated individuals.
#'  \item 'po_ind_1' and 'po_dam_1' are an offspring and its parent (mother).
#'  \item 'sib_ind_1' and 'sib_ind_2' are full siblings.
#'  \item 'halfsib_ind_1' are 'halfsib_ind_2' are half siblings.
#'  \item 'cous_ind_1' and 'cous_ind_2' are cousins.
#'  \item 'halfcous_ind_1' and 'halfcous_ind_2' are half cousins.
#' }
#'
#' Note, that within the sth simulation, all pairs are unrelated to each other as
#' well. For example, in simulation S1, sib_ind_1 and sib_ind_2 are completely unrelated
#' to halfsib_ind_1 and halfsib_ind_2. Within an sth simulation there are 66 pairs,
#' 5 of these are related (the focal parent-offspring, siblings, half-siblings,
#' cousins, and half-cousins), and the rest are unrelated.
#' Individuals from different simulation are also completely unrelated. For example,
#' all individuals from simulation S1 are completely unrelated to all individuals
#' from simulation S2.
#'
#' If requested with \code{returnParents==TRUE}, then the genotypes of all parental
#' indiviudals will also be returned for the parent-offspring, sibling, half-sibling,
#' cousin, and half-cousin focal pairs. This includes grandparents to generate the
#' cousin and half-cousin focal pairs.
#'
#' If requested with \code{returnPedigree==TRUE}, the pedigree will be returned.
#'
#' @returns Returns a list with indexes \code{$focal.pairs}, \code{$parents}, and
#' \code{$pedigree}.
#'
#' The \code{$focal.pairs} index will always contain a data.table with the following columns:
#' \enumerate{
#'    \item \code{$PLOIDY}, the ploidy for the locus.
#'    \item \code{$LOCUS}, the locus ID.
#'    \item \code{$FREQ}, the individual-level allele frequency.
#'    \item \code{$SIM}, the simulation number.
#'    \item \code{$SAMPLE}, the sample ID.
#'    \item \code{$FAMILY}, the familial relationship this individual was used to simualte.
#'    \item \code{$GT}, the genotype as counts of the Alt alleles.
#' }
#'
#' If parent genotypes are requested, the \code{$pedigree} index will contain
#' a data.table.
#'
#' If the pedigree is requested,
#' \enumerate{
#'    \item \code{$SAMPLE}, the sample ID.
#'    \item \code{$DAM}, the dam's ID.
#'    \item \code{$SIRE}, the sire's ID.
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
#' # Simulate 100 familial relationships of each class
#' simFamily <- family_sim_genos(
#'    freqData=freqsPop1,
#'    locusCol='LOCUS',
#'    freqCol='FREQ',
#'    numSims=100,
#'    returnParents=TRUE,
#'    returnPedigree=TRUE
#' )
#'
#' # Take a look at the focal pairs
#' simFamily$focal.pairs
#'
#' # Take a look at the parentals
#' simFamily$parents
#'
#' # Take a look at the pedigree
#' simFamily$pedigree
#'
#' ### THE OBSERVED GENOMIC RELATIONSHIPS MATRIX
#' library(AGHmatrix)
#'
#' # A genotype matrix for the focal pairs
#' obsGenosMat <- simFamily$focal.pairs %>% DT2Mat_genos()
#'
#' # Calculate the GRM
#' obsGRM <- Gmatrix(obsGenosMat, method='Yang', ploidy=2)
#'
#' ### THE SIMULATED GENOMIC RELATIONSHIPS MATRIX
#' # Convert simulated families into a genotype matrix
#' simGenosMat <- DT2Mat_genos(simFamily$focal.pairs)
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
#' # and siblings/parent-offspring, respectively.
#' # You will note a large variance in the expected values, which
#' # is not surprising for this very small SNP dataset (200 loci).
#' relComp$plot
#'
#' ### SIMULATE A QUANTITATIVE TRAIT
#'
#' # Combine the focal pairs and parentals, and simualte a trait controlled by
#' # 100 loci with Va = 1, and Ve = 1.
#' simQTL <- family_sim_qtl(
#'   famGenos=rbind(simFamily$focal.pairs, simFamily$parents),
#'   numLoci=100, additiveVar=1, environVar=1
#'   )
#'
#' # The trait values
#' simQTL$trait
#'
#' # The locus values
#' simQTL$loci
#'
#' @export

family_sim_genos <- function(
    freqData, locusCol='LOCUS', freqCol='FREQ', numSims=100L,
    returnParents=FALSE, returnPedigree=FALSE){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('tidyr', 'data.table')){require(lib, character.only=TRUE)}

  # Check if data.table/can be converted to data.table
  freqData <- as.data.table(freqData)
  if(!'data.table' %in% class(freqData)){
    stop('Argument `freqData` must be a data.table class. See ?sim_family_data.')
  }

  # Check that columns are specified correctly
  column.check <- sum(c(locusCol,freqCol) %in% colnames(freqData))
  if(column.check!=2){
    stop(
      'Arguments `locusCol` and `freqCol` must be column names in the argument
    freqData. See ?sim_family_data.')
  }

  # Reassign column names.
  freqData <- freqData %>%
    copy %>%
    setnames(., old=c(locusCol,freqCol), new=c('LOCUS','FREQ'))

  # Check that numSims is >0 and is an integer
  numSims <- as.integer(numSims)
  if(!numSims>0){
    stop('Argument `numSims` must be an integer >0. See ?sim_family_data.')
  }

  # Return values for parents and pedigree data
  if(class(returnParents)!='logical'){
    stop('Argument `returnParents` must be TRUE or FALSE. See ?sim_family_data.')
  }

  if(class(returnPedigree)!='logical'){
    stop('Argument `returnPedigree` must be TRUE or FALSE. See ?sim_family_data.')
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
  focPairsLs <- list()
  parentsLs <- list()
  pedigreeLs <- list()

  for(sim in 1:numSims){
    # Note, 'g[x]' indicates generation in each group.

    # Random completely unrelated individuals
    unrel_ind_1 <- FUN_draw_genos(freqData, ploidy=2)
    unrel_ind_2 <- FUN_draw_genos(freqData, ploidy=2)

    unrel.names <- paste0('S',sim,c('_unrel_ind_1','_unrel_ind_2'))

    # Create parent-offspring pair
    po_dam_1 <- FUN_draw_genos(freqData, ploidy=2)
    po_sire_1 <- FUN_draw_genos(freqData, ploidy=2)

    po_ind_1 <- FUN_make_diploid_offspring(po_dam_1, po_sire_1)

    po.names <- paste0('S',sim,c('_po_ind_1','_po_dam_1','_po_sire_1'))

    # Create sibling pair
    sib_dam_1 <- FUN_draw_genos(freqData, ploidy=2)
    sib_sire_1 <- FUN_draw_genos(freqData, ploidy=2)

    sib_ind_1 <- FUN_make_diploid_offspring(sib_dam_1, sib_sire_1)
    sib_ind_2 <- FUN_make_diploid_offspring(sib_dam_1, sib_sire_1)

    sib.names <- paste0('S',sim,c('_sib_ind_1','_sib_ind_2','_sib_dam_1','_sib_sire_1'))

    # Create half sibling pair
    hsib_dam_1 <- FUN_draw_genos(freqData, ploidy=2)
    hsib_sire_1 <- FUN_draw_genos(freqData, ploidy=2)
    hsib_sire_2 <- FUN_draw_genos(freqData, ploidy=2)

    hsib_ind_1 <- FUN_make_diploid_offspring(hsib_dam_1, hsib_sire_1)
    hsib_ind_2 <- FUN_make_diploid_offspring(hsib_dam_1, hsib_sire_2)

    hsib.names <- paste0(
      'S',sim,
      c('_hsib_ind_1','_hsib_ind_2',
        '_hsib_dam_1','_hsib_sire_1','_hsib_sire_2')
    )

    # Create cousin pair
    cuz_gdam_1 <- FUN_draw_genos(freqData, ploidy=2)
    cuz_gsire_1 <- FUN_draw_genos(freqData, ploidy=2)

    cuz_dam_1 <- FUN_make_diploid_offspring(cuz_gdam_1, cuz_gsire_1)
    cuz_dam_2 <- FUN_make_diploid_offspring(cuz_gdam_1, cuz_gsire_1)
    cuz_sire_1 <- FUN_draw_genos(freqData, ploidy=2)
    cuz_sire_2 <- FUN_draw_genos(freqData, ploidy=2)

    cuz_ind_1 <- FUN_make_diploid_offspring(cuz_dam_1, cuz_sire_1)
    cuz_ind_2 <- FUN_make_diploid_offspring(cuz_dam_2, cuz_sire_2)

    cuz.names <- paste0(
      'S',sim,
      c('_cuz_ind_1','_cuz_ind_2',
        '_cuz_dam_1','_cuz_dam_2','_cuz_sire_1','_cuz_sire_2',
        '_cuz_gdam_1','_cuz_gsire_1')
    )

    # Create half-cousin pair
    hcuz_gdam_1 <- FUN_draw_genos(freqData, ploidy=2)
    hcuz_gsire_1 <- FUN_draw_genos(freqData, ploidy=2)
    hcuz_gsire_2 <- FUN_draw_genos(freqData, ploidy=2)

    hcuz_dam_1 <- FUN_make_diploid_offspring(hcuz_gdam_1, hcuz_gsire_1)
    hcuz_dam_2 <- FUN_make_diploid_offspring(hcuz_gdam_1, hcuz_gsire_2)
    hcuz_sire_1 <- FUN_draw_genos(freqData, ploidy=2)
    hcuz_sire_2 <- FUN_draw_genos(freqData, ploidy=2)

    hcuz_ind_1 <- FUN_make_diploid_offspring(hcuz_dam_1, hcuz_sire_1)
    hcuz_ind_2 <- FUN_make_diploid_offspring(hcuz_dam_2, hcuz_sire_2)

    hcuz.names <- paste0(
      'S',sim,
      c('_hcuz_ind_1','_hcuz_ind_2',
        '_hcuz_dam_1','_hcuz_dam_2','_hcuz_sire_1','_hcuz_sire_2',
        '_hcuz_gdam_1','_hcuz_gsire_1','_hcuz_gsire_2')
    )

    # Focal pairs of individuals
    foc.pairs.tab <- rbind(
      unrel_ind_1 %>% data.table(SIM=sim, SAMPLE=unrel.names[1]),
      unrel_ind_2 %>% data.table(SIM=sim, SAMPLE=unrel.names[2]),
      po_ind_1 %>% data.table(SIM=sim, SAMPLE=po.names[1]),
      po_dam_1 %>% data.table(SIM=sim, SAMPLE=po.names[2]),
      sib_ind_1 %>% data.table(SIM=sim, SAMPLE=sib.names[1]),
      sib_ind_2 %>% data.table(SIM=sim, SAMPLE=sib.names[2]),
      hsib_ind_1 %>% data.table(SIM=sim, SAMPLE=hsib.names[1]),
      hsib_ind_2 %>% data.table(SIM=sim, SAMPLE=hsib.names[2]),
      cuz_ind_1 %>% data.table(SIM=sim, SAMPLE=cuz.names[1]),
      cuz_ind_2 %>% data.table(SIM=sim, SAMPLE=cuz.names[2]),
      hcuz_ind_1 %>% data.table(SIM=sim, SAMPLE=hcuz.names[1]),
      hcuz_ind_2 %>% data.table(SIM=sim, SAMPLE=hcuz.names[2])
    ) %>%
      data.table(
        ., FAMILY=c(
          rep('Unrelated',2), rep('Parent-offspring',2), rep('Siblings',2),
          rep('Half-siblings',2), rep('Cousins',2), rep('Half-cousins',2)
        )
      ) %>%
      .[, GT:=FREQ*PLOIDY]

    # All other dams, sires, grandams and grandsires.
    dam.sire.tab <- rbind(
      po_sire_1 %>% data.table(SIM=sim, SAMPLE=po.names[3], FAMILY='Parent-offspring'),
      sib_dam_1 %>% data.table(SIM=sim, SAMPLE=sib.names[3], FAMILY='Siblings'),
      sib_sire_1 %>% data.table(SIM=sim, SAMPLE=sib.names[4], FAMILY='Siblings'),
      hsib_dam_1 %>% data.table(SIM=sim, SAMPLE=hsib.names[3], FAMILY='Half-siblings'),
      hsib_sire_1 %>% data.table(SIM=sim, SAMPLE=hsib.names[4], FAMILY='Half-siblings'),
      hsib_sire_2 %>% data.table(SIM=sim, SAMPLE=hsib.names[5], FAMILY='Half-siblings'),
      cuz_dam_1 %>% data.table(SIM=sim, SAMPLE=cuz.names[3], FAMILY='Cousins'),
      cuz_dam_2 %>% data.table(SIM=sim, SAMPLE=cuz.names[4], FAMILY='Cousins'),
      cuz_sire_1 %>% data.table(SIM=sim, SAMPLE=cuz.names[5], FAMILY='Cousins'),
      cuz_sire_2 %>% data.table(SIM=sim, SAMPLE=cuz.names[6], FAMILY='Cousins'),
      cuz_gdam_1 %>% data.table(SIM=sim, SAMPLE=cuz.names[7], FAMILY='Cousins'),
      cuz_gsire_1 %>% data.table(SIM=sim, SAMPLE=cuz.names[8], FAMILY='Cousins'),
      hcuz_dam_1 %>% data.table(SIM=sim, SAMPLE=hcuz.names[3], FAMILY='Half-cousins'),
      hcuz_dam_2 %>% data.table(SIM=sim, SAMPLE=hcuz.names[4], FAMILY='Half-cousins'),
      hcuz_sire_1 %>% data.table(SIM=sim, SAMPLE=hcuz.names[5], FAMILY='Half-cousins'),
      hcuz_sire_2 %>% data.table(SIM=sim, SAMPLE=hcuz.names[6], FAMILY='Half-cousins'),
      hcuz_gdam_1 %>% data.table(SIM=sim, SAMPLE=hcuz.names[7], FAMILY='Half-cousins'),
      hcuz_gsire_1 %>% data.table(SIM=sim, SAMPLE=hcuz.names[8], FAMILY='Half-cousins'),
      hcuz_gsire_2 %>% data.table(SIM=sim, SAMPLE=hcuz.names[9], FAMILY='Half-cousins')
    ) %>%
      .[, GT:=FREQ*PLOIDY]

    # The pedigree
    ped.tab <- rbind(
      data.table(
        SAMPLE=rev(unrel.names),
        DAM=rep(NA,2), SIRE=rep(NA,2)
      ),
      data.table(
        SAMPLE=rev(po.names),
        DAM=c(NA,NA,po.names[2]),
        SIRE=c(NA,NA,po.names[3])
      ),
      data.table(
        SAMPLE=rev(sib.names),
        DAM=c(NA,NA,sib.names[3],sib.names[3]),
        SIRE=c(NA,NA,sib.names[4],sib.names[4])
      ),
      data.table(
        SAMPLE=rev(hsib.names),
        DAM=c(NA,NA,NA,hsib.names[3],hsib.names[3]),
        SIRE=c(NA,NA,NA,hsib.names[5],hsib.names[4])
      ),
      data.table(
        SAMPLE=rev(cuz.names),
        DAM=c(NA,NA,NA,NA,cuz.names[7],cuz.names[7],cuz.names[4],cuz.names[3]),
        SIRE=c(NA,NA,NA,NA,cuz.names[8],cuz.names[8],cuz.names[6],cuz.names[5])
      ),
      data.table(
        SAMPLE=rev(hcuz.names),
        DAM=c(NA,NA,NA,NA,NA,hcuz.names[7],hcuz.names[7],hcuz.names[4],hcuz.names[3]),
        SIRE=c(NA,NA,NA,NA,NA,hcuz.names[9],hcuz.names[8],hcuz.names[6],hcuz.names[5])
      )
    )

    # Add in where relevant
    focPairsLs[[sim]] <- foc.pairs.tab

    if(returnParents==TRUE){
      parentsLs[[sim]] <- dam.sire.tab
    }

    if(returnPedigree==TRUE){
      pedigreeLs[[sim]] <- ped.tab
    }
  }

  # Output
  list(
    focal.pairs=do.call('rbind',focPairsLs),
    parents=do.call('rbind',parentsLs),
    pedigree=do.call('rbind',pedigreeLs)
  ) %>% return
}
