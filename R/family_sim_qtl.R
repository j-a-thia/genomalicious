#' Simulate quantitative trait from simulated families
#'
#' This function is designed to pair with the \code{family_sim_genos} function
#' to simulate a quantitative trait from a set of simulated families.
#' However, it could also be used to simulate a quantitative trait from observed data.
#' The function is currently only designed to simulate a single trait controlled
#' by \emph{n} loci with phenotypic variance, \emph{Vp}, the sum of the
#' additive, \emph{Va}, and environmental, \emph{Ve}, variances.
#'
#' @param famGenos Data.table: Long-format, must contain the columns
#' \code{$SAMPLE} of sample IDs, \code{$LOCUS} of locus IDs, and \code{$GT} of
#' genotypes.
#'
#' @param numLoci Integer: The number of loci contributing to the trait. These loci
#' are drawn from random in the dataset, \code{famGenos}. Therefore, it is
#' essential the \code{numLoci} is <= the total unique loci in \code{famGenos}.
#' Alternate parameterisation is to specify the names of the loci to use with \code{qtlLoci}.
#'
#' @param qtlLoci Character: A character vector of loci contributing to the trait.
#' These loci must be present in the dataset, \code{famGenos}.
#' Alternate parameterisation is to specify the number of  loci to use with \code{numLoci}.
#'
#' @param additiveVar Numeric: The additive genetic variance.
#'
#' @param environVar Numeric: The environmental variance.
#'
#' @details
#' Note, the current implementation of this function generates a specific
#' additive + environmental variance for a trait, given a set of allele frequencies
#' and the number/identify of loci of interest. The phenotypes and genotypic values
#' simulated are only relevant to the input sample only.
#'
#' A total of \emph{n} == \code{numLoci} loci from \code{famGenos} are drawn
#' at random to underpin the quantitative trait. It is assumed that all loci contribute
#' equally to the trait, i.e., \emph{Va}/\emph{n}. The simulation starts by first
#' fitting genotype AA with a genetic value of +1, Aa with a value of 0,
#' and aa with a value of -1. Because the additive genetic variance
#' is also dependent on alelle frequencies: \cr\cr
#' \emph{Va} = 2\emph{pq}[\emph{a} + \emph{d}(\emph{q} - \emph{p})]^2 \cr\cr
#' The genetic values are then modified relative to the allele frequencies to ensure
#' the total additivie genetic variance in the sample sum to the specified \emph{Va}.
#' The phenotype per individual is the sum of their genotypic values plus randomly
#' drawn environmental deviation (\code{rnorm(.., mean=0, sd=sqrt(Ve))}).
#'
#' @returns
#' A list with the indexes \code{$trait} and \code{$loci} is returned, containing
#' information on the traits and the underpinning loci.
#'
#' The index \code{$trait} contains a data.table with the following columns:
#' \enumerate{
#'    \item \code{$SAMPLE}, the sample IDs.
#'    \item \code{$G}, the additive genetic values, with variance \emph{Va}.
#'    \item \code{$E}, the environmental values, with variance \emph{Ve}.
#'    \item \code{$P}, the phenoypic values, \emph{G} + \emph{E}.
#' }
#'
#' The index \code{$loci} contains a data.table with the following columns:
#' \enumerate{
#'    \item \code{$LOCUS}, the locus IDs.
#'    \item \code{$FREQ}, the frequency of the focal allele.
#'    \item \code{$A.VAL}, the additive genetic value for this locus.
#' }
#'
#' @references Falconer and MacKay (1996) Introduction to Quantitative Genetics, 3rd Ed., chapter 8, page 129.
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
#' obsGenosMat <- genosPop1 %>% DT2Mat_genos()
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
family_sim_qtl <- function(famGenos, numLoci=NULL, qtlLoci=NULL, additiveVar, environVar){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+

  # Checks for famGenos
  req.cols <- c('GT','LOCUS','SAMPLE')

  if(!'data.table' %in% class(famGenos)){
    stop('Argument `famGenos` must be a data.table class. See ?family_sim_qtl.')
  }

  if(sum(req.cols %in% colnames(famGenos))!=3){
    stop('Not all required columns are in argument `famGenos`. See ?family_sim_qtl.')
  }

  # Make sure only one of numLoci and qtlLoci are specified
  if(is.null(numLoci) & is.null(qtlLoci)){
    stop('Neither arugment `numLoci` or `qtlLoci` has been specified, but you must specify one of these. See ?family_sim_qtl.')
  }

  if(!is.null(numLoci) & !is.null(qtlLoci)){
    stop('Both arguments `numLoci` and `qtlLoci` have been specified, but you must only specifcy one of these. See ?family_sim_qtl.')
  }

  # Make sure that numLoci is <= the total number of loci
  loc.uniq <- unique(famGenos$LOCUS)

  if(!is.null(numLoci)){
    if(numLoci>length(loc.uniq) ){
      stop('Argument `numLoci` must be <= the total number of unique loci in `famGenos$LOCUS`. See ?family_sim_qtl.')
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Unique loci and their frequencies
  uniqLocTab <- unique(famGenos[,.(FREQ=sum(GT)/(length(GT)*2)),LOCUS])

  # Reassign values for variance components
  Va_t <- additiveVar
  Ve_t <- environVar

  # Subset a set of loci to underpin the trait
  if(is.null(numLoci)==FALSE){
    n <- numLoci
    traitLocTab <- uniqLocTab[sample(1:nrow(uniqLocTab), size=n, replace=FALSE),]
  } else if(is.null(qtlLoci)==FALSE){
    traitLocTab <- uniqLocTab[LOCUS %in% qtlLoci]
    n <- length(qtlLoci)
  }

  # Frequency of trait loci
  p_vec <- traitLocTab$FREQ
  q_vec <- 1 - p_vec

  # The effect of the AA genotype (a) and the Aa genotype (d)
  a <- 1
  d <- 0

  # Falconer and MacKay: Introduction to Quantitative genetics, 4th Ed., page 129
  # Va == 2pq[a + d(q-p)]^2

  # An initial calculation of the variance per locus. This won't sum to the
  # desired amount, so have to scale the per locus variance
  Va_i_vec <- 2 * p_vec * q_vec * (a + d * (q_vec - p_vec))^2
  sum(Va_i_vec)

  Va_i_vec <- Va_i_vec * (Va_t / sum(Va_i_vec))
  sum(Va_i_vec)

  # We also need to scale the AA effect per locus after scaling the variance
  # per locus. All loci should have the same effect, unless the locus is
  # fixed for one of the two alleles.
  a_vec <- lapply(1:n, function(i){
    # If: Va_i == 2pq[a + d(q - p)]^2, and if d = 0, then the solution is
    # Va_i == 2pq[a]^2, so a == sqrt(Va_i/2pq)

    v_i <- Va_i_vec[i]
    p_i <- p_vec[i]
    q_i <- q_vec[i]

    a_i <- sqrt(v_i / (2*p_i*q_i))

    if_else(is.na(a_i), 0, a_i)
  }) %>% unlist()

  # Using the scaled effects
  sum(2 * p_vec * q_vec * (a_vec + d * (q_vec - p_vec))^2)

  # Add into the table of trait loci
  traitLocTab$A.VAL <- a_vec

  # Combine the trait loci into the genotype table and calculate the
  # additive genotypic effects
  simFamTraits <- left_join(famGenos, traitLocTab[, c('LOCUS','A.VAL')]) %>%
    # Remove loci with no contributions
    .[!is.na(A.VAL)] %>%
    # If genotype is 0, breeding value is +a.
    .[GT==0, G:=A.VAL] %>%
    # If genotype is 1, breeding value is 0.
    .[GT==1, G:=0] %>%
    # If genotype is 2, breeding value is -a.
    .[GT==2, G:=-A.VAL] %>%
    # Sum the breeding values
    .[, .(G=sum(G)), by=SAMPLE]

  simFamTraits$G %>% var

  # Add in an environmental deviation per individual.
  simFamTraits$E <- rnorm(nrow(simFamTraits), mean=0, sd=sqrt(Ve_t))
  simFamTraits$E <- simFamTraits$E * (sqrt(Ve_t) / sd(simFamTraits$E))

  simFamTraits$E %>% var

  # The phenotype as a combination of genotypic and environmental effects
  simFamTraits[, P:=G+E, by=SAMPLE]

  # simFamTraits$P %>%  var

  # simFamTraits$P %>% hist

  # Correlations between genotype and environmental contributions to phenotype
  # pairs(simFamTraits[, c('G', 'E','P')] )

  # Output
  list(trait=simFamTraits, loci=traitLocTab) %>%
  return()
}
