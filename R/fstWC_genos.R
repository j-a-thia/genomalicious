#' Calculate Weir and Cockerham's FST from genotypes
#'
#' Function takes a long format genotype data table and calculates
#' Weir & Cockerhams FST, i.e. theta (Weir & Cockerham, 1984). Can
#' also conduct permutation testing and calculate bootstrap
#' confidence intervals.
#'
#' @param dat Data table: A long format data table of biallelic genotypes,
#' coded as '/' separated alleles ('0/0', '0/1', '1/1'), or counts
#' of the Alt alleles (0, 1, 2, repsectively).
#' Four columns are required:
#' \enumerate{
#'    \item The population ID (see param \code{popCol}).
#'    \item The sampled individual ID (see param \code{sampCol}).
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The genotype (see param \code{genoCol}).
#' }
#'
#' @param popCol Character: The column name with the population information.
#' Default = \code{'POP'}.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default = \code{'GT'}.
#'
#' @param iters Integer: The number of permutations or bootstrap replicates.
#' Default = 100.
#'
#' @param permTest Logical: Should a permutation test be conducted?
#' Default = \code{FALSE}.
#'
#' @param bootCI Logical: Should bootstrap confidence intervals be estimated?
#' Default = \code{FALSE}.
#'
#' @param doPairs Logical: Should pairwise FST be calculated? Default = \code{FALSE},
#' which calculates the among population FST.
#'
#' @param doDist Logical: Should a a distance matrix of FST be returned?
#' Default = \code{FALSE}. Only applied when \code{doPairs==TRUE}.
#'
#' @param perLocus Logical: Should the per locus FST be returned?
#' Default = \code{FALSE}.
#'
#' @details Permutation tests involve random shuffling of individuals
#' among populations, recalculating FST, and testing the hypothesis
#' that the permuted FST > observed FST. The p-value represents the
#' proportion of iterations that were "TRUE" to this expression. \cr \cr
#' Bootstrapping randomly samples the loci with replacement
#' to estimate the distribution of FST values around the multilocus
#' FST. Confidence intervals are calculated as the 2.5 and
#' 97.5 percentiles of the obtained bootstrapped FST values.
#'
#' @return Returns a list, the contents vary depending on argument choices. \cr
#' When \code{doPairs==FALSE}, the analysis is among populations.
#' \enumerate{
#'    \item \code{$multilocus.mean}: the multilocus FST.
#'    \item \code{$multilocus.perms}: a vector of permuted FST values.
#'    \item \code{$multilocus.pval}: the p-value from permutation tests.
#'    \item \code{$multilocus.ci}: a vector of bootstrap confidence intervals.
#'    \item \code{$perlocus}: a data table of the FST for each locus.
#' }
#' When \code{doPairs==TRUE}, the analysis is between population pairs.
#' \enumerate{
#'    \item \code{$multilocus.mean}: a data table of multilocus FST.
#'    \item \code{$multilocus.ci}: a data table of bootstrap confidence intervals.
#'    \item \code{$multilocus.dist}: a distance matrix of multilocus FST values.
#'    \item \code{$perlocus}: a data table of the FST for each locus.
#' }
#'
#' @examples
#' data(data_4pops)
#'
#' # Among popuation FST, with permutation tests
#' # and bootstrap confidence intervals.
#' fst_among <- fstWC_genos(data_4pops
#'                 , permTest=TRUE
#'                 , bootCI=TRUE
#'                 , iters=50)
#'
#' # Pairwise FST, with per locus estimates
#' # and a distance matrix
#' fst_pairs <- fstWC_genos(data_4pops
#'                 , doPairs=TRUE
#'                 , doDist=TRUE
#'                 , perLocus=TRUE)
#'
#' @export
fstWC_genos <- function(dat, popCol='POP', sampCol='SAMPLE'
                        , locusCol='LOCUS', genoCol='GT'
                        , permTest=FALSE, bootCI=FALSE, iters=100
                        , doPairs=FALSE, doDist=FALSE, perLocus=FALSE){
  # BEGIN ...........

  # --------------------------------------------+
  # Assertions and environment
  # --------------------------------------------+
  # Check the class of dat
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Check that all columns are specified correctly
  if(sum(c(popCol, sampCol, locusCol, genoCol) %in% colnames(dat)) != 4){
    stop('Arguments popCol, sampCol, locusCol, and genoCol must be
         columns in dat. See ?genos2freqs.')
  }

  # Rename columns
  colnames(dat)[
    match(c(sampCol, popCol, locusCol, genoCol), colnames(dat))
    ] <- c('SAMPLE', 'POP', 'LOCUS', 'GT')

  # Get genotype class
  genoClass <- class(dat$GT)

  # Convert to integer
  if(genoClass=='character'){
    dat[, GT:=genoscore_converter(GT)]
  } else if(genoClass=='numeric'){
    dat[, GT:=as.integer(GT)]
  }

  # --------------------------------------------+
  # Code: Set up key parameters
  # --------------------------------------------+
  # Get the unique population and sample combinations
  popsampDt <- unique(dat[, c('POP', 'SAMPLE')])

  # Create a genotype matrix
  genoMat <- DT2Mat_genos(dat, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT')

  # Match genoMat rows to populations in popsampDt to
  # create a population ID
  pop_id <- popsampDt$POP[match(rownames(genoMat), popsampDt$SAMPLE)]

  # Create an empty list to store values in
  fstList <- list()
  fstList <- list()

  # --------------------------------------------+
  # Code: FST among all populations
  # --------------------------------------------+
  if(doPairs==FALSE){
    # START AMONG

    # Observed variance components
    obsVarcomps <- fstWC_varcomps(dat=genoMat, input_type='genos', pop_id=pop_id)

    # The observed FST
    fstList$multilocus.mean <- sum(obsVarcomps$NUMER) / sum(obsVarcomps$DENOM)

    if(permTest==TRUE){
      fstList$multilocus.perms <- fstWC_perm(genoMat=genoMat, pop_id=pop_id, perms=iters)

      # p-value for permutations
      fstList$multilocus.pval <- sum(fstList$multilocus.perms > fstList$multilocus.mean) / iters
    }

    if(bootCI==TRUE){
      fstList$multilocus.boots <- fstWC_boot(dat=genoMat, pop_id=pop_id, boots=iters)

      # Percentiles
      fstList$multilocus.ci <- quantile(fst_boot, c(0.025, 0.975))
    }

    if(perLocus==TRUE){
      obsVarcomps[, FST:=NUMER/DENOM]
      fstList$perlocus <- obsVarcomps[, c('LOCUS', 'FST')]
    }

    return(fstList)
    # END AMONG
  }

  # --------------------------------------------+
  # Code: FST between population pairs
  # --------------------------------------------+
  if(doPairs==TRUE){
    # START PAIRWISE

    # Unique populations
    uniq_pops <- unique(pop_id)

    # Population pair combinations
    pairCombos <- combn(uniq_pops, 2)

    # Get the pairwise variance components
    cat('Performing pairwise FST calculations', '\n')
    pairCalcs <- apply(pairCombos, 2, function(pops){
      cat('....', pops, '\n')
      # Create a list to return
      popList <- list()

      # Get the population ID indices
      pop1 <- which(pop_id==pops[1])
      pop2 <- which(pop_id==pops[2])

      # Pair name
      pair_id <- paste(pops[1], pops[2], sep='/')

      # Subset the genotypes and ID
      genoSub <- genoMat[c(pop1, pop2),]
      pop_id_sub <- pop_id[c(pop1, pop2)]

      # Calculate the variance components
      cat('........ Variance components', '\n')
      popsVarcomps <- fstWC_varcomps(dat=genoSub, input_type='genos', pop_id=pop_id_sub)

      popsVarcomps$POP1 <- pops[1]
      popsVarcomps$POP2 <- pops[2]
      popsVarcomps$PAIR <- pair_id

      # Multilocus FST
      popList$multilocus.mean <- data.table(
        popsVarcomps[1, c('POP1', 'POP2', 'PAIR')]
        , FST=sum(popsVarcomps$NUMER)/sum(popsVarcomps$DENOM))

      # Permutation tests
      if(permTest==TRUE){
        cat('........ Permutation tests', '\n')
        boot_fst <- fstWC_perm(genoMat=genoSub, pop_id=pop_id_sub, perms=iters)
        boot_pval <- sum(boot_fst > popList$multilocus.mean) / iters
        popList$multilocus.perms <- data.table(PAIR=pair_id, POP1=pops[1], POP2=pops[2], PERM=1:iters, FST=boot_fst)
        popList$multilocus.pval <- data.table(PAIR=pair_id, POP1=pops[1], POP2=pops[2], PVAL=boot_pval)
      }

      # Bootstrap confidence intervals
      if(bootCI==TRUE){
        cat('........ Bootstrap confidence intervals', '\n')
        boot_fst <- fstWC_boot(dat=genoSub, pop_id=pop_id_sub, boots=iters)
        popList$multilocus.boots <- data.table(PAIR=pair_id, POP1=pops[1], POP2=pops[2], BOOT=1:iters, FST=boot_fst)
        boot_ci <- quantile(boot_fst, c(0.025, 0.975))
        popList$multilocus.ci <- data.table(PAIR=pair_id, POP1=pops[1], POP2=pops[2], CI025=boot_ci[1], CI975=boot_ci[2])
      }

      # Calculate perLocus FST
      if(perLocus==TRUE){
        cat('........ Per locus FST', '\n')
        popsVarcomps[, FST:=NUMER/DENOM]
        popList$perlocus <- popsVarcomps[, c('PAIR', 'POP1', 'POP2', 'FST')]
      }

      # Return all the calculations for this population pair
      return(popList)
    })

    # Combine multilocus FST pairwise data
    fstList$multilocus.mean <- do.call('rbind', lapply(pairCalcs, function(pops){ pops$multilocus.mean }))

    # If permutation test conducted, combine pairwise data
    if(doPairs==TRUE){
      fstList$multilocus.perms <- do.call('rbind', lapply(pairCalcs, function(pops){ pops$multilocus.perms }))
      fstList$multilocus.pvals <- do.call('rbind', lapply(pairCalcs, function(pops){ pops$multilocus.pval }))
    }

    # If bootstrapped CI calculated, combine pairwise data
    if(bootCI==TRUE){
      fstList$multilocus.boots <- do.call('rbind', lapply(pairCalcs, function(pops){ pops$multilocus.boots }))
      fstList$multilocus.ci <- do.call('rbind', lapply(pairCalcs, function(pops){ pops$multilocus.ci }))
    }

    # If a distance matrix is required, make it.
    if(doDist==TRUE){
      distMat <- matrix(0, nrow=length(uniq_pops), ncol=length(uniq_pops), dimnames=list(uniq_pops, uniq_pops))
      for(j in 1:ncol(pairCombos)){
        pop1 <- pairCombos[1,j]
        pop2 <- pairCombos[2,j]
        distMat[pop2, pop1] <- fstList$multilocus.mean[POP1==pop1 & POP2==pop2]$FST
        fstList$multilocus.dist <- as.dist(distMat, diag=TRUE)
      }
    }

    # If per locus FST estimates calculated, combine pairwise data
    if(perLocus==TRUE){
      fstList$perlocus <- do.call('rbind', lapply(pairCalcs, function(pops){ pops$perlocus }))
    }

    return(fstList)
    # END PAIRWISE
  }

  # ............ END
  }
