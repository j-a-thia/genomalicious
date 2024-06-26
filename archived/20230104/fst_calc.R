#' Calculate Weir and Cockerham's FST from genotypes or allele frequencies
#'
#' Function takes a genotypes or allele frequencies in a long-format data table
#' and calculates Weir & Cockerhams FST, i.e. theta (Weir & Cockerham, 1984).
#' Permutations can be used to test statistical significance of FST in
#' genotype data sets.
#'
#' @param dat Data table: A long format data table of biallelic genotypes,
#' coded as '/' separated alleles ('0/0', '0/1', '1/1'), or counts
#' of the Alt alleles (0, 1, 2, respectively).
#' Columns required for \strong{both} genotypes and allele frequencies:
#' \enumerate{
#'    \item The population ID (see param \code{popCol}).
#'    \item The locus ID (see param \code{locusCol}).
#' }
#' Columns required only for \strong{genotypes}:
#' \enumerate{
#'    \item The sample ID (see param \code{sampCol}).
#'    \item The genotypes (see param \code{genoCol}).
#' }
#' Columns required only for \strong{allele frequencies}:
#' \enumerate{
#'    \item The allele frequencies (see param \code{freqCol}).
#'    \item The number of individuals used to obtain the allele frequency
#'    estimate (see param \code{indsCol}).
#' }
#'
#' @param type Character: One of \code{'genos'} or \code{'freqs'}, to calculate
#' FST on genotype or allele frequency data, respectively.
#'
#' @param popCol Character: The column name with the population information.
#' Default is \code{'POP'}.
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
#' @param freqCol Character: The column name with the allele freuqency information.
#' Default is \code{'FREQ'}.
#'
#' @param indsCol Character: The column name with the number of individuals
#' contributing to the allele freuqency estimate. Default is \code{indsCol}.
#'
#' @param global Logical: Calculate the global FST across all populations?
#' Default is \code{TRUE}. Note, \code{global==TRUE} and \code{pairwise==TRUE}
#' is not permitted.
#'
#' @param pairwise Logical: Calculate pairwise FST between populations?
#' Default is \code{FALSE}. Note, \code{global==TRUE} and \code{pairwise==TRUE}
#' is not permitted.
#'
#' @param permute Logical: Should permutations be performed to test the
#' statistical significance of FST? Default is \code{FALSE}. Can only be performed
#' on genotype data, i.e., \code{type=='genos'}.
#'
#' @param numIters Integer: The number of permutations to perform. Default is 100.
#'
#' @details Permutation tests involve random shuffling of individuals
#' among populations, recalculating FST, and testing the hypothesis
#' that the permuted FST > observed FST. The p-value represents the
#' proportion of iterations that were "TRUE" to this expression. That is,
#' if no permuted values are greater than the observed, p=0. Likewise, if all
#' the permuted values are greater than the observed, p=1.
#'
#' @return A list is returned with three indexes.
#'
#' The first index is \code{$locus}, the locus-specific FST. This is a data table
#' object with the columns \code{$LOCUS}, the locus ID, \code{$NUMER}, the
#' numerator variance component, \code{$DENOM}, the denominator variance component,
#' and \code{$FST}, the locus specific FST. Additionally, if pairwise population
#' estimates have been requested, \code{pairwise==TRUE}, then there are
#' \code{$POP1} and \code{$POP2}, which represent the two populations tested.
#'
#' The second index is \code{$genome}, the genome-wide FST. If the global estimate
#' has been requested, \code{global==TRUE}, then this is just a single value;
#' the estimate across all populations. If pairwise esimates were requested,
#' \code{pairwise==TRUE}, then there are \code{$POP1} and \code{$POP2}, which
#' represent two populations tested.
#'
#' The third index is \code{$permute}, the permutation results. This index will
#' be \code{NULL} when frequencies are used, i.e., \code{type=='freqs'}, and
#' will only contain data if \code{type=='genos'} and \code{permute==TRUE}.
#' \code{$permute} is itself a list, with two subindexes:
#' \enumerate{
#'    \item \code{$fst}: The permuted FST values. If \code{global==TRUE}, then
#'    this will simply be a numeric vector of permuted global FST values.
#'    If \code{pairwise==TRUE}, then this will be a data table with columns
#'    \code{$POP1}, \code{$POP2}, and \code{$FST}.
#'
#'    \item \code{$pval}: The permuted p-values. If \code{global==TRUE}, then
#'    this will single p-value for the global FST. If \code{pairwise==TRUE},
#'    then this will be a data table with columns \code{$POP1}, \code{$POP2},
#'    and \code{$PVAL}.
#' }
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#' data(data_PoolFreqs)
#'
#' ### Genotypes, global, with 50 permutations
#' genoGlobal <- fst_calc(
#'    data_Genos,
#'    type='genos',
#'    global=TRUE,
#'    permute=TRUE,
#'    numIters=50)
#'
#' # Locus specific estimates
#' genoGlobal$locus
#'
#' # Genome-wide estimate
#' genoGlobal$genome
#'
#' # Permutations
#' genoGlobal$permute$fst
#' genoGlobal$permute$pval
#'
#' ### Genotypes, pairwise, no permutations
#' genoPairs <- fst_calc(
#'    data_Genos,
#'    type='genos',
#'    global=FALSE,
#'    pairwise=TRUE
#'    )
#'
#' # Locus specific estimates
#' genoPairs$locus
#'
#' # Genome-wide estimates
#' genoPairs$genome
#'
#' # Permutations (there are none!)
#' genoPairs$permute
#'
#' ### Genotypes, pairwise, with 10 permutations
#' genoPairsPerms <- fst_calc(
#'    data_Genos,
#'    type='genos',
#'    global=FALSE,
#'    pairwise=TRUE,
#'    permute=TRUE,
#'    numIters=10
#'    )
#'
#' # Permuted FST
#' genoPairsPerms$permute$fst
#'
#' # Permuted p-value
#' genoPairsPerms$permute$pval
#'
#' ### Frequencies, pairwise, with 10 permutations
#' # Note columns in data_PoolFreqs
#' colnames(data_PoolFreqs)
#'
#' # There is no $POP column, but there is a $POOL column that contains
#' # the population pool ID. Need to manually specify this.
#' freqPairs <- fst_calc(
#'    data_PoolFreqs,
#'    type='freqs',
#'    global=FALSE,
#'    pairwise=TRUE,
#'    popCol='POOL'
#' )
#'
#' @export

fst_calc <- function(
  dat, type, popCol='POP', sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT',
  freqCol='FREQ', indsCol='INDS', global=TRUE, pairwise=FALSE,
  permute=FALSE, numIters=100){

  # --------------------------------------------+
  # Assertions and environment
  # --------------------------------------------+
  require(data.table); require(tidyverse)

  # Make sure dat is a data.table
  dat <- as.data.table(dat)

  # Has type been specified correctly?
  if(!type %in% c('genos','freqs')){
    stop('Argument `type` must be one of "genos" or "freqs".')
  }

  # Check that global and pairwise are logicals
  if(class(global)!='logical'){
    stop('Argument `global` must be a logical value. See ?fst_calc.')
  }

  if(class(pairwise)!='logical'){
    stop('Argument `pairwise` must be a logical value. See ?fst_calc.')
  }

  # Cannot specify bpth global and pairwise
  if(global==TRUE & pairwise==TRUE){
    stop('Arguments `global` and `pairwise` cannot both be TRUE. See ?fst_calc.')
  }

  # Must specify one of global or pairwise
  if(global==FALSE & pairwise==FALSE){
    stop('One of arguments `global` and `pairwise` must be TRUE. See ?fst_calc.')
  }

  # Check that permute is logical
  if(class(permute)!='logical'){
    stop('Argument `permute` must be a logical value. See ?fst_calc.')
  }

  # If permute is specified, then check numIters is >0 and is an integer
  if(permute==TRUE){
    if(numIters<1){
      stop('Argument `numIters` must be an integer >0. See ?fst_calc.')
    }

    numIters <- as.integer(numIters)
  }

  # Check all necessary columns are present and reassign values.
  # If genotypes are character, convert to integers.
  if(type=='genos'){
    column.orig <- c(popCol,sampCol,locusCol,genoCol)

    if(sum(column.orig %in% colnames(dat))!=4){
      stop(
        'Arguments `popCol`, `sampCol`, `locusCol`, and `genoCol` must all be
      column names in `dat`. See ?fst_calc.')
    }

    col.index <- match(column.orig,colnames(dat))
    colnames(dat)[col.index] <- c('POP','SAMPLE','LOCUS','GT')

    if(class(dat$GT)=='character'){
      dat$GT <- genoscore_converter(dat$GT)
    }
  }

  if(type=='freqs'){
    column.orig <- c(popCol,locusCol,freqCol,indsCol)

    if(sum(column.orig %in% colnames(dat))!=4){
      stop(
        'Arguments `popCol`, `locusCol`, `freqCol`, and `indsCol`
        must all be column names in `dat`. See ?fst_calc.')
    }

    col.index <- match(column.orig,colnames(dat))
    colnames(dat)[col.index] <- c('POP','LOCUS','FREQ','INDS')
  }

  # --------------------------------------------+
  # Internal functions
  # --------------------------------------------+

  # Variance components for FST from frequencies
  FUN_varcomps_freqs <- function(pi, ni, r){
    # Mean weighted allele frequency
    p.mean <- sum(ni * pi)/sum(ni)

    # Mean squares variance components
    msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
    msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

    # Sample size correction factor
    nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

    # Output as data.table
    data.table(NUMER=msp-msg, DENOM=msp+(nc-1)*msg)
  }

  FUN_varcomps_genos <- function(ni, pi, hi, r){
    # The mean sample size
    n.mean <- sum(ni/r)

    # The sample size scaling parameter
    nc <- (r*n.mean - sum((ni^2)/(r*n.mean))) / (r-1)

    # The average sample allele frequency
    p.mean <- sum((ni*pi)/(r*n.mean))

    # The variance in allele frequencies
    s2 <- sum( (ni*(pi-p.mean)^2)/((r-1)*n.mean) )

    # The average heterozygosity
    h.mean <- sum( (ni*hi)/(r*n.mean) )

    # The a, b, and c components
    a <- (nc/n.mean) * (s2 - (1/(n.mean-1))*((p.mean*(1-p.mean)) - (s2*(r-1)/r) - (0.25*h.mean)))

    b <- (n.mean/(n.mean-1)) * ((p.mean*(1-p.mean)) - (s2*(r-1)/r) - (h.mean*((2*n.mean-1)/(4*n.mean))))

    c <- 0.5 * h.mean

    # Return as numerator and denominator
    return(data.table(NUMER=a, DENOM=a+b+c))
  }

  # Variance components for FST from genotypes
  FUN_fst_genos <- function(Dx){
    left_join(
      Dx[, .(FREQ=sum(GT)/(length(GT)*2)), by=c('POP','LOCUS')],
      Dx[, .(HO=sum(GT==1)/length(GT)), by=c('POP','LOCUS')],
      by=c('POP','LOCUS')
    ) %>%
      left_join(
        ., Dx[, .(N=length(unique(SAMPLE))), by=c('POP','LOCUS')], by=c('POP','LOCUS')
      ) %>%
      .[, FUN_varcomps_genos(ni=N, pi=FREQ, hi=HO, r=length(unique(POP))), by=LOCUS] %>%
      .[, FST:=NUMER/DENOM]
  }

  # Permute FST
  FUN_fst_perm <- function(Dx, fst.genome, type, numIters){
    # Dx is a data table with $POP, $SAMPLE, $LOCUS, $GT.

    # fst.genome is the genome-wide FST to test.

    # numIters is the number of iterations to perform.

    perm.list <- list()

    # Create a table to reshuffle populations
    pop.perm.tab <- Dx[, c('POP','SAMPLE')] %>% unique

    perm.results <- list()

    # Iterate through i iterations
    for(i in 1:numIters){
      # Reshuffle the populations in reference tbale
      pop.perm.tab$POP <- sample(pop.perm.tab$POP, nrow(pop.perm.tab), replace=FALSE)

      # Permute
      perm.results[[i]] <- Dx %>%
        # Subset the samples and their genotypes
        .[, c('SAMPLE','LOCUS','GT')] %>%
        # Add in the reshuffled populations
        left_join(., pop.perm.tab, by='SAMPLE') %>%
        # Calculate FST
        FUN_fst_genos(.) %>%
        .[, sum(NUMER)/sum(DENOM)]

      rm(i)
    }
    # Results
    perm.list$fst <- unlist(perm.results)
    perm.list$pval <- sum(perm.list$fst>fst.genome)/numIters
    return(perm.list)
  }

  # --------------------------------------------+
  # Code for genotypes
  # --------------------------------------------+

  if(type=='genos'){
    # START: FST among all populations
    if(global==TRUE){
      cat('FST calculation on genotype data, global estimate', '\n')
      fstOutput <- list()

      # Per locus and genome-wide
      fstOutput$locus <- FUN_fst_genos(dat)
      fstOutput$genome <- fstOutput$locus[, sum(NUMER)/sum(DENOM)]

      # For genotype permutations
      if(permute == TRUE){
        cat('Performing permutations', '\n')
        fstOutput$permute <- FUN_fst_perm(
          Dx=dat, fst.genome=fstOutput$genome, type='genos', numIters=numIters)
      } else{
        fstOutput$permute <- NULL
      }
      # END: FST among all populations
    } else if(pairwise==TRUE){
      # START: FST between pairwise populations
      cat('FST calculation on genotype data, pairwise estimates', '\n')
      pair.comp <- CJ(POP1=unique(dat$POP), POP2=unique(dat$POP)) %>%
        .[POP1!=POP2]
      num.pairs <- nrow(pair.comp)

      fstOutput <- list()

      fstOutput$locus <- list()
      fstOutput$genome <- list()
      fstOutput$permute <- list(fst=list(), pval=list())

      # START: Iterate over i pairs
      for(i in 1:nrow(pair.comp)){
        pop1 <- pair.comp$POP1[i]
        pop2 <- pair.comp$POP2[i]

        fstOutput$locus[[i]] <- dat[POP %in% c(pop1, pop2)] %>%
          FUN_fst_genos(.) %>%
          data.table(POP1=pop1, POP2=pop2, .)
        fstOutput$genome[[i]] <- fstOutput$locus[[i]] %>%
          .[, sum(NUMER)/sum(DENOM)] %>%
          data.table(POP1=pop1, POP2=pop2, FST=.)

        # If permuting
        if(permute==TRUE){
          cat('Performing permutation:', i, '/', num.pairs, '\n')

          pair.perm <- FUN_fst_perm(
            Dx=dat[POP %in% c(pop1, pop2)],
            fst.genome=fstOutput$genome[[i]]$FST,
            type='genos',
            numIters=numIters
          )

          fstOutput$permute$fst[[i]] <- pair.perm$fst %>%
            data.table(POP1=pop1, POP2=pop2, FST=.)

          fstOutput$permute$pval[[i]] <- data.table(
            POP1=pop1, POP2=pop2, PVAL=pair.perm$pval
          )
        } else{
          fstOutput$permute[[i]] <- NULL
        }
        rm(i)

        # END: Iterate over i pairs
      }

      # Combine results
      fstOutput$locus <- do.call('rbind', fstOutput$locus)

      fstOutput$genome <- do.call('rbind', fstOutput$genome)

      if(!is.null(fstOutput$permute$fst)){
        fstOutput$permute$fst <- do.call('rbind', fstOutput$permute$fst)
      }

      if(!is.null(fstOutput$permute$pval)){
        fstOutput$permute$pval <- do.call('rbind', fstOutput$permute$pval)
      }

      # END: FST between population pairs
    }
  }

  # --------------------------------------------+
  # Code for frequencies
  # --------------------------------------------+
  if(type=='freqs'){
    # START: FST among all populations
    if(global==TRUE){
      cat('FST calculation on frequency data, global estimate', '\n')
      fstOutput <- list()

      # Per locus and genome-wide
      fstOutput$locus <- dat %>%
        .[, FUN_varcomps_freqs(pi=FREQ, ni=INDS, r=length(unique(POP))), by=LOCUS] %>%
        .[, FST:=NUMER/DENOM]
      fstOutput$genome <- fstOutput$locus[, sum(NUMER)/sum(DENOM)]

      # No pemutations
      fstOutput$permute <- NULL

      # END: FST among all populations

    } else if(pairwise==TRUE){
      # START: FST between pairwise populations

      cat('FST calculation on frequency data, global estimate', '\n')

      pair.comp <- CJ(POP1=unique(dat$POP), POP2=unique(dat$POP)) %>%
        .[POP1!=POP2]
      num.pairs <- nrow(pair.comp)

      fstOutput <- list()

      fstOutput$locus <- list()
      fstOutput$genome <- list()

      # START: Iterate over i population pairs
      for(i in 1:nrow(pair.comp)){
        pop1 <- pair.comp$POP1[i]
        pop2 <- pair.comp$POP2[i]

        fstOutput$locus[[i]] <- dat[POP %in% c(pop1, pop2)] %>%
          .[, FUN_varcomps_freqs(pi=FREQ, ni=INDS, r=length(unique(POP))), by=LOCUS] %>%
          .[, FST:=NUMER/DENOM] %>%
          data.table(POP1=pop1, POP2=pop2, .)
        fstOutput$genome[[i]] <- fstOutput$locus[[i]] %>%
          .[, sum(NUMER)/sum(DENOM)] %>%
          data.table(POP1=pop1, POP2=pop2, FST=.)
        # END: Iteration over i population pairs
      }

      # Combine results
      fstOutput$locus <- do.call('rbind', fstOutput$locus)
      fstOutput$genome <- do.call('rbind', fstOutput$genome)
      fstOutput$permute <- NULL

      # END: FST between pairwise populations
    }
  }

  # --------------------------------------------+
  # Output results
  # --------------------------------------------+
  return(fstOutput)
}

