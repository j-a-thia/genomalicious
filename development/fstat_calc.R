#' Calculate Weir and Cockerham's FST from genotypes or allele frequencies
#'
#' Function takes a genotypes or allele frequencies in a long-format data table
#' and calculates Weir & Cockerhams FST, i.e. theta (Weir & Cockerham, 1984).
#' Permutations can be used to test statistical significance of FST in
#' genotype data sets.
#'
#' @param dat Data table: A long-format data table of biallelic genotypes,
#' coded as '/' separated alleles ('0/0', '0/1', '1/1'), or counts
#' of the Alt alleles (0, 1, 2, respectively). Alternatively, a long-format
#' data table of allele frequencies.
#' Columns required for \strong{both} genotypes and allele frequencies:
#' \enumerate{
#'    \item \code{$POP}, the population ID.
#'    \item \code{$LOCUS}, the locus ID.
#' }
#' Columns required only for \strong{genotypes}:
#' \enumerate{
#'    \item \code{$SAMPLE}, the sample ID (a diploid individual).
#'    \item \code{$GT}, the genotypes.
#' }
#' Columns required only for \strong{allele frequencies}:
#' \enumerate{
#'    \item \code{$FREQ}, the allele frequencies.
#'    \item \code{$INDS}, the number of pooled diploid individuals used to
#'    obtain the allele frequency estimate.
#' }
#'
#' @param type Character: One of \code{'genos'} or \code{'freqs'}, to calculate
#' FST on genotype or allele frequency data, respectively.
#'
#' @param fstat Character: A vector of F-statistics to calculate. This is only
#' applicable for genotype data. Must include one of \code{"FST"},
#' \code{"FIS"}, or \code{"FIT"}.
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
#' @param numPerms Integer: The number of permutations to perform. Default is 100.
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
#' genoGlobal <- fstat_calc(
#'    data_Genos,
#'    type='genos',
#'    global=TRUE,
#'    permute=TRUE,
#'    numPerms=50)
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
#' genoPairs <- fstat_calc(
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
#' genoPairsPerms <- fstat_calc(
#'    data_Genos,
#'    type='genos',
#'    global=FALSE,
#'    pairwise=TRUE,
#'    permute=TRUE,
#'    numPerms=10
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
#' freqPairs <- fstat_calc(
#'    data_PoolFreqs,
#'    type='freqs',
#'    global=FALSE,
#'    pairwise=TRUE,
#'    popCol='POOL'
#' )
#'
#' @export

fstat_calc <- function(
    dat, type, fstat=NULL, global=TRUE, pairwise=FALSE,
    permute=FALSE, numPerms=100, numCores=1
  ){

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
    stop('Argument `global` must be a logical value. See ?fstat_calc.')
  }

  if(class(pairwise)!='logical'){
    stop('Argument `pairwise` must be a logical value. See ?fstat_calc.')
  }

  # Cannot specify bpth global and pairwise
  if(global==TRUE & pairwise==TRUE){
    stop('Arguments `global` and `pairwise` cannot both be TRUE. See ?fstat_calc.')
  }

  # Must specify one of global or pairwise
  if(global==FALSE & pairwise==FALSE){
    stop('One of arguments `global` and `pairwise` must be TRUE. See ?fstat_calc.')
  }

  # Check that permute is logical
  if(class(permute)!='logical'){
    stop('Argument `permute` must be a logical value. See ?fstat_calc.')
  }

  # If permute is specified, then check numPerms is >0 and is an integer
  if(permute==TRUE){
    if(numPerms<1){
      stop('Argument `numPerms` must be an integer >0. See ?fstat_calc.')
    }

    numPerms <- as.integer(numPerms)
  }

  # Checks specific to each type
  if(type=='genos'){
    # Check all necessary columns are present and reassign values.
    if(sum(c("POP", "SAMP", "LOCUS", "GT") %in% colnames(dat))!=4){
      stop(
        'Argument `dat` requires columns "POP", "SAMP", "LOCUS", "GT" to
        estimate F-statistics from genotypes. See ?fstat_calc.')
    }

    # If genotypes are characters, convert to integers.
    if(class(dat$GT)=='character'){
      dat$GT <- genoscore_converter(dat$GT)
    }

    # If no F-statistics have been specified for genotypes, stop.
    if(is.null(fstat)){
      stop(
      'Argument `fstat` cannot be NULL when requesting F-statistics
      from genotypes. See ?fstat.')
    }

    # Make sure the specified value in fstat are correct
    if(sum(!fstat %in% c('FST','FIS','FIT'))!=0){
      stop(
        'Argument `fstat` must only contain "FST", "FIS", and/or "FIT" when
        requesting F-statistics from genotypes. See ?fstat.'
      )
    }
  }

  if(type=='freqs'){
    if(sum(c("POP", "LOCUS", "FREQ", "INDS") %in% colnames(dat))!=4){
      stop(
        'Argument `dat` requires columns "POP", "LOCUS", "FREQ", "INDS" to
        estimate F-statistics from allele frequencies. See ?fstat_calc.')
    }

    if(is.null(fstat)!=TRUE){
      stop(
        'Argument `fstat` must be NULL when requesting FST from allele
        frequencies. See ?fstat_calc.'
      )
    }
  }

  # --------------------------------------------+
  # Code for genotypes
  # --------------------------------------------+

  if(type=='genos'){
    # START: FST among all populations
    if(global==TRUE){
      cat('FST calculation on genotype data, global estimate', '\n')

      # Basic result
      result <- fstat_subfun(dat, type='genos', fstat=fstat)

      # For genotype permutations
      if(permute == TRUE){
        cat('Performing permutations', '\n')
        result$permute <- fstat_permute(
          dat=dat, fstat=fstat, numPerms=numPerms, numCores=numCores
          )
      } else{
        result$permute <- NULL
      }
      # END: FST among all populations
    }

    # START: FST between pairwise populations
    if(pairwise==TRUE){
      cat('FST calculation on genotype data, pairwise estimates', '\n')
      pair.comp <- CJ(POP1=unique(dat$POP), POP2=unique(dat$POP)) %>%
        .[POP1!=POP2]
      num.pairs <- nrow(pair.comp)

      result <- list(genome=list(), locus=list(), permute=list(fstat=NULL, pval=NULL))

      # START: Iterate over i pairs
      for(i in 1:nrow(pair.comp)){
        cat('Performing permutation:', i, '/', num.pairs, '\n')
        pop1 <- pair.comp$POP1[i]
        pop2 <- pair.comp$POP2[i]

        # Calculate genome-wide and locus-wise F-statistics
        pair.calc <- dat[POP %in% c(pop1, pop2)] %>%
          fstat_subfun(., type='genos', fstat=fstat) %>%
          lapply(., function(D){
            data.table(POP1=pop1, POP2=pop2, D)
          })

        # Add genome-wide values into main result
        result$genome[[i]] <- pair.calc$genome
        result$locus[[i]] <- pair.calc$locus

        # If permuting
        if(permute==TRUE){
          # Perform the permutation
          pair.perm <- fstat_permute(
            dat[POP %in% c(pop1, pop2)],
            fstat=fstat, numPerms=numPerms, numCores=numCores
          ) %>%
            data.table(POP1=pop1, POP2=pop2, .)

          # Get the p-value
          pair.pval <- lapply(fstat, function(fx){
            p <- sum(pair.perm[[fx]] > pair.calc$genome[[fx]]) / numPerms
            data.table(STAT=fx, PVAL=p)
          }) %>%
            do.call('rbind', .) %>%
            data.table(POP1=pop1, POP2=pop2, .)

          # Add permutations into main results
          result$permute$fstat[[i]] <- pair.perm
          result$permute$pval[[i]] <- pair.pval
        } else{
          result$permute <- NULL
        }

        rm(i)
        # END: Iterate over i pairs
      }

      # Combine results across indexed pairs into single data table
      result$genome <- do.call('rbind', result$genome)

      result$locus <- do.call('rbind', result$locus)

      if(!is.null(result$permute$fstat)){
        result$permute$fstat <- do.call('rbind', result$permute$fstat)
      }

      if(!is.null(result$permute$pval)){
        result$permute$pval <- do.call('rbind', result$permute$pval) %>%
          setorder(., 'STAT')
      }

      # END: FST between population pairs
    }
  }

  # --------------------------------------------+
  # Code for allele frequencies
  # --------------------------------------------+
  if(type=='freqs'){
    # START: FST among all populations
    if(global==TRUE){
      cat('FST calculation on frequency data, global estimate', '\n')
      result <- fstat_subfun(dat, type='freqs')

      # No pemutations
      result$permute <- NULL

      # END: FST among all populations
    }

    # START: FST between pairwise populations
    if(pairwise==TRUE){
      cat('FST calculation on frequency data, global estimate', '\n')

      pair.comp <- CJ(POP1=unique(dat$POP), POP2=unique(dat$POP)) %>%
        .[POP1!=POP2]
      num.pairs <- nrow(pair.comp)

      result <- list(genome=list(), locus=list())

      # START: Iterate over i population pairs
      for(i in 1:nrow(pair.comp)){
        pop1 <- pair.comp$POP1[i]
        pop2 <- pair.comp$POP2[i]

        pair.calc <- dat[POP %in% c(pop1, pop2)] %>%
          fstat_subfun(., type='freqs') %>%
          lapply(., function(D){
            data.table(POP1=pop1, POP2=pop2, D)
          })

        result$genome[[i]] <- pair.calc$genome
        result$locus[[i]] <- pair.calc$locus
        # END: Iteration over i population pairs
      }

      # Combine results
      result$genome <- do.call('rbind', result$genome)
      result$locus <- do.call('rbind', result$locus)
      result$permute <- NULL
      # END: FST between pairwise populations
    }
  }

  # --------------------------------------------+
  # Output results
  # --------------------------------------------+
  return(result)
}

#' Calculate variance components for F-statistics
#'
#' Takes a list of values of allele frequencies, sample sizes, number of
#' populations, and heterozygosity, and returns variance components, as
#' per Weir & Cockerham (1984). Assumes that all values are for a single
#' biallelic SNP locus. \cr\cr Note, this function is not exported and is an
#' internal function for \code{fstat_calc}.
#'
#' @keywords internal
#'
#' @param pi Numeric: A vector of allele frequencies for each population.
#'
#' @param ni Numeric: A vector of sample sizes (number of diploid individuals)
#' for each population.
#'
#' @param r Integer: A single value, the number of populations.
#'
#' @param hetStand Logical: Should the estimates be standardised for observed
#' heterozygosity?
#'
#' @param hi Numeric: A vector of observed heterozygosities for each population.
#'
#' @returns Returns a data table. If \code{hetStand==FALSE}, then the list has three
#' columns: \code{$MSP}, mean sqaures for populations; \code{$MSG}, mean squares
#' for gametes; \code{$Nc}, a sample size constant. If \code{hetStand==TRUE},
#' then the three columns are \code{$A}, \code{$B}, \code{$C}, which
#' correspond to the 'a', 'b' and 'c' variance components described in
#' Weir & Cockheram (1984).
#'
#' @references Weir & Cockerham (1984) Evolution. DOI: 10.1111/j.1558-5646.1984.tb05657.x
#' Weir et al. (2002) Annals of Human Genetics. DOI: 10.1146/annurev.genet.36
#'
#' @examples
#' fstat_varcomps(pi=c(0.5,0.25), ni=c(20,20), r=2, hetStand=FALSE)
#'
#' fstat_varcomps(pi=c(0.5,0.25), ni=c(20,20), r=2, hi=c(0.05,0.375), hetStand=TRUE)
fstat_varcomps <- function(pi, ni, r, hi=NULL, hetStand=FALSE){
  require(data.table); require(tidyverse)

  if(hetStand==FALSE){
    # Mean weighted allele frequency
    p.mean <- sum(ni * pi)/sum(ni)

    # Mean squares variance components
    msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
    msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

    # Sample size correction factor
    nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

    # Output as data.table
    return(data.table(MSP=msp, MSG=msg, Nc=nc))
  }

  if(hetStand==TRUE){
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

    # Return as list
    return(data.table(A=a, B=b, C=c))
  }
}

#' Subfunction for calculating F-statistics
#'
#' Calculates F-statistics from genotype or allele frequency data in a
#' long-format data table as per Weir & Cockerham (1984). Assumes that all
#' values are for a single biallelic SNP locus. \cr\cr Note, this function is
#' not exported and is an internal function for \code{fstat_calc}.
#'
#' @keywords internal
#'
#' @param dat Data table: In long-format, requires columns \code{$POP},
#' \code{$SAMPLE}, \code{$LOCUS} and \code{$GT} for genotypes;
#' requires \code{$POP}, \code{$LOCUS}, \code{$FREQ}, and \code{$INDS} for
#' allele frequencies
#'
#' @param type Character: either \code{"genos"} or \code{"freqs"}.
#'
#' @param fstat Character: a vector containing the F-statistics to calculate
#' when the data are genotypes. Must include one of \code{"FST"},
#' \code{"FIS"}, or \code{"FIT"}.
#'
#' @returns Returns a list with \code{$genome}, the genome-wide F-statistic, and
#' \code{$locus}, the per locus F-statistic.
#'
#' @examples
#' data_Genos %>%
#'    fstat_subfun(., type='genos', fstat=c('FST','FIS','FIT'))
#'
#' left_join(data_PoolFreqs, data_PoolInfo) %>%
#'    setnames(., 'POOL', 'POP') %>%
#'    fstat_subfun(., type='freqs')
fstat_subfun <- function(dat, type, fstat=NULL){
  require(data.table); require(tidyverse)

  # Empty list
  result <- list(genome=NULL, locus=NULL)

  # If genotype data
  if(type=='genos'){
    # Summarise
    D <- list(
      dat[, .(P=sum(GT)/(length(GT)*2)), by=c('LOCUS','POP')],
      dat[, .(N=length(unique(SAMPLE))), by=c('LOCUS','POP')],
      dat[, .(H=sum(GT==1)/length(GT)), by=c('LOCUS','POP')]
    ) %>%
      reduce(full_join, by=c('LOCUS','POP'))

    num.pops <- D$POP %>% unique %>% length

    # Variance components
    D.varcomp <- D[, fstat_varcomps(pi=P, ni=N, r=num.pops, hi=H, hetStand=TRUE), by=LOCUS]

    # Iterate through requested F-statistics
    for(f in fstat){
      if(f=='FST'){
        result$genome[[f]] <- D.varcomp[, sum(A)/sum(A+B+C)]
        result$locus[[f]] <- D.varcomp[, .(FST=A/(A+B+C)), by=LOCUS]
      } else if(f=='FIS'){
        result$genome[[f]] <- D.varcomp[, 1-(sum(C)/sum(B+C))]
        result$locus[[f]] <- D.varcomp[, .(FIS=1-(C/sum(B+C))), by=LOCUS]
      } else if(f=='FIT'){
        result$genome[[f]] <- D.varcomp[, 1-(sum(C)/sum(A+B+C))]
        result$locus[[f]] <- D.varcomp[, .(FIT=1-(C/sum(A+B+C))), by=LOCUS]
      }
    }

    # Output
    result <- list(
      genome=as.data.table(do.call('cbind', result$genome)),
      locus=result$locus %>% reduce(full_join, by='LOCUS')
    )
  }

  # If allele frequency data
  if(type=='freqs'){
    # Number of populations
    num.pops <- dat$POP %>% unique %>% length

    # Variance components
    D.varcomp <- dat[, fstat_varcomps(pi=FREQ, ni=INDS, r=num.pops, hetStand=FALSE), by=LOCUS]

    # Calculations
    result$genome <- D.varcomp[, .(FST=sum(MSP-MSG)/sum(MSP+((Nc-1)*MSG)))]
    result$locus <- D.varcomp[, .(FST=(MSP-MSG)/(MSP+((Nc-1)*MSG))), by=LOCUS]
  }

  # Output
  return(result)
}

#' Permute the genome-wide F-statistics for genotype data
#'
#' Uses permutation of individuals among populations to generate a null
#' distribution of F-statistics. Assumes biallelic SNP genotypes. \cr\cr
#' Note, this function is not exported and is an internal function for
#' \code{fstat_calc}.
#'
#' @keywords internal
#'
#' @param dat Data table: In long-format, requires columns \code{$POP},
#' \code{$SAMPLE}, \code{$LOCUS} and \code{$GT} for genotypes;
#' requires \code{$POP}, \code{$LOCUS}, \code{$FREQ}, and \code{$INDS} for
#' allele frequencies
#'
#' @param fstat Character: A vector containing the F-statistics to calculate
#' when the data are genotypes. Must include one of \code{"FST"},
#' \code{"FIS"}, or \code{"FIT"}.
#'
#' @param numPerms Integer: The number of permutations to perform.
#'
#' @param numCores Integer: The number of cores to use. Default = 1. When 1,
#' this will be a single core processes using lapply. If >1, then this will
#' utilise the doParallel package to run on a mini cluster.
#'
#' @returns Returns a data table of permuted F-statistics for each permutation.
#'
#' @examples
#' data(data_Genos)
#'
#' fstat_permute(data_Genos, 'FST', 100, 1)
#'
#' fstat_permute(data_Genos, c('FST','FIT','FIS'), 100, 3)
fstat_permute <- function(dat, fstat, numPerms, numCores=1){
  require('tidyverse'); require('doParallel'); require('data.table')

  pop.tab <- dat[, c('POP','SAMPLE')] %>% unique
  num.samps <- pop.tab %>% nrow

  # If single core
  if(numCores==1){
    result <- lapply(1:numPerms, function(i){
      pop.tab.i <- pop.tab %>% copy
      pop.tab.i$SAMPLE <- sample(pop.tab.i$SAMPLE, size=num.samps, replace=FALSE)

      dat[, c('SAMPLE','LOCUS','GT')] %>%
        left_join(., pop.tab.i, by='SAMPLE') %>%
        fstat_subfun(., type='genos', fstat=fstat) %>%
        .[['genome']] %>%
        data.table(PERM=i, .)
    }) %>%
      do.call('rbind',.)
  }

  # If multi core
  if(numCores>1){
    my.cluster <- makeCluster(numCores)
    registerDoParallel(my.cluster)

    result <- foreach(i=1:numPerms) %dopar% {
      require('genomalicious')

      pop.tab.i <- pop.tab %>% copy
      pop.tab.i$SAMPLE <- sample(pop.tab.i$SAMPLE, size=num.samps, replace=FALSE)

      dat[, c('SAMPLE','LOCUS','GT')] %>%
        left_join(., pop.tab.i, by='SAMPLE') %>%
        fstat_subfun(., type='genos', fstat=fstat) %>%
        .[['genome']] %>%
        data.table(PERM=i, .)
    } %>%
      do.call('rbind',.)

    stopCluster(my.cluster)
  }

  return(result)
}


